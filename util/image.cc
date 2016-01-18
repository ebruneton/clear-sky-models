/*
 * Copyright © 2009 Jeff Muizelaar
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of Jeff Muizelaar not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.  Jeff Muizelaar makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 *
 * JEFF MUIZELAAR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT
 * SHALL JEFF MUIZELAAR BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
 * PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 * SOFTWARE.
 */

#include "util/image.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace {

/* PNG specification Annex D */

/* Table of CRCs of all 8-bit messages. */
uint64_t crc_table[256];

/* Flag: has the table been computed? Initially false. */
int crc_table_computed = 0;

/* Make the table for a fast CRC. */
void make_crc_table(void) {
  uint64_t c;
  int n, k;

  for (n = 0; n < 256; n++) {
    c = (uint64_t) n;
    for (k = 0; k < 8; k++) {
      if (c & 1) {
        c = 0xedb88320L ^ (c >> 1);
      } else {
        c = c >> 1;
      }
    }
    crc_table[n] = c;
  }
  crc_table_computed = 1;
}

/* Update a running CRC with the bytes buf[0..len-1]--the CRC
   should be initialized to all 1's, and the transmitted value
   is the 1's complement of the final running CRC (see the
   crc() routine below). */
uint64_t update_crc(uint64_t crc, unsigned char* buf, int len) {
  uint64_t c = crc;
  int n;

  if (!crc_table_computed) {
    make_crc_table();
  }
  for (n = 0; n < len; n++) {
    c = crc_table[(c ^ buf[n]) & 0xff] ^ (c >> 8);
  }
  return c;
}

/* Return the CRC of the bytes buf[0..len-1]. */
uint64_t crc(unsigned char* buf, int len) {
  return update_crc(0xffffffffL, buf, len) ^ 0xffffffffL;
}

unsigned char* write_be32(unsigned char* buf, int a) {
  buf[0] = (a >> 24) & 0xff;
  buf[1] = (a >> 16) & 0xff;
  buf[2] = (a >> 8)  & 0xff;
  buf[3] = (a >> 0)  & 0xff;
  return &buf[4];
}

struct buf {
  unsigned char *data;
  int len;
};

struct buf chunk(const char* type, struct buf q) {
  struct buf b;
  unsigned char* t = (unsigned char*) malloc(4 + 4 + q.len + 4);
  unsigned char* chunk_start;
  b.data = t;
  t = write_be32(t, q.len);
  chunk_start = t;

  /* copy over type */
  t[0] = type[0];
  t[1] = type[1];
  t[2] = type[2];
  t[3] = type[3];
  t = &t[4];

  /* copy over data */
  memcpy(t, q.data, q.len);
  t += q.len;
  free(q.data);

  t = write_be32(t, crc(chunk_start, q.len + 4));

  b.len = 4 + 4 + q.len + 4;
  return b;
}

struct buf buf_cat_str(struct buf b, unsigned char* d, int len) {
  struct buf r;
  r.data = (unsigned char*) realloc(b.data, b.len + len);
  memcpy(r.data + b.len, d, len);
  r.len = b.len + len;
  return r;
}

typedef struct buf (*data_cat_fn)(struct buf b, void* src, int len);

struct buf buf_cat_str_argb(struct buf b, void* abstract_src, int len) {
  unsigned int* src = (unsigned int*) abstract_src;
  struct buf r;
  unsigned char* dest;
  r.data = (unsigned char*) realloc(b.data, b.len + len);
  dest = r.data + b.len;
  r.len = b.len + len;
  while (len) {
    /* we probably need to unpremultiply here */
    unsigned char alpha = (*src >> 24) & 0xff;
    unsigned char red   = (*src >> 16) & 0xff;
    unsigned char green = (*src >> 8) & 0xff;
    unsigned char blue  = (*src >> 0) & 0xff;
    *dest++ = red;
    *dest++ = green;
    *dest++ = blue;
    *dest++ = alpha;
    len -= 4;
    src++;
  }
  return r;
}

struct buf buf_cat(struct buf b, struct buf c) {
  return buf_cat_str(b, c.data, c.len);
}

struct buf be32(int a) {
  struct buf r;
  r.data = (unsigned char*) malloc(4);
  r.len = 4;
  write_be32(r.data, a);
  return r;
}

struct adler {
  int a;  // 1
  int b;  // 0
};

#define MOD_ADLER 65521

/* data: Pointer to the data to be summed; len is in bytes */
unsigned int adler32(unsigned char* data, size_t len) {
  unsigned int a = 1, b = 0;

  while (len != 0) {
    a = (a + *data++) % MOD_ADLER;
    b = (b + a) % MOD_ADLER;
    len--;
  }

  return (b << 16) | a;
}

/* From Wikipedia */
/* data: Pointer to the data to be summed; len is in bytes */
struct adler adler32_str(struct adler ctx, unsigned char* data, size_t len) {
  while (len != 0) {
    ctx.a = (ctx.a + *data++) % MOD_ADLER;
    ctx.b = (ctx.b + ctx.a) % MOD_ADLER;
    len--;
  }

  return ctx;
}

/* data: Pointer to the data to be summed; len is in bytes */
struct adler adler32_buf(struct adler ctx, struct buf b)  {
  int i = 0;
  while (i != b.len) {
    ctx.a = (ctx.a + b.data[i++]) % MOD_ADLER;
    ctx.b = (ctx.b + ctx.a) % MOD_ADLER;
  }

  return ctx;
}

struct adler adler32_init() {
  struct adler ret;
  ret.a = 1;
  ret.b = 0;
  return ret;
}

unsigned int adler32_fin(struct adler ctx) {
  return (ctx.b << 16) | ctx.a;
}

/* 16bit length in little endian followed by
 * ones compliment of length in little endian */
struct buf zlib_block_length(int length) {
  struct buf r;
  r.data = (unsigned char *)malloc(4);
  r.len = 4;
  r.data[0] = length & 0xff;
  r.data[1] = (length >> 8) & 0xff;
  length = ~length;
  r.data[2] = length & 0xff;
  r.data[3] = (length >> 8) & 0xff;
  return r;
}

/* inspired by "A use for uncompressed PNGs"
 * http://drj11.wordpress.com/2007/11/20/a-use-for-uncompressed-pngs/
 * by David Jones */
struct buf make_png(void* d, int width, int height, int stride,
    data_cat_fn buf_cat_str_data) {
  struct buf r = {0};
  unsigned char predictor[] = {0x0};
  unsigned char zlib_prefix[] = {0x78, 0x9c};
  unsigned char zlib_final_block_prefix[] = {0x01};
  unsigned char zlib_block_prefix[] = {0x00};
  unsigned char hdr_tail[] = {0x08, 0x06, 0x00, 0x00, 0x00};
  unsigned char png_start[] = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
  struct buf ihdr = {0};
  int block_length = (width * 4 + 1);
  struct buf idat = {0}, iend = {0};
  int i;
  assert(block_length <= 65535);

  r = buf_cat_str(r, png_start, sizeof(png_start));

  ihdr = buf_cat(ihdr, be32(width));
  ihdr = buf_cat(ihdr, be32(height));
  ihdr = buf_cat_str(ihdr, hdr_tail, sizeof(hdr_tail));

  r = buf_cat(r, chunk("IHDR", ihdr));

  idat = buf_cat_str(idat, zlib_prefix, sizeof(zlib_prefix));

  struct adler chksum = adler32_init();
  for (i = 0 ; i < height; i++) {
    struct buf row_data = {0};
    if (i == height - 1) {
      idat = buf_cat_str(
          idat, zlib_final_block_prefix, sizeof(zlib_final_block_prefix));
    } else {
      idat = buf_cat_str(idat, zlib_block_prefix, sizeof(zlib_block_prefix));
    }
    idat = buf_cat(idat, zlib_block_length(block_length));
    chksum = adler32_str(chksum, predictor, sizeof(predictor));
    idat = buf_cat_str(idat, predictor, sizeof(predictor));

    row_data = buf_cat_str_data(row_data, d, width * 4);
    chksum = adler32_buf(chksum, row_data);
    idat = buf_cat_str_data(idat, d, width * 4);
    free(row_data.data);
    d = (unsigned char*) d + stride;
  }
  // be32 leaks
  idat = buf_cat(idat, be32(adler32_fin(chksum)));

  r = buf_cat(r, chunk("IDAT", idat));
  r = buf_cat(r, chunk("IEND", iend));
  return r;
}

}  // anonymous namespace

void WritePngArgb(const std::string& name, void* pixels, int width,
    int height) {
  FILE *f = fopen(name.c_str(), "wb+");
  struct buf png = make_png(pixels, width, height, 4 * width, buf_cat_str_argb);
  fwrite(png.data, png.len, 1, f);
  free(png.data);
  fclose(f);
}
