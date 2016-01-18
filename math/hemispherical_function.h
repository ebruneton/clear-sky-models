/**
 Copyright (c) 2015 Eric Bruneton
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef MATH_HEMISPHERICAL_FUNCTION_H_
#define MATH_HEMISPHERICAL_FUNCTION_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "math/angle.h"

// A function from the hemisphere to values of type T represented by its values
// at 9x9 sample points mapped from the unit square to the hemisphere (using
// "A low distortion map between disk and square.", Hirley et al. 1997), and
// interpolated with bicubic interpolation in between.
template<class T>
class HemisphericalFunction {
 public:
  inline const T& Get(int i, int j) const {
    assert(i >= 0 && i < 9 && j >= 0 && j < 9);
    return value_[i][j];
  }

  inline void Set(int i, int j, const T& value) {
    assert(i >= 0 && i < 9 && j >= 0 && j < 9);
    value_[i][j] = value;
  }

  static void GetSampleDirection(int i, int j, Angle* view_zenith,
      Angle* view_azimuth) {
    double r;
    double phi;
    ToUnitDisk((i + 0.5) / 9, (j + 0.5) / 9, &r, &phi);
    *view_zenith = acos(1.0 - r * r) * rad;
    *view_azimuth = phi *rad;
  }

  T operator()(Angle view_zenith, Angle view_azimuth) const {
    // Computes the nearest ring indices corresponding to view_zenith and
    // returns the cubic interpolation of the function values on these rings
    // (themselves interpolated from the sample values using cubic
    // interpolation, see below).
    Number z = cos(view_zenith);
    double r = std::sqrt(1.0 - z());
    double u = 4.5 * r;
    int i = std::floor(u);
    u -= i;
    return Cubic(u, GetRing(i - 1, view_azimuth), GetRing(i, view_azimuth),
        GetRing(i + 1, view_azimuth), GetRing(i + 2, view_azimuth));
  }

  void Load(const std::string& filename) {
    std::ifstream file(filename, std::ifstream::binary | std::ifstream::in);
    for (int i = 0; i < 9; ++i) {
      file.read(reinterpret_cast<char*>(&(value_[i])), 9 * sizeof(T));
    }
    file.close();
  }

  void Save(const std::string& filename) const {
    std::ofstream file(filename, std::ofstream::binary);
    for (int i = 0; i < 9; ++i) {
      file.write(reinterpret_cast<const char*>(&(value_[i])), 9 * sizeof(T));
    }
    file.close();
  }

 private:
  // From "A low distortion map between disk and square.", Hirley et al. 1997.
  static void ToUnitDisk(double x, double y, double *r, double *phi) {
    double a = 2 * x - 1;
    double b = 2 * y - 1;
    if (a > -b) {
      if (a > b) {
        *r = a;
        *phi = (PI / 4) * (b / a);
      } else {
        *r = b;
        *phi = (PI / 4) * (2 - (a / b));
      }
    } else {
      if (a < b) {
        *r = -a;
        *phi = (PI / 4) * (4 + (b / a));
      } else {
        *r = -b;
        if (b != 0.0) {
          *phi = (PI / 4) * (6 - (a / b));
        } else {
          *phi = 0;
        }
      }
    }
  }

  // Returns the interpolated value of the hemispherical function on the given
  // ring (see GetRing below), using cubic interpolation of the nearest samples
  // on this ring. Also extrapolate the results to ring index -1 such that the
  // interpolated function is differentiable at the zenith, and to ring indices
  // larger than 4 to extrapolate the function up to the horizon.
  T GetRing(int ring_index, Angle view_azimuth) const {
    assert(ring_index >= -1);
    if (ring_index == -1) {
      // Computes the value such that the finite difference derivative
      //   0.5 * (GetRing(1, view_azimuth) - GetRing(-1, view_azimuth))
      // is equal to GetZenithDerivative(view_azimuth).
      return GetRing(1, view_azimuth) - GetZenithDerivative(view_azimuth) * 2.0;
    } else if (ring_index > 4) {
      // Noting R(i) = GetRing(ring_index + i, view_azimuth), computes the value
      // R(0) such that the second derivative at -1, R(-2) - 2R(-1) + R(0) is
      // equal to the second derivative at -2, R(-3) - 2R(-2) + R(-1).
      return GetRing(ring_index - 3, view_azimuth) -
          GetRing(ring_index - 2, view_azimuth) * 3.0 +
          GetRing(ring_index - 1, view_azimuth) * 3.0;
    } else {
      // Normal case: cubic interpolation of the nearest samples along the ring.
      double u = view_azimuth.to(pi) * (4 * ring_index);
      int i = std::floor(u);
      u -= i;
      return Cubic(u, GetRing(ring_index, i - 1), GetRing(ring_index, i),
          GetRing(ring_index, i + 1), GetRing(ring_index, i + 2));
    }
  }

  // Returns the value of the sample point designated by the given ring index
  // (0 is for the zenith sample 4,4, 1 is for samples 5,4 5,5 4,5 3,5 3,4 etc,
  // and so on up to ring index 4) and index in ring (0 is for azimuth 0,
  // ring_index for azimuth 45, k * ring_index for azimuth k * 45).
  const T& GetRing(int ring_index, int index_in_ring) const {
    assert(ring_index >= 0 && ring_index < 5);
    if (ring_index == 0) {
      return value_[4][4];
    } else {
      while (index_in_ring < 0) {
        index_in_ring += 8 * ring_index;
      }
      int dx = std::max(-ring_index, std::min(ring_index,
          std::abs(index_in_ring % (8 * ring_index) - 4 * ring_index) -
              2 * ring_index));
      int dy = std::max(-ring_index, std::min(ring_index,
          std::abs((index_in_ring + 6 * ring_index) % (8 * ring_index) -
              4 * ring_index) - 2 * ring_index));
      return value_[4 + dx][4 + dy];
    }
  }

  // Returns the derivative of the interpolated function at the zenith, in
  // direction view_azimuth. For this, estimates the derivatives in the x and y
  // directions from the 8 samples around the zenith, and returns a linear
  // combination of those depending on the azimuth angle.
  T GetZenithDerivative(Angle view_azimuth) const {
    constexpr double norm = 1.0 / sqrt(2.0);
    constexpr double scale = 0.5 / (1.0 + 2.0 * norm);
    T derivative_x = ((value_[5][4] - value_[3][4]) +
        (value_[5][5] - value_[3][3]) * norm +
        (value_[5][3] - value_[3][5]) * norm) * scale;
    T derivative_y = ((value_[4][5] - value_[4][3]) +
        (value_[5][5] - value_[3][3]) * norm +
        (value_[3][5] - value_[5][3]) * norm) * scale;
    return
        derivative_x * cos(view_azimuth)() + derivative_y * sin(view_azimuth)();
  }

  static T Cubic(double u, T am1, T a0, T a1, T a2) {
    return (a0 * 2.0 + ((-am1 + a1) + ((am1 * 2.0 - a0 * 5.0 + a1 * 4.0 - a2) +
        (-am1 + a0 * 3.0 - a1 * 3.0 + a2) * u) * u) * u) * 0.5;
  }

  // The function values in the unit square, mapped to the hemisphere using
  // "A low distortion map between disk and square.", Hirley et al. 1997. The
  // positive x axis corresponds to azimuth 0 (North), and the positive y axis
  // to azimuth 90 (East). More precisely, element [4][4] corresponds to the
  // zenith, element [8][4] corresponds to azimuth 0 and elevation 12.1151 from
  // the horizon, and element[4][8] to azimuth 90 and the same elevation.
  T value_[9][9];
};

#endif  // MATH_HEMISPHERICAL_FUNCTION_H_
