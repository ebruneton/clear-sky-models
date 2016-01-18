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
#include "atmosphere/comparisons.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "math/hemispherical_function.h"
#include "physics/cie.h"
#include "util/image.h"

namespace {

const std::string kOutputDirectory = "output/comparisons/";

const int kColorMap[162] = {
  109,  90, 203, 126, 124, 217, 138, 148, 226, 148, 169, 232, 157, 187, 238,
  165, 203, 242, 173, 218, 246, 179, 231, 250, 186, 244, 253,   0, 102, 181,
    0, 131, 201,   0, 151, 213,   0, 167, 223,   0, 180, 230,   0, 192, 237,
    0, 203, 242,   0, 213, 247,   0, 222, 251,  72, 185, 171,  90, 201, 189,
  105, 213, 203, 119, 222, 215, 132, 231, 225, 143, 238, 234, 154, 244, 242,
  165, 250, 249, 174, 255, 255, 136, 196, 145, 151, 209, 158, 162, 219, 168,
  172, 227, 176, 180, 234, 183, 187, 241, 189, 193, 246, 195, 199, 251, 200,
  204, 255, 204, 196, 167, 125, 209, 185, 143, 219, 199, 159, 227, 211, 173,
  234, 222, 186, 240, 231, 198, 246, 240, 209, 251, 248, 220, 255, 255, 230,
  212, 112,  84, 222, 134, 105, 229, 154, 123, 235, 172, 139, 240, 190, 154,
  245, 207, 168, 249, 224, 181, 252, 240, 193, 255, 255, 204
};

const int kErrorColorMap[123] = {
    8,  54, 106,  15,  67, 123,  23,  82, 144,  29,  95, 162,  36, 106, 174,
  +46, 119, 181,  54, 129, 186,  63, 142, 192,  76, 153, 198,  95, 165, 205,
  117, 178, 212, 135, 190, 218, 155, 201, 224, 169, 209, 229, 184, 216, 233,
  202, 225, 238, 213, 231, 241, 224, 236, 243, 233, 240, 244, 242, 245, 246,
  248, 243, 240, 249, 237, 229, 251, 229, 216, 252, 222, 205, 252, 213, 191,
  249, 198, 172, 247, 185, 156, 245, 170, 137, 240, 156, 123, 233, 139, 110,
  225, 120,  96, 218, 104,  83, 208,  85,  72, 200,  68,  64, 191,  51,  56,
  182,  31,  46, 168,  21,  41, 147,  14,  38, 129,   8,  35, 112,   3,  32,
  103,   0,  31
};

void GetColor(double x, int* red, int* green, int* blue) {
  int n = floor(log(x) / log(10.0));
  int i = floor(x / pow(10.0, n)) - 1;
  int color_index = 3 * (i + 9 * n);
  if (color_index < 0) {
    *red = 255;
    *green = 0;
    *blue = 0;
  } else if (color_index > 159) {
    *red = 255;
    *green = 255;
    *blue = 0;
  } else {
    *red = kColorMap[color_index];
    *green = kColorMap[color_index + 1];
    *blue = kColorMap[color_index + 2];
  }
}

void GetErrorColor(double x, int* red, int* green, int* blue) {
  int n = floor(x * 10.0) + 20;
  int color_index = 3 * std::max(0, std::min(40, n));
  *red = kErrorColorMap[color_index];
  *green = kErrorColorMap[color_index + 1];
  *blue = kErrorColorMap[color_index + 2];
}

uint32_t Blend(uint32_t background, uint32_t foreground) {
  int red1 = (background >> 16) & 0xFF;
  int green1 = (background >> 8) & 0xFF;
  int blue1 = background & 0xFF;
  int red2 = (foreground >> 16) & 0xFF;
  int green2 = (foreground >> 8) & 0xFF;
  int blue2 = foreground & 0xFF;
  double blend = (foreground >> 24) / 255.0;
  int red = red1 * (1.0 - blend) + red2 * blend;
  int green = green1 * (1.0 - blend) + green2 * blend;
  int blue = blue1 * (1.0 - blend) + blue2 * blend;
  return (0xFF << 24) | (red << 16) | (green << 8) | blue;
}

void DrawCross(Angle sun_zenith, Angle sun_azimuth, int size, uint32_t color,
    int width, int height, uint32_t* pixels) {
  Number radius = sun_zenith / (pi / 2.0);
  Number x = radius * sin(sun_azimuth);
  Number y = radius * cos(sun_azimuth);
  int i = static_cast<int>((x * (width / 2.0) + width / 2.0 - 0.5)());
  int j = static_cast<int>((y * (height / 2.0) + height / 2.0 - 0.5)());
  for (int x = -size; x <= size; ++x) {
    pixels[i + x + j * width] = Blend(pixels[i + x + j * width], color);
  }
  for (int x = -size; x <= size; ++x) {
    if (x != 0) {
      pixels[i + (j + x) * width] = Blend(pixels[i + (j + x) * width], color);
    }
  }
}

void DrawCrosses(Angle sun_zenith, Angle sun_azimuth, int width, int height,
    uint32_t* pixels) {
  DrawCross(sun_zenith, sun_azimuth, 1, 0xFFFF0000, width, height, pixels);
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      Angle view_zenith;
      Angle view_azimuth;
      HemisphericalFunction<Number>::GetSampleDirection(
          i, j, &view_zenith, &view_azimuth);
      DrawCross(
          view_zenith, view_azimuth, 0, 0x40000000, width, height, pixels);
    }
  }
}

}  // anonymous namespace

Comparisons::Comparisons(const std::string& name, const Atmosphere& atmosphere,
    const MeasuredAtmospheres& reference, Wavelength min_wavelength,
    Wavelength max_wavelength)
    : name_(name), atmosphere_(atmosphere), reference_(reference),
      min_wavelength_(min_wavelength), max_wavelength_(max_wavelength) {}

void Comparisons::PlotSunIlluminanceAttenuation() const {
  const DimensionlessSpectrum& y_bar = cie_y_bar_function();
  Irradiance extra_terrestrial_illuminance = Integral(SolarSpectrum() * y_bar);

  std::ofstream file(kOutputDirectory + "sun_illuminance_attenuation_" +
      name_ + ".txt");
  for (int i = 0; i <= 90; ++i) {
    Angle sun_zenith = i * deg;
    Irradiance ground_illuminance =
        Integral(atmosphere_.GetSunIrradiance(0.0 * m, sun_zenith) * y_bar);
    file << sun_zenith.to(deg) << " "
         << (ground_illuminance / extra_terrestrial_illuminance)() << std::endl;
  }
  file.close();
}

void Comparisons::PlotRadiance(const std::string& name, Angle sun_zenith,
    Angle sun_azimuth, Angle view_zenith, Angle view_azimuth) const {
  RadianceSpectrum radiance = atmosphere_.GetSkyRadiance(
      0.0 * m, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
  std::stringstream filename;
  filename << kOutputDirectory << "radiance_" << name << "_"
           << view_zenith.to(deg) << "_" << view_azimuth.to(deg) << "_" << name_
           << ".txt";
  std::ofstream file(filename.str());
  for (unsigned int i = 0; i < radiance.size(); ++i) {
    file << radiance.GetWavelength(i).to(nm) << " "
         << radiance[i].to(watt_per_square_meter_per_sr_per_nm) << std::endl;
  }
  file.close();
}

void Comparisons::RenderLuminanceAndImage(const std::string& name,
    Angle sun_zenith, Angle sun_azimuth) const {
  const auto zenith_luminance = 683.0 * Integral(cie_y_bar_function() *
      atmosphere_.GetSkyRadiance(
          0.0 * m, sun_zenith, sun_azimuth, 0.0 * deg, 0.0 * deg));

  // Compute the normalization constants k_r,k_g,k_b for non spectral rendering.
  WavelengthFunction<Scalar<0, -3, 0, 0>> f;
  for (unsigned int i = 0; i < f.size(); ++i) {
    Wavelength lambda = f.GetWavelength(i);
    f[i] = 1.0 / (lambda * lambda * lambda);
  }
  const Wavelength lambda_r = 680.0 * nm;
  const Wavelength lambda_g = 550.0 * nm;
  const Wavelength lambda_b = 440.0 * nm;
  const auto X0 = Integral(cie_x_bar_function() * SolarSpectrum() * f);
  const auto Y0 = Integral(cie_y_bar_function() * SolarSpectrum() * f);
  const auto Z0 = Integral(cie_z_bar_function() * SolarSpectrum() * f);
  const auto k_r = (3.2404542 * X0 - 1.5371385 * Y0 - 0.4985314 * Z0) /
      SolarSpectrum()(lambda_r) * (lambda_r * lambda_r * lambda_r);
  const auto k_g = (-0.9692660 * X0 + 1.8760108 * Y0 + 0.0415560 * Z0) /
      SolarSpectrum()(lambda_g) * (lambda_g * lambda_g * lambda_g);
  const auto k_b = (0.0556434 * X0 - 0.2040259 * Y0 + 1.0572252 * Z0) /
      SolarSpectrum()(lambda_b) * (lambda_b * lambda_b * lambda_b);

  const int width = 256;
  const int height = 256;
  std::unique_ptr<uint32_t[]> pixels1(new uint32_t[width * height]);
  std::unique_ptr<uint32_t[]> pixels2(new uint32_t[width * height]);
  std::unique_ptr<uint32_t[]> pixels3(new uint32_t[width * height]);
  std::unique_ptr<uint32_t[]> pixels4(new uint32_t[width * height]);
  std::unique_ptr<uint32_t[]> pixels5(new uint32_t[width * height]);
  std::unique_ptr<uint32_t[]> pixels6(new uint32_t[width * height]);

  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      Number x = (i + 0.5 - width / 2.0) / (width / 2.0);
      Number y = (height / 2.0 - 0.5 - j) / (height / 2.0);
      Number radius = sqrt(x * x + y * y);
      if (radius >= 1.0) {
        pixels1[i + j * width] = 0;
        pixels2[i + j * width] = 0;
        pixels3[i + j * width] = 0;
        pixels4[i + j * width] = 0;
        pixels5[i + j * width] = 0;
        pixels6[i + j * width] = 0;
        continue;
      }
      Angle view_zenith = radius * pi / 2.0;
      Angle view_azimuth = atan2(x, -y);
      RadianceSpectrum radiance = atmosphere_.GetSkyRadiance(
          0.0 * m, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
      const auto view_dir_luminance =
          683.0 * Integral(cie_y_bar_function() * radiance);
      Number relative_luminance = view_dir_luminance / zenith_luminance;

      int red, green, blue;
      GetColor(view_dir_luminance.to(watt_per_square_meter_per_sr),
          &red, &green, &blue);
      pixels1[i + j * width] = (0xFF << 24) | (red << 16) | (green << 8) | blue;

      GetColor(relative_luminance() * 100.0, &red, &green, &blue);
      pixels2[i + j * width] = (0xFF << 24) | (red << 16) | (green << 8) | blue;

      const Radiance X = Integral(cie_x_bar_function() * radiance);
      const Radiance Y = Integral(cie_y_bar_function() * radiance);
      const Radiance Z = Integral(cie_z_bar_function() * radiance);
      // XYZ to sRGB matrix.
      const Radiance R = std::max(
          3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z,
          0.0 * watt_per_square_meter_per_sr);
      const Radiance G = std::max(
          -0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z,
          0.0 * watt_per_square_meter_per_sr);
      const Radiance B = std::max(
          0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z,
          0.0 * watt_per_square_meter_per_sr);

      const Radiance RGB_max = std::max(R, std::max(G, B));
      red = 255.0 * (R / RGB_max)();
      green = 255.0 * (G / RGB_max)();
      blue = 255.0 * (B / RGB_max)();
      pixels3[i + j * width] = (0xFF << 24) | (red << 16) | (green << 8) | blue;

      constexpr Radiance c = 5.0 * watt_per_square_meter_per_sr;
      red = 255.0 * (1.0 - exp(-R / c))();
      green = 255.0 * (1.0 - exp(-G / c))();
      blue = 255.0 * (1.0 - exp(-B / c))();
      pixels4[i + j * width] = (0xFF << 24) | (red << 16) | (green << 8) | blue;

      int r = 255.0 * (1.0 - exp(-k_r * radiance(lambda_r) / c))();
      int g = 255.0 * (1.0 - exp(-k_g * radiance(lambda_g) / c))();
      int b = 255.0 * (1.0 - exp(-k_b * radiance(lambda_b) / c))();
      pixels5[i + j * width] = (0xFF << 24) | (r << 16) | (g << 8) | b;

      r = std::abs(r - red) * 10;
      g = std::abs(g - green) * 10;
      b = std::abs(b - blue) * 10;
      pixels6[i + j * width] = (0xFF << 24) | (r << 16) | (g << 8) | b;
    }
  }
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels1.get());
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels2.get());
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels3.get());
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels4.get());
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels5.get());
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels6.get());

  std::string filename1 =
      "output/figures/absolute_luminance_" + name + "_" + name_ + ".png";
  WritePngArgb(filename1, pixels1.get(), width, height);

  std::string filename2 =
      "output/figures/relative_luminance_" + name + "_" + name_ + ".png";
  WritePngArgb(filename2, pixels2.get(), width, height);

  std::string filename3 =
      "output/figures/chromaticity_" + name + "_" + name_ + ".png";
  WritePngArgb(filename3, pixels3.get(), width, height);

  std::string filename4 =
      "output/figures/image_" + name + "_" + name_ + ".png";
  WritePngArgb(filename4, pixels4.get(), width, height);

  std::string filename5 =
      "output/figures/image_rgb_" + name + "_" + name_ + ".png";
  WritePngArgb(filename5, pixels5.get(), width, height);

  std::string filename6 =
      "output/figures/image_rgb_diff_" + name + "_" + name_ + ".png";
  WritePngArgb(filename6, pixels6.get(), width, height);
}

void Comparisons::PlotLuminanceProfile(const std::string& name,
    Angle sun_zenith, Angle sun_azimuth) const {
  std::string filename =
      "output/comparisons/luminance_profile_" + name + "_" + name_ + ".txt";
  std::ofstream file(filename);
  if (&atmosphere_ == &reference_) {
    for (int i = 0; i < 9; ++i) {
      Angle view_zenith;
      Angle view_azimuth;
      HemisphericalFunction<Number>::GetSampleDirection(
          8 - i, i, &view_zenith, &view_azimuth);
      Radiance measured = 683.0 * Integral(cie_y_bar_function() *
          reference_.GetSkyRadianceMeasurement(0.0 * m, sun_zenith, sun_azimuth,
              view_zenith, view_azimuth));
      file << (i < 4 ? -view_zenith.to(deg) : view_zenith.to(deg))
           << " " << measured.to(watt_per_square_meter_per_sr) << std::endl;
    }
  } else {
    for (int i = 2; i <= 62; ++i) {
      Angle view_zenith = std::abs(i - 32) / 32.0 * pi / 2.0;
      Angle view_azimuth = i < 32 ? -45 * deg : 135 * deg;
      Radiance model = 683.0 * Integral(cie_y_bar_function() *
          atmosphere_.GetSkyRadiance(0.0 * m, sun_zenith, sun_azimuth,
              view_zenith, view_azimuth));
      file << (i < 32 ? -view_zenith.to(deg) : view_zenith.to(deg))
           << " " << model.to(watt_per_square_meter_per_sr) << std::endl;
    }
  }
  file.close();
}

void Comparisons::PlotRelativeError(const std::string& name, Angle sun_zenith,
    Angle sun_azimuth) const {
  // Relative error for the measured directions, and interpolated in between.
  HemisphericalFunction<Number> error_function;
  int count = 0;
  SpectralRadiance sum = 0.0 * watt_per_square_meter_per_sr_per_nm;
  auto sum_square_error = sum * sum;
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      Angle view_zenith;
      Angle view_azimuth;
      HemisphericalFunction<Number>::GetSampleDirection(
          i, j, &view_zenith, &view_azimuth);
      RadianceSpectrum model_spectrum = atmosphere_.GetSkyRadiance(
          0.0 * m, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
      RadianceSpectrum measured_spectrum = reference_.GetSkyRadianceMeasurement(
          0.0 * m, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
      Radiance model =
          Integral(model_spectrum, min_wavelength_, max_wavelength_);
      Radiance measured =
          Integral(measured_spectrum, min_wavelength_, max_wavelength_);
      Number relative_error = (model - measured) / measured;
      for (unsigned int k = 0; k < model_spectrum.size(); ++k) {
        if (model_spectrum.GetWavelength(k) >= min_wavelength_ &&
            model_spectrum.GetWavelength(k) <= max_wavelength_) {
          count += 1;
          sum += model_spectrum[k];
          sum_square_error += (model_spectrum[k] - measured_spectrum[k]) *
              (model_spectrum[k] - measured_spectrum[k]);
        }
      }
      error_function.Set(i, j, relative_error());
    }
  }
  double avg = (sum / count).to(watt_per_square_meter_per_sr_per_nm);
  double rmse = (sqrt(sum_square_error / count)).to(
      watt_per_square_meter_per_sr_per_nm);
  double cvrmse = static_cast<int>(round((rmse / avg) * 100.0));
  std::ofstream file(kOutputDirectory + "error_" + name + "_" + name_ + ".txt");
  file << avg << std::endl << rmse << std::endl;
  file.close();
  std::ofstream file_percent(
      "paper/figures/error_percent_" + name + "_" + name_ + ".txt");
  file_percent << cvrmse << std::endl;
  file_percent.close();

  const int width = 256;
  const int height = 256;
  std::unique_ptr<uint32_t[]> pixels(new uint32_t[width * height]);
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      Number x = (i + 0.5 - width / 2.0) / (width / 2.0);
      Number y = (height / 2.0 - 0.5 - j) / (height / 2.0);
      Number radius = sqrt(x * x + y * y);
      if (radius >= 1.0) {
        pixels[i + j * width] = 0;
        continue;
      }
      Angle view_zenith = radius * pi / 2.0;
      Angle view_azimuth = atan2(x, -y);
      double error = error_function(view_zenith, view_azimuth)();
      int red, green, blue;
      GetErrorColor(error, &red, &green, &blue);
      pixels[i + j * width] = (0xFF << 24) | (red << 16) | (green << 8) | blue;
    }
  }
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels.get());
  std::string filename =
      "output/figures/relative_error_" + name + "_" + name_ + ".png";
  WritePngArgb(filename, pixels.get(), width, height);
}

void Comparisons::PlotRelativeError(const std::string& name,
    const Atmosphere& reference, Angle sun_zenith, Angle sun_azimuth) const {
  const int width = 256;
  const int height = 256;
  std::unique_ptr<uint32_t[]> pixels(new uint32_t[width * height]);
  int count = 0;
  SpectralRadiance sum = 0.0 * watt_per_square_meter_per_sr_per_nm;
  auto sum_square_error = sum * sum;
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      Number x = (i + 0.5 - width / 2.0) / (width / 2.0);
      Number y = (height / 2.0 - 0.5 - j) / (height / 2.0);
      Number radius = sqrt(x * x + y * y);
      if (radius >= 1.0) {
        pixels[i + j * width] = 0;
        continue;
      }
      Angle view_zenith = radius * pi / 2.0;
      Angle view_azimuth = atan2(x, -y);
      RadianceSpectrum model_spectrum = atmosphere_.GetSkyRadiance(
          0.0 * m, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
      RadianceSpectrum reference_spectrum = reference.GetSkyRadiance(
          0.0 * m, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
      Radiance model =
          Integral(model_spectrum, min_wavelength_, max_wavelength_);
      Radiance ref =
          Integral(reference_spectrum, min_wavelength_, max_wavelength_);
      Number relative_error = (model - ref) / ref;
      for (unsigned int k = 0; k < model_spectrum.size(); ++k) {
        if (model_spectrum.GetWavelength(k) >= min_wavelength_ &&
            model_spectrum.GetWavelength(k) <= max_wavelength_) {
          count += 1;
          sum += model_spectrum[k];
          sum_square_error += (model_spectrum[k] - reference_spectrum[k]) *
              (model_spectrum[k] - reference_spectrum[k]);
        }
      }
      int red, green, blue;
      GetErrorColor(relative_error(), &red, &green, &blue);
      pixels[i + j * width] = (0xFF << 24) | (red << 16) | (green << 8) | blue;
    }
  }
  DrawCrosses(sun_zenith, sun_azimuth, width, height, pixels.get());
  std::string filename =
      "output/figures/relative_error_" + name + "_" + name_ + ".png";
  WritePngArgb(filename, pixels.get(), width, height);

  double avg = (sum / count).to(watt_per_square_meter_per_sr_per_nm);
  double rmse = (sqrt(sum_square_error / count)).to(
      watt_per_square_meter_per_sr_per_nm);
  double cvrmse = static_cast<int>(round((rmse / avg) * 100.0));
  std::ofstream file(kOutputDirectory + "error_" + name + "_" + name_ + ".txt");
  file << avg << std::endl << rmse << std::endl;
  file.close();
  std::ofstream file_percent(
      "paper/figures/error_percent_" + name + "_" + name_ + ".txt");
  file_percent << cvrmse << std::endl;
  file_percent.close();
}

Radiance Comparisons::ComputeZenithLuminance(Angle sun_zenith) const {
  return 683.0 * Integral(cie_y_bar_function() *
      atmosphere_.GetSkyRadiance(0.0 * m, sun_zenith, 0.0 * deg, 0.0 * deg));
}

void Comparisons::ComputeIrradiance(Angle sun_zenith, Irradiance* sun,
      Irradiance* sky) const {
  Number cosine = cos(sun_zenith);
  *sky = Integral(atmosphere_.GetSkyIrradiance(0.0 * m, sun_zenith),
      min_wavelength_, max_wavelength_);
  *sun = Integral(atmosphere_.GetSunIrradiance(0.0 * m, sun_zenith) * cosine,
      min_wavelength_, max_wavelength_);
}

void Comparisons::PlotDayZenithLuminance(
    const std::vector<Angle>& sun_zenith,
    const std::vector<Radiance>& zenith_luminance) const {
  std::ofstream file(
      kOutputDirectory + "day_zenith_luminance_" + name_ + ".txt");
  for (unsigned int i = 0; i < sun_zenith.size(); ++i) {
    file << i << " " << zenith_luminance[i].to(watt_per_square_meter_per_sr)
         << std::endl;
  }
  file.close();
}

void Comparisons::PlotDayIrradiance(
    const std::vector<Angle>& sun_zenith,
    const std::vector<Irradiance>& sun,
    const std::vector<Irradiance>& sky) const {
  std::ofstream file(kOutputDirectory + "day_irradiance_" + name_ + ".txt");
  for (unsigned int i = 0; i < sun_zenith.size(); ++i) {
    file << i << " "
         << (sun[i] + sky[i]).to(watt_per_square_meter) << " "
         << sky[i].to(watt_per_square_meter) << std::endl;
  }
  file.close();
}

void Comparisons::SaveErrorCaption(bool png_output) {
  std::ofstream palette("output/comparisons/error_palette.txt");
  for (int i = 0; i < 41; ++i) {
    palette << kErrorColorMap[3 * i] << " " << kErrorColorMap[3 * i + 1]
        << " " << kErrorColorMap[3 * i + 2] << std::endl;
  }
  palette.close();

  std::ofstream caption("output/comparisons/error_caption.txt");
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 41; ++i) {
      caption << i << " ";
    }
    caption << std::endl;
  }
  caption.close();

  std::ofstream plot(png_output ? "output/comparisons/error_caption.plot" :
      "output/comparisons/error_caption_paper.plot");
  plot << "reset\n";
  if (png_output) {
    plot << "set terminal png size 512,80\n"
        "set output \"output/figures/error_caption.png\"\n";
  } else {
    plot << "set terminal postscript eps size 15cm,1.2cm enhanced\n"
        "set output \"paper/figures/error_caption.eps\"\n";
  }
  plot << "set xtics border out nomirror\n"
      "set xtics (\"-200\" -0.5, \"-150\" 4.5, \"-100\" 9.5, "
      "\"-50\" 14.5, \"0\" 19.5, \"50\" 24.5, \"100\" 29.5, \"150\" 34.5, "
      "\"200\" 39.5)\n"
      "unset ytics\n"
      "set palette model RGB file \"output/comparisons/error_palette.txt\" "
      "using ($1/255):($2/255):($3/255)\n"
      "unset key\n"
      "unset title\n"
      "unset colorbox\n"
      "plot \"output/comparisons/error_caption.txt\" matrix with image\n";
  plot.close();
}

void Comparisons::SaveScaleCaption(bool png_output) {
  std::ofstream palette("output/comparisons/scale_palette.txt");
  for (int i = -12; i < 6 * 64 + 12; ++i) {
    int red, green, blue;
    GetColor(pow(10.0, i / 64.0), &red, &green, &blue);
    palette << red << " " << green << " " << blue << std::endl;
  }
  palette.close();

  std::ofstream caption("output/comparisons/scale_caption.txt");
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 6 * 64 + 24; ++i) {
      caption << i << " ";
    }
    caption << std::endl;
  }
  caption.close();

  std::ofstream plot(png_output ? "output/comparisons/scale_caption.plot" :
      "output/comparisons/scale_caption_paper.plot");
  plot << "reset\n";
  if (png_output) {
    plot << "set terminal png size 768,80\n"
        "set output \"output/figures/scale_caption.png\"\n"
        "set xtics (\"1\" 11, \"10\" 75, \"10^2\" 139, \"10^3\" 203, "
        "\"10^4\" 267, \"10^5\" 331, \"10^6\" 395)\n";
  } else {
    plot << "set terminal postscript eps size 15cm,1.2cm enhanced\n"
        "set output \"paper/figures/scale_caption.eps\"\n"
        "set xtics (\"1\" 13, \"10\" 77, \"10^2\" 140, \"10^3\" 204, "
        "\"10^4\" 268, \"10^5\" 331, \"10^6\" 395)\n";
  }
  plot << "set xtics border out nomirror\n"
      "unset ytics\n"
      "set palette model RGB file \"output/comparisons/scale_palette.txt\" "
      "using ($1/255):($2/255):($3/255)\n"
      "unset key\n"
      "unset title\n"
      "unset colorbox\n"
      "plot \"output/comparisons/scale_caption.txt\" matrix with image\n";
  plot.close();
}
