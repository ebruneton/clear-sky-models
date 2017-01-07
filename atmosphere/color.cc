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
#include "atmosphere/color.h"

#include "atmosphere/atmosphere.h"
#include "math/matrix.h"

using dimensional::Matrix3;
using dimensional::Vector3;
using dimensional::Scalar;
using dimensional::vec3;

namespace {

// The 3 sample wavelengths used when simulating RGB rendering. All the models
// are run spectrally, but we can simulate RGB rendering by simply dropping the
// values at all wavelengths but 3.
constexpr Wavelength lambda_r = 680.0 * nm;
constexpr Wavelength lambda_g = 550.0 * nm;
constexpr Wavelength lambda_b = 440.0 * nm;

}  // namespace

Color GetSrgbColorNaive(const RadianceSpectrum& spectrum) {
  // Extract 3 samples of the given spectrum, to be used as RGB components. Note
  // that 'spectrum' is proportional to the input solar spectrum, whereas the
  // original O'Neal and Bruneton model implementations implicitly use a
  // wavelength independent solar spectrum. To simulate this, we need to cancel
  // the effect of the input solar spectrum when extracting the RGB values:
  static const RadianceSpectrum solar_spectrum = SolarSpectrum() * (1.0 / sr);
  Number r = spectrum(lambda_r) / solar_spectrum(lambda_r);
  Number g = spectrum(lambda_g) / solar_spectrum(lambda_g);
  Number b = spectrum(lambda_b) / solar_spectrum(lambda_b);
  // Scale the resulting color with a constant factor such that the result can
  // be compared quantitatively with the correct way of converting 'spectrum' to
  // a color, as in GetSrgbColor().
  static const auto k =
      dot(GetSrgbColor(solar_spectrum), vec3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
  return Color(k * r, k * g, k * b);
}

Color GetSrgbColorFrom3SpectrumSamples(const RadianceSpectrum& spectrum) {
  // Compute the normalization vector to properly convert 3 spectrum samples to
  // a perceptual, linear sRGB color.
  static const auto k = []() {
    WavelengthFunction<0, -3, 0, 0, 0> f;
    for (unsigned int i = 0; i < f.size(); ++i) {
      Wavelength lambda = f.GetSample(i);
      f[i] = 1.0 / (lambda * lambda * lambda);
    }
    const auto RGB = XYZ_to_sRGB * Vector3<Scalar<-2, -3, 0, 1, 0>>(
        Integral(cie_x_bar_function() * SolarSpectrum() * f),
        Integral(cie_y_bar_function() * SolarSpectrum() * f),
        Integral(cie_z_bar_function() * SolarSpectrum() * f)) *
            MaxLuminousEfficacy;
    return Vector3<Scalar<0, 1, 0, -1, 1>>(
        RGB.x / SolarSpectrum()(lambda_r) * (lambda_r * lambda_r * lambda_r),
        RGB.y / SolarSpectrum()(lambda_g) * (lambda_g * lambda_g * lambda_g),
        RGB.z / SolarSpectrum()(lambda_b) * (lambda_b * lambda_b * lambda_b));
  }();
  return Color(k.x * spectrum(lambda_r), k.y * spectrum(lambda_g),
      k.z * spectrum(lambda_b));
}

// Returns an approximate version of the given spectrum, reconstructed from 3 of
// its samples with the CIE S0, S1 and S2 components.
RadianceSpectrum GetApproximateSpectrumFrom3SpectrumSamples(
    const RadianceSpectrum& spectrum) {
  static const Matrix3<Number> sRGB_to_XYZ = inverse(XYZ_to_sRGB);
  Color XYZ = sRGB_to_XYZ * GetSrgbColorFrom3SpectrumSamples(spectrum);
  Number x = XYZ.x / (XYZ.x + XYZ.y + XYZ.z);
  Number y = XYZ.y / (XYZ.x + XYZ.y + XYZ.z);
  Radiance Y = XYZ.y / MaxLuminousEfficacy;
  Number M1 = (-1.3515 - 1.7703 * x + 5.9114 * y) /
      (0.0241 + 0.2562 * x - 0.7341 * y);
  Number M2 = (0.0300 - 31.4424 * x + 30.0717 * y) /
      (0.0241 + 0.2562 * x - 0.7341 * y);
  DimensionlessSpectrum S =
      S0_function() + S1_function() * M1 + S2_function() * M2;
  static const auto Y_S0 = Integral(S0_function() * cie_y_bar_function());
  static const auto Y_S1 = Integral(S1_function() * cie_y_bar_function());
  static const auto Y_S2 = Integral(S2_function() * cie_y_bar_function());
  SpectralRadiance k = Y / (Y_S0 + Y_S1 * M1 + Y_S2 * M2);
  return S * k;
}

Color WhiteBalanceNaive(const Color& c) {
  const Color sun_color = GetSrgbColor(SolarSpectrum() * (1.0 / sr));
  const Luminance sun_luminance =
      dot(sun_color, vec3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
  const auto white_point = sun_color / sun_luminance;
  return Color(c.x / white_point.x, c.y / white_point.y, c.z / white_point.z);
}

