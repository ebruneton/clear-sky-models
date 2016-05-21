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
  // Extracting 3 samples from a spectrum of the form c0 S0 + c1 S1 + c2 S2, at
  // the lambda_r, lambda_g, lambda_b wavelengths, gives the values
  //   sample_r = c0 S0(lambda_r) + c1 S1(lambda_r) + c2 S2(lambda_r)
  //   sample_g = c0 S0(lambda_g) + c1 S1(lambda_g) + c2 S2(lambda_g)
  //   sample_b = c0 S0(lambda_b) + c1 S1(lambda_b) + c2 S2(lambda_b)
  // Posing that sample_i is equal to spectrum(lambda_i), for each i, gives 3
  // equations for the 3 unknowns c0, c1 and c2. We can then compute c0, c1, c2
  // with the inverse of the above matrix. From this, we can compute the XYZ
  // color corresponding to the reconstructed spectrum c0 S0 + c1 S1 + c2 S2:
  //   X = c0 X0 + c1 X1 + c2 X2
  //   Y = c0 Y0 + c1 Y1 + c2 Y2
  //   Z = c0 Z0 + c1 Z1 + c2 Z2
  // where X0,Y0,Z0 is the XYZ color corresponding to the S0 spectrum, and
  // similarly for X1,Y1,Z1 and X2,Y2,Z2. Finally, using the XYZ_to_sRGB matrix,
  // we get the linear sRGB color corresponding to the reconstructed spectrum.
  // Putting this all together, we see that the sRGB color is obtained from the
  // 3 spectrum samples with a matrix vector product, where the matrix can be
  // precomputed and is the product of XYZ_to_sRGB and the two above matrices:
  /*static const auto conversion_matrix = XYZ_to_sRGB * Matrix3<Wavelength>(
      Integral(S0_function() * cie_x_bar_function()),
      Integral(S1_function() * cie_x_bar_function()),
      Integral(S2_function() * cie_x_bar_function()),
      Integral(S0_function() * cie_y_bar_function()),
      Integral(S1_function() * cie_y_bar_function()),
      Integral(S2_function() * cie_y_bar_function()),
      Integral(S0_function() * cie_z_bar_function()),
      Integral(S1_function() * cie_z_bar_function()),
      Integral(S2_function() * cie_z_bar_function())
  ) * inverse(Matrix3<Number>(
      S0_function()(lambda_r), S1_function()(lambda_r), S2_function()(lambda_r),
      S0_function()(lambda_g), S1_function()(lambda_g), S2_function()(lambda_g),
      S0_function()(lambda_b), S1_function()(lambda_b), S2_function()(lambda_b)
  ));
  Vector3<SpectralRadiance> samples(
      spectrum(lambda_r), spectrum(lambda_g), spectrum(lambda_b));
  return (conversion_matrix * samples) * MaxLuminousEfficacy;*/

  // Compute the normalization vector to properly convert 3 spectrum samples to
  // a perceptual, linear sRGB color.
  static const auto k = []() {
    WavelengthFunction<Scalar<0, -3, 0, 0, 0>> f;
    for (unsigned int i = 0; i < f.size(); ++i) {
      Wavelength lambda = f.GetWavelength(i);
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
  // Extracting 3 samples from a spectrum of the form c0 S0 + c1 S1 + c2 S2, at
  // the lambda_r, lambda_g, lambda_b wavelengths, gives the values
  //   sample_r = c0 S0(lambda_r) + c1 S1(lambda_r) + c2 S2(lambda_r)
  //   sample_g = c0 S0(lambda_g) + c1 S1(lambda_g) + c2 S2(lambda_g)
  //   sample_b = c0 S0(lambda_b) + c1 S1(lambda_b) + c2 S2(lambda_b)
  // Posing that sample_i is equal to spectrum(lambda_i), for each i, gives 3
  // equations for the 3 unknowns c0, c1 and c2. We can then compute c0, c1, c2
  // with the inverse of the above matrix:
  /*static const Matrix3<Number> conversion_matrix = inverse(Matrix3<Number>(
      S0_function()(lambda_r), S1_function()(lambda_r), S2_function()(lambda_r),
      S0_function()(lambda_g), S1_function()(lambda_g), S2_function()(lambda_g),
      S0_function()(lambda_b), S1_function()(lambda_b), S2_function()(lambda_b)
  ));
  Vector3<SpectralRadiance> samples(
      spectrum(lambda_r), spectrum(lambda_g), spectrum(lambda_b));
  auto coefficients = conversion_matrix * samples;
  // and return the reconstructed, approximate spectrum c0 S0 + c1 S1 + c2 S2:
  return S0_function() * coefficients.x + S1_function() * coefficients.y +
      S2_function() * coefficients.z;*/

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
