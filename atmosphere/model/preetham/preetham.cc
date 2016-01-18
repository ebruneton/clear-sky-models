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
#include "atmosphere/model/preetham/preetham.h"

#include <vector>

#include "physics/cie.h"

namespace {

class PerezFunction {
 public:
  PerezFunction(Number a, Number b, Number c, Number d, Number e)
      : a_(a), b_(b), c_(c), d_(d), e_(e) {}

  Number operator()(Angle sun_zenith, Angle view_zenith, Angle view_sun) const {
    return (*this)(view_zenith, view_sun) / (*this)(0.0 * rad, sun_zenith);
  }

 private:
  Number operator()(Angle view_zenith, Angle view_sun) const {
    return (1.0 + a_ * exp(b_ / cos(view_zenith))) *
        (1.0 + c_ * exp(d_ * view_sun.to(rad)) +
         e_ * cos(view_sun) * cos(view_sun));
  }

  const Number a_;
  const Number b_;
  const Number c_;
  const Number d_;
  const Number e_;
};

struct Tables {
  std::vector<Dimensionless> S0_samples;
  std::vector<Dimensionless> S1_samples;
  std::vector<Dimensionless> S2_samples;
  std::vector<ScatteringCoefficient> ozone_samples;
  std::vector<ScatteringCoefficient> vapor_samples;
  std::vector<Dimensionless> gas_samples;

  Tables() {
    // Values from Table 2 in Preetham 1999 and
    // CIE TR "Colorimetry" 15:2004 Third Edition (ISBN 3 901 906 33 9)
    Add(+61.5,  38.0,  5.3,  0.0,   0.0,  0.0);
    Add(+68.8,  42.4,  6.1,  0.0,   0.0,  0.0);
    Add(+63.4,  38.5,  3.0,  0.0,   0.0,  0.0);
    Add(+65.8,  35.0,  1.2,  0.0,   0.0,  0.0);
    Add(+94.8,  43.4, -1.1,  0.0,   0.0,  0.0);
    Add(104.8,  46.3, -0.5,  0.0,   0.0,  0.0);
    Add(105.9,  43.9, -0.7,  0.0,   0.0,  0.0);
    Add(+96.8,  37.1, -1.2,  0.0,   0.0,  0.0);
    Add(113.9,  36.7, -2.6,  0.0,   0.0,  0.0);
    Add(125.6,  35.9, -2.9,  0.3,   0.0,  0.0);
    Add(125.5,  32.6, -2.8,  0.6,   0.0,  0.0);
    Add(121.3,  27.9, -2.6,  0.9,   0.0,  0.0);
    Add(121.3,  24.3, -2.6,  1.4,   0.0,  0.0);
    Add(113.5,  20.1, -1.8,  2.1,   0.0,  0.0);
    Add(113.1,  16.2, -1.5,  3.0,   0.0,  0.0);
    Add(110.8,  13.2, -1.3,  4.0,   0.0,  0.0);
    Add(106.5,   8.6, -1.2,  4.8,   0.0,  0.0);
    Add(108.8,   6.1, -1.0,  6.3,   0.0,  0.0);
    Add(105.3,   4.2, -0.5,  7.5,   0.0,  0.0);
    Add(104.4,   1.9, -0.3,  8.5,   0.0,  0.0);
    Add(100.0,   0.0,  0.0, 10.3,   0.0,  0.0);
    Add(+96.0,  -1.6,  0.2, 12.0,   0.0,  0.0);
    Add(+95.1,  -3.5,  0.5, 12.0,   0.0,  0.0);
    Add(+89.1,  -3.5,  2.1, 11.5,   0.0,  0.0);
    Add(+90.5,  -5.8,  3.2, 12.5,   0.0,  0.0);
    Add(+90.3,  -7.2,  4.1, 12.0,   0.0,  0.0);
    Add(+88.4,  -8.6,  4.7, 10.5,   0.0,  0.0);
    Add(+84.0,  -9.5,  5.1,  9.0,   0.0,  0.0);
    Add(+85.1, -10.9,  6.7,  7.9,   0.0,  0.0);
    Add(+81.9, -10.7,  7.3,  6.7,   0.0,  0.0);
    Add(+82.6, -12.0,  8.6,  5.7,   0.0,  0.0);
    Add(+84.9, -14.0,  9.8,  4.8,   0.0,  0.0);
    Add(+81.3, -13.6, 10.2,  3.6,   0.0,  0.0);
    Add(+71.9, -12.0,  8.3,  2.8,   1.6,  0.0);
    Add(+74.3, -13.3,  9.6,  2.3,   2.4,  0.0);
    Add(+76.4, -12.9,  8.5,  1.8,  1.25,  0.0);
    Add(+63.3, -10.6,  7.0,  1.4, 100.0,  0.0);
    Add(+71.7, -11.6,  7.6,  1.1,  87.0,  0.0);
    Add(+77.0, -12.2,  8.0,  1.0,   6.1,  0.0);
    Add(+65.2, -10.2,  6.7,  0.9,   0.1,  0.0);
    Add(+47.7,  -7.8,  5.2,  0.7, 0.001,  3.0);
    Add(+68.6, -11.2,  7.4,  0.4, 0.001, 0.21);
    Add(+65.0, -10.4,  6.8,  0.0,  0.06,  0.0);
    Add(+66.0, -10.6,  7.0,  0.0,   0.0,  0.0);
    Add(+61.0,  -9.7,  6.4,  0.0,   0.0,  0.0);
    Add(+53.3,  -8.3,  5.5,  0.0,   0.0,  0.0);
    Add(+58.9,  -9.3,  6.1,  0.0,   0.0,  0.0);
    Add(+61.9,  -9.8,  6.5,  0.0,   0.0,  0.0);
  }

 private:
  void Add(double S0_sample, double S1_sample, double S2_sample,
      double ozone_sample, double vapor_sample, double gas_sample) {
    S0_samples.push_back(S0_sample);
    S1_samples.push_back(S1_sample);
    S2_samples.push_back(S2_sample);
    ozone_samples.push_back(ozone_sample / m);
    vapor_samples.push_back(vapor_sample / m);
    gas_samples.push_back(gas_sample);
  }
};

const Tables tables = Tables();
const DimensionlessSpectrum S0(360.0 * nm, 830.0 * nm, tables.S0_samples);
const DimensionlessSpectrum S1(360.0 * nm, 830.0 * nm, tables.S1_samples);
const DimensionlessSpectrum S2(360.0 * nm, 830.0 * nm, tables.S2_samples);
const ScatteringSpectrum k_ozone(360.0 * nm, 830.0 * nm, tables.ozone_samples);
const ScatteringSpectrum k_vapor(360.0 * nm, 830.0 * nm, tables.vapor_samples);
const DimensionlessSpectrum k_gas(360.0 * nm, 830.0 * nm, tables.gas_samples);

// Since we don't have a separate dimension for candela or lumen, we reuse the
// power dimension here for lumens. Be careful when using these photometric
// units together with the radiometric units used everewhere else!
typedef Scalar<0, 0, 0, 1> Lumen;
typedef Scalar<0, 0, -1, 1> Candela;
typedef Scalar<-2, 0, -1, 1> Luminance;
constexpr Lumen lm = Lumen::Unit();
constexpr Candela cd = lm / sr;
constexpr Candela kcd = 1000.0 * cd;

}  // anonymous namespace

Preetham::Preetham(double turbidity, bool use_lbl_zenith_luminance)
    : turbidity_(turbidity),
      use_lbl_zenith_luminance_(use_lbl_zenith_luminance) {}

IrradianceSpectrum Preetham::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  assert(altitude == 0.0 * m);
  IrradianceSpectrum result = SolarSpectrum();
  // Formulas from Appendix 1 in Preetham 1999.
  Number air_mass =
      1.0 / (cos(sun_zenith) + 0.15 * pow(93.885 - sun_zenith.to(deg), -1.253));
  for (unsigned int i = 0; i < result.size(); ++i) {
    Wavelength lambda = result.GetWavelength(i);
    constexpr Wavelength lambda0 = 1000.0 * nm;
    // Transmittance due to Rayleigh scattering.
    Number rayleigh_transmittance =
        exp(-0.008735 * air_mass * pow(lambda / lambda0, -4.08));
    // Transmittance due to Mie scattering.
    double beta = 0.04608 * turbidity_ - 0.04586;
    constexpr double alpha = 1.3;
    Number mie_transmittance =
        exp(-beta * air_mass * pow(lambda / lambda0, -alpha));
    // Transmittance of the ozone layer.
    constexpr Length l = 0.0035 * m;
    Number ozone_transmittance = exp(-air_mass * l *  k_ozone[i]);
    // Transmittance of the other gases.
    Number gas_transmittance = exp(-1.41 * air_mass * k_gas[i] /
        pow(1.0 + 118.93 * k_gas[i] * air_mass, 0.45));
    // Transmittance of the water vapor.
    constexpr Length w = 0.02 * m;
    Number vapor_transmittance = exp(-0.2385 * k_vapor[i] * w * air_mass /
        pow(1.0 + 20.07 * k_vapor[i] * w * air_mass, 0.45));
    // The ground irradiance is the extraterrestrial solar irradiance times the
    // total transmittance.
    result[i] = result[i] * rayleigh_transmittance * mie_transmittance *
        ozone_transmittance * gas_transmittance * vapor_transmittance;
  }
  return result;
}

RadianceSpectrum Preetham::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  assert(altitude == 0.0 * m);
  // Values and formulas from Appendix 2 of Preetham 1999.
  PerezFunction Y_function(
       0.1787 * turbidity_ - 1.4630,
      -0.3554 * turbidity_ + 0.4275,
      -0.0227 * turbidity_ + 5.3251,
       0.1206 * turbidity_ - 2.5771,
      -0.0670 * turbidity_ + 0.3703);
  PerezFunction x_function(
      -0.0193 * turbidity_ - 0.2592,
      -0.0665 * turbidity_ + 0.0008,
      -0.0004 * turbidity_ + 0.2125,
      -0.0641 * turbidity_ - 0.8989,
      -0.0033 * turbidity_ + 0.0452);
  PerezFunction y_function(
      -0.0167 * turbidity_ - 0.2608,
      -0.0950 * turbidity_ + 0.0092,
      -0.0079 * turbidity_ + 0.2102,
      -0.0441 * turbidity_ - 1.6537,
      -0.0109 * turbidity_ + 0.0529);
  double t = sun_zenith.to(rad);
  double t2 = t * t;
  double t3 = t * t2;
  double T = turbidity_;
  double T2 = T * T;
  // Absolute zenith luminance and zenith chromaticity values:
  Luminance Y_zenith;
  if (use_lbl_zenith_luminance_) {
    // See Eq. 1 in "Zenith luminance and sky luminance distributions for
    // daylighting calculations", Karayel et al. in Energy and Buildings, 1984.
    Y_zenith = ((1.376 * T - 1.81) * tan(PI / 2.0 - t) + 0.38) * kcd / m2;
  } else {
    // Formula from Preetham 1999 paper.
    double khi = (4.0 / 9.0 - T / 120.0) * (PI - 2.0 * t);
    Y_zenith =
        ((4.0453 * T - 4.9710) * tan(khi) - 0.2155 * T + 2.4192) * kcd / m2;
  }
  double x_zenith =
      (+0.00166 * t3 - 0.00375 * t2 + 0.00209 * t + 0.00000) * T2 +
      (-0.02903 * t3 + 0.06377 * t2 - 0.03202 * t + 0.00394) * T +
      (+0.11693 * t3 - 0.21196 * t2 + 0.06052 * t + 0.25886);
  double y_zenith =
      (+0.00275 * t3 - 0.00610 * t2 + 0.00317 * t + 0.00000) * T2 +
      (-0.04214 * t3 + 0.08970 * t2 - 0.04153 * t + 0.00516) * T +
      (+0.15346 * t3 - 0.26756 * t2 + 0.06670 * t + 0.26688);
  // Absolute luminance and chromaticity values in the view direction:
  Angle view_sun = GetViewSunAngle(sun_zenith, view_zenith, view_sun_azimuth);
  Luminance Y = Y_zenith * Y_function(sun_zenith, view_zenith, view_sun);
  Number x = x_zenith * x_function(sun_zenith, view_zenith, view_sun);
  Number y = y_zenith * y_function(sun_zenith, view_zenith, view_sun);
  // Relative spectral power distribution (unitless).
  // Values and formulas from Appendix 6 of Preetham 1999.
  Number M1 = (-1.3515 - 1.7703 * x + 5.9114 * y) /
      (0.0241 + 0.2562 * x - 0.7341 * y);
  Number M2 = (0.0300 - 31.4424 * x + 30.0717 * y) /
      (0.0241 + 0.2562 * x - 0.7341 * y);
  DimensionlessSpectrum S = S0 + S1 * M1 + S2 * M2;
  // The absolute spectral power distribution L (in W/m^2/sr/nm) we want to
  // return is the relative spectral power distribution S times some unknown
  // scalar k (also in W/m^2/sr/nm): L = k*S. In order to compute k, we use the
  // relation between the spectral radiance L and the luminance Y, namely:
  //
  //   Y = 683.0 lm / W * integral of y_bar(lambda)*L(lambda) dlambda
  //
  // where y_bar is the CIE photoptic luminosity function (see
  // https://en.wikipedia.org/wiki/Luminosity_function#Details). Substituting
  // L = k*S = k*(S0 + S1 * M1 + S2 * M2) in the above equation, and noting:
  //
  static const auto Y_S0 =
      683.0 * (lm / watt) * Integral(S0 * cie_y_bar_function());
  static const auto Y_S1 =
      683.0 * (lm / watt) * Integral(S1 * cie_y_bar_function());
  static const auto Y_S2 =
      683.0 * (lm / watt) * Integral(S2 * cie_y_bar_function());
  //
  // we get:
  //
  SpectralRadiance k = Y / (Y_S0 + Y_S1 * M1 + Y_S2 * M2);
  return S * k;
}
