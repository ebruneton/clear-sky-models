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
#include "atmosphere/model/libradtran/libradtran.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {

constexpr int kNumLayers = 100;

double AltitudeKm(int k) {
  constexpr Number M = (1.0 - exp(-(AtmosphereRadius - EarthRadius) /
      RayleighScaleHeight)) / kNumLayers;
  return k == kNumLayers ? 0.0 :
      (-RayleighScaleHeight * log(1.0 - (kNumLayers - k) * M)).to(km);
}

double Legendre(int k, double x) {
  if (k == 0) {
    return 1.0;
  } else if (k == 1) {
    return x;
  } else {
    return ((2 * k - 1) * x * Legendre(k - 1, x) -
        (k - 1) * Legendre(k - 2, x)) / k;
  }
}

double MiePhaseFunctionMoment(int k) {
  constexpr int kNumSamples = 1000;
  double sum = 0.0;
  for (int i = 0; i < kNumSamples; ++i) {
    double mu = (i + 0.5) / kNumSamples;
    sum += Legendre(k, mu) * MiePhaseFunction(mu).to(1.0 / sr);
  }
  return sum / (2.0 * kNumSamples);
}

}  // anonymous namespace

constexpr Angle LibRadtran::kDeltaPhi;

LibRadtran::LibRadtran(const std::string& libradtran_path, bool interpolate)
    : libradtran_path_(libradtran_path), interpolate_(interpolate) {
  // Source: solar spectrum.
  IrradianceSpectrum solar = SolarSpectrum();
  std::ofstream source("output/libradtran/solar_spectrum.txt");
  for (unsigned int i = 0; i < solar.size(); ++i) {
    source << solar.GetWavelength(i).to(nm) << " "
        << solar[i].to(watt_per_square_meter_per_nm) << std::endl;
  }
  source.close();

  // Atmosphere profile (altitude, pressure, temperature, density).
  std::ofstream atmosphere("output/libradtran/atmosphere.txt");
  for (int i = 0; i < kNumLayers + 1; ++i) {
    double x = exp(-AltitudeKm(i) * km / RayleighScaleHeight)();
    atmosphere << AltitudeKm(i) << " "  << 1013.0 * x << " 288 "
        << 2.545818e19 * x << std::endl;
  }
  atmosphere.close();

  // Rayleigh: scattering and absorption dtau per layer.
  ScatteringSpectrum rayleigh = RayleighScattering();
  std::ofstream molecular_scattering(
      "output/libradtran/molecular_scattering.txt");
  std::ofstream molecular_absorption(
      "output/libradtran/molecular_absorption.txt");
  for (int i = 0; i < kNumLayers + 1; ++i) {
    molecular_scattering << AltitudeKm(i) << " ";
    molecular_absorption << AltitudeKm(i) << " ";
  }
  molecular_scattering << std::endl;
  molecular_absorption << std::endl;
  for (unsigned int i = 0; i < rayleigh.size(); ++i) {
    molecular_scattering << rayleigh.GetWavelength(i).to(nm) << " ";
    molecular_absorption << rayleigh.GetWavelength(i).to(nm) << " ";
    for (int j = 0; j < kNumLayers; ++j) {
      Number dtau = rayleigh[i] * RayleighScaleHeight * (
          exp(-AltitudeKm(j + 1) * km / RayleighScaleHeight) -
          exp(-AltitudeKm(j) * km / RayleighScaleHeight));
      molecular_scattering << dtau() << " ";
      molecular_absorption << "0.0 ";
    }
    molecular_scattering << std::endl;
    molecular_absorption << std::endl;
  }
  molecular_scattering.close();
  molecular_absorption.close();

  // Aerosols: optical properties per layer.
  constexpr int kNumMoments = 20;
  double moments[kNumMoments];
  for (int i = 0; i < kNumMoments; ++i) {
    moments[i] = MiePhaseFunctionMoment(i);
  }
  ScatteringSpectrum mie_scattering = MieScattering();
  ScatteringSpectrum mie_extinction = MieExtinction();
  std::ofstream aerosol_properties("output/libradtran/aerosol_properties.txt");
  for (int i = 0; i < kNumLayers + 1; ++i) {
    std::stringstream os;
    os << "output/libradtran/aerosol_properties_" << i << ".txt";
    aerosol_properties << AltitudeKm(i) << " " << os.str() << std::endl;
    std::ofstream layer_properties(os.str());
    for (unsigned int j = 0; j < mie_scattering.size(); ++j) {
      ScatteringCoefficient extinction = i == 0 ? 0.0 / m :
          mie_extinction[j] * exp(-AltitudeKm(i) * km / MieScaleHeight);
      layer_properties << mie_scattering.GetWavelength(j).to(nm) << " ";
      layer_properties << extinction.to(1.0 / km) << " ";
      layer_properties << (mie_scattering[j] / mie_extinction[j])();
      for (int k = 0; k < kNumMoments; ++k) {
        layer_properties << " " << moments[k];
      }
      layer_properties << std::endl;
    }
    layer_properties.close();
  }
  aerosol_properties.close();

  // Ground: albedo.
  std::ofstream albedo("output/libradtran/albedo.txt");
  DimensionlessSpectrum ground_albedo = GroundAlbedo();
  for (unsigned int i = 0; i < ground_albedo.size(); ++i) {
    albedo << ground_albedo.GetWavelength(i).to(nm) << " " << ground_albedo[i]()
        << std::endl;
  }
  albedo.close();

  // libRadtran input model file.
  std::ofstream model("output/libradtran/model.txt");
  model << "#source irradiance\n";
  model << "source solar output/libradtran/solar_spectrum.txt per_nm\n";
  model << "wavelength " << solar.GetWavelength(0).to(nm) << " "
      << solar.GetWavelength(solar.size() - 1).to(nm) << std::endl;
  model << "\n#atmosphere profile\n";
  model << "atmosphere_file output/libradtran/atmosphere.txt\n";
  model << "earth_radius " << EarthRadius.to(km) << std::endl;
  model << "\n#molecular properties\n";
  model << "mol_tau_file sca output/libradtran/molecular_scattering.txt\n";
  model << "mol_tau_file abs output/libradtran/molecular_absorption.txt\n";
  model << "mol_abs_param none\n";
  model << "no_absorption mol\n";
  model << "rayleigh_depol 0.0\n";
  model << "crs_model rayleigh Penndorf\n";
  model << "\n#aerosol properties\n";
  model << "aerosol_default\n";
  model << "aerosol_file explicit output/libradtran/aerosol_properties.txt\n";
  model << "disort_intcor moments\n";
  model << "\n#ground properties\n";
  model << "albedo_file " << "output/libradtran/albedo.txt" << std::endl;
  model << "\n#solver and output options\n";
  model << "pseudospherical\n";
  model << "output_process per_nm\n";
  model << "output_user lambda uu\n";
  model << "quiet\n";

  model.close();
}

IrradianceSpectrum LibRadtran::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  IrradianceSpectrum result(0.0 * watt_per_square_meter_per_nm);
  return result;
}

RadianceSpectrum LibRadtran::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  assert(altitude == 0.0 * m);
  if (interpolate_) {
    MaybeComputeSkyDome(sun_zenith);
    if (view_sun_azimuth < 0.0 * deg) {
      view_sun_azimuth = -view_sun_azimuth;
    }
    while (view_sun_azimuth > 2.0 * pi) {
      view_sun_azimuth = view_sun_azimuth - 2.0 * pi;
    }
    if (view_sun_azimuth > pi) {
      view_sun_azimuth = 2.0 * pi - view_sun_azimuth;
    }
    Number u = view_zenith / (kNumTheta * kDeltaPhi);
    Number v = view_sun_azimuth / (kNumPhi / 2 * kDeltaPhi);
    return sky_dome_(u(), v());
  } else {
    std::cout << "sun_zenith: " << sun_zenith.to(deg) << ", view_zenith: "
        << view_zenith.to(deg) << ", azimuth: "
        << std::abs(view_sun_azimuth.to(deg)) << std::endl;
    RadianceSpectrum result(0.0 * watt_per_square_meter_per_sr_per_nm);

    std::ofstream input("output/libradtran/input.txt");
    input << "include output/libradtran/model.txt" << std::endl;
    input << "sza " << sun_zenith.to(deg) << std::endl;
    input << "phi0 0.0" << std::endl;
    input << "umu " << -cos(view_zenith)() << std::endl;
    input << "phi " << std::abs(view_sun_azimuth.to(deg)) << std::endl;
    input.close();

    const std::string cmd =
        libradtran_path_ + "bin/uvspec -i output/libradtran/input.txt";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
      return result;
    }
    char buffer[256];
    std::string output = "";
    while (!feof(pipe)) {
      if (fgets(buffer, 256, pipe) != NULL) {
        output += buffer;
      }
    }
    pclose(pipe);

    std::stringstream string_stream(output);
    for (unsigned int i = 0; i < result.size(); ++i) {
      double lambda;
      double radiance;
      string_stream >> lambda >> radiance;
      assert(lambda == result.GetWavelength(i).to(nm));
      result[i] = radiance * watt_per_square_meter_per_sr_per_nm;
    }
    return result;
  }
}

void LibRadtran::MaybeComputeSkyDome(Angle sun_zenith) const {
  if (current_sun_zenith_ == sun_zenith) {
    return;
  }
  current_sun_zenith_ = sun_zenith;

  std::ofstream input("output/libradtran/input.txt");
  input << "include output/libradtran/model.txt" << std::endl;
  input << "sza " << sun_zenith.to(deg) << std::endl;
  input << "phi0 0.0" << std::endl;
  input << "umu";
  for (int i = 0; i < kNumTheta; ++i) {
    Angle view_zenith = (i + 0.5) * kDeltaPhi;
    input << " " << -cos(view_zenith)();
  }
  input << std::endl << "phi";
  for (int i = 0; i < kNumPhi / 2; ++i) {
    Angle view_azimuth = (i + 0.5) * kDeltaPhi;
    input << " " << view_azimuth.to(deg);
  }
  input << std::endl;
  input.close();

  const std::string cmd =
      libradtran_path_ + "bin/uvspec -i output/libradtran/input.txt";
  FILE* pipe = popen(cmd.c_str(), "r");
  if (!pipe) {
    return;
  }
  char buffer[256];
  std::string output = "";
  while (!feof(pipe)) {
    if (fgets(buffer, 256, pipe) != NULL) {
      output += buffer;
    }
  }
  pclose(pipe);

  const auto& spectrum = sky_dome_.Get(0, 0);
  std::stringstream string_stream(output);
  for (unsigned int i = 0; i < spectrum.size(); ++i) {
    double lambda;
    double radiance;
    string_stream >> lambda;
    assert(lambda == spectrum.GetWavelength(i).to(nm));
    for (int j = 0; j < kNumTheta; ++j) {
      for (int k = 0; k < kNumPhi / 2; ++k) {
        string_stream >> radiance;
        sky_dome_.Get(j, k)[i] = radiance * watt_per_square_meter_per_sr_per_nm;
      }
    }
  }
}
