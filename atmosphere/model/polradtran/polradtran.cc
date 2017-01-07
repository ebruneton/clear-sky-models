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
#include "atmosphere/model/polradtran/polradtran.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

constexpr Angle PolRadtran::kDeltaPhi;

PolRadtran::PolRadtran(const std::string& libradtran_uvspec,
    const std::string& libradtran_data, bool polarization)
    : libradtran_uvspec_(libradtran_uvspec), libradtran_data_(libradtran_data),
      polarization_(polarization) {
  // Source: solar spectrum.
  IrradianceSpectrum solar = SolarSpectrum();
  std::ofstream source("output/libradtran/solar_spectrum.txt");
  for (unsigned int i = 0; i < solar.size(); ++i) {
    source << solar.GetSample(i).to(nm) << " "
        << solar[i].to(watt_per_square_meter_per_nm) << std::endl;
  }
  source.close();

  // Wavelengths.
  std::ofstream grid("output/libradtran/wavelength_grid.txt");
  for (unsigned int i = 0; i < solar.size(); ++i) {
    grid << solar.GetSample(i).to(nm) << std::endl;
  }
  grid.close();

  // libRadtran input model file.
  std::ofstream model("output/libradtran/model.txt");
  model << "data_files_path " << libradtran_data << "\n";
  model << "\n#source irradiance\n";
  model << "source solar output/libradtran/solar_spectrum.txt per_nm\n";
  model << "wavelength " << solar.GetSample(0).to(nm) << " "
      << solar.GetSample(solar.size() - 1).to(nm) << std::endl;
  model << "wavelength_grid_file output/libradtran/wavelength_grid.txt\n";
  model << "\n#molecular properties\n";
  model << "mol_abs_param crs\n";
  model << "\n#aerosol properties\n";
  model << "aerosol_default\n";
  model << "aerosol_species_file continental_clean\n";
  model << "\n#ground properties\n";
  model << "albedo 0.1\n";
  model << "\n#solver and output options\n";
  model << "rte_solver polradtran\n";
  model << "polradtran nstokes " << (polarization ? 3 : 1) << "\n";
  model << "output_process per_nm\n";
  model << "output_user lambda uu\n";
  model << "quiet\n";
  model.close();
}

IrradianceSpectrum PolRadtran::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  IrradianceSpectrum result(0.0 * watt_per_square_meter_per_nm);
  return result;
}

RadianceSpectrum PolRadtran::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  assert(altitude == 0.0 * m);
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
}

void PolRadtran::MaybeComputeSkyDome(Angle sun_zenith) const {
  if (current_sun_zenith_ == sun_zenith) {
    return;
  }
  current_sun_zenith_ = sun_zenith;

  std::stringstream name;
  name << "output/cache/polradtran/sky_dome_" << sun_zenith.to(deg) << "_"
      << (polarization_ ? "vector.dat" : "scalar.dat");
  std::ifstream f;
  f.open(name.str());
  if (f.good()) {
    f.close();
    sky_dome_.Load(name.str());
    return;
  }

  const auto& spectrum = sky_dome_.Get(0, 0);
  for (int i = 0; i < kNumTheta; ++i) {
    std::cout << i << " / " << kNumTheta << std::endl;
    Angle view_zenith = (i + 0.5) * kDeltaPhi;
    std::ofstream input("output/libradtran/input.txt");
    input << "include output/libradtran/model.txt" << std::endl;
    input << "sza " << sun_zenith.to(deg) << std::endl;
    input << "phi0 0.0" << std::endl;
    input << "umu " << -cos(view_zenith)();
    input << std::endl << "phi";
    for (int j = 0; j < kNumPhi / 2; ++j) {
      Angle view_azimuth = (j + 0.5) * kDeltaPhi;
      input << " " << view_azimuth.to(deg);
    }
    input << std::endl;
    input.close();

    const std::string cmd =
        libradtran_uvspec_ + " -i output/libradtran/input.txt";
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

    std::stringstream string_stream(output);
    for (unsigned int j = 0; j < spectrum.size(); ++j) {
      double lambda;
      double radiance;
      double ignored;
      string_stream >> lambda;
      assert(lambda == spectrum.GetSample(j).to(nm));
      for (int k = 0; k < kNumPhi / 2 + (polarization_ ? 6 : 2); ++k) {
        string_stream >> ignored;
      }
      std::string comment;
      std::getline(string_stream, comment);
      std::getline(string_stream, comment);
      string_stream >> ignored;
      string_stream >> ignored;
      for (int k = 0; k < kNumPhi / 2; ++k) {
        string_stream >> radiance;
        sky_dome_.Get(i, k)[j] = radiance * watt_per_square_meter_per_sr_per_nm;
      }
      if (polarization_) {
        std::getline(string_stream, comment);
        std::getline(string_stream, comment);
        for (int k = 0; k < kNumPhi / 2 + 2; ++k) {
          string_stream >> ignored;
        }
        std::getline(string_stream, comment);
        std::getline(string_stream, comment);
        for (int k = 0; k < kNumPhi / 2 + 2; ++k) {
          string_stream >> ignored;
        }
      }
    }
  }
  sky_dome_.Save(name.str());
}
