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
#include "atmosphere/model/bruneton/bruneton.h"

#include <sstream>
#include <string>

#include "atmosphere/model/bruneton/core.h"

namespace {

std::string ToString(int i) {
  std::stringstream f;
  f << i;
  return f.str();
}

}  // anonymous namespace

Bruneton::Bruneton(ScatteringType scattering_type,
    int original_number_of_wavelength)
        : original_number_of_wavelength_(original_number_of_wavelength) {
  std::ifstream f;
  std::string name;
  const std::string cache_directory = "output/cache/bruneton/";

  name = cache_directory + "transmittance.dat";
  f.open(name);
  if (f.good()) {
    f.close();
    std::cout << "Loading..." << std::endl;
    transmittance_sampler_.Load(name);
  } else {
    ComputeTransmittance(&transmittance_sampler_);
    transmittance_sampler_.Save(name);
  }

  if (scattering_type == DOUBLE_SCATTERING_ONLY) {
    inscatter1R_sampler_ = IrradianceTexture(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
    inscatter1M_sampler_ = IrradianceTexture(
        IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm));
  } else {
    name = cache_directory + "inscatter1R.dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      inscatter1R_sampler_.Load(name);
      inscatter1M_sampler_.Load(cache_directory + "inscatter1M.dat");
    } else {
      ComputeInscatter1(transmittance_sampler_, &inscatter1R_sampler_,
          &inscatter1M_sampler_);
      inscatter1R_sampler_.Save(name);
      inscatter1M_sampler_.Save(cache_directory + "inscatter1M.dat");
    }
  }

  if (scattering_type == SINGLE_SCATTERING_ONLY) {
    inscatterN_sum_sampler_ = RadianceTexture(
        RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm));
    name = cache_directory + "irradiance2.dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      sky_irradiance_sum_sampler_.Load(name);
    } else {
      std::cerr << name << " must be precomputed. Run with ALL_ORDERS first."
                << std::endl;
      exit(-1);
    }
    return;
  } else if (scattering_type == DOUBLE_SCATTERING_ONLY) {
    name = cache_directory + "inscatter2.dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      inscatterN_sum_sampler_.Load(name);
      sky_irradiance_sum_sampler_.Load(cache_directory + "irradiance3.dat");
    } else {
      std::cerr << name << " must be precomputed. Run with ALL_ORDERS first."
                << std::endl;
      exit(-1);
    }
    return;
  } else {
    name = cache_directory + "inscatterNSum.dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      inscatterN_sum_sampler_.Load(name);
      sky_irradiance_sum_sampler_.Load(cache_directory + "irradianceNSum.dat");
      return;
    }
  }

  SkyIrradianceTexture sky_irradiance_sampler;
  name = cache_directory + "irradiance1.dat";
  f.open(name);
  if (f.good()) {
    f.close();
    std::cout << "Loading..." << std::endl;
    sky_irradiance_sampler.Load(name);
  } else {
    ComputeSkyIrradiance1(transmittance_sampler_, &sky_irradiance_sampler);
    sky_irradiance_sampler.Save(name);
  }

  RadianceDensityTexture inscatterS_sampler;
  RadianceTexture inscatterN_sampler;
  for (int i = 2; i <= 4; ++i) {
    const std::string iteration = ToString(i);
    bool first_iteration = i == 2;

    name = cache_directory + "inscatterS" + iteration + ".dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      inscatterS_sampler.Load(name);
    } else {
      ComputeInscatterS(transmittance_sampler_, sky_irradiance_sampler,
          inscatter1R_sampler_, inscatter1M_sampler_, inscatterN_sampler,
          first_iteration, &inscatterS_sampler);
      inscatterS_sampler.Save(name);
    }

    name = cache_directory + "irradiance" + iteration + ".dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      sky_irradiance_sampler.Load(name);
    } else {
      ComputeSkyIrradianceN(inscatter1R_sampler_, inscatter1M_sampler_,
          inscatterN_sampler, first_iteration, &sky_irradiance_sampler);
      sky_irradiance_sampler.Save(name);
    }

    name = cache_directory + "inscatter" + iteration + ".dat";
    f.open(name);
    if (f.good()) {
      f.close();
      std::cout << "Loading..." << std::endl;
      inscatterN_sampler.Load(name);
    } else {
      ComputeInscatterN(transmittance_sampler_, inscatterS_sampler,
          &inscatterN_sampler);
      inscatterN_sampler.Save(name);
    }

    if (first_iteration) {
      sky_irradiance_sum_sampler_ = sky_irradiance_sampler;
      inscatterN_sum_sampler_ = inscatterN_sampler;
    } else {
      sky_irradiance_sum_sampler_ += sky_irradiance_sampler;
      inscatterN_sum_sampler_ += inscatterN_sampler;
    }
  }
  inscatterN_sum_sampler_.Save(cache_directory + "inscatterNSum.dat");
  sky_irradiance_sum_sampler_.Save(cache_directory + "irradianceNSum.dat");
}

IrradianceSpectrum Bruneton::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  return GetTransmittance(transmittance_sampler_, EarthRadius + altitude,
      cos(sun_zenith)) * SolarSpectrum();
}

RadianceSpectrum Bruneton::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  Number mu_s = cos(sun_zenith);
  Number mu = cos(view_zenith);
  Number nu =
      cos(view_sun_azimuth) * sin(view_zenith) * sin(sun_zenith) + mu * mu_s;
  InverseSolidAngle rayleigh_phase = RayleighPhaseFunction(nu);
  InverseSolidAngle mie_phase = MiePhaseFunction(nu);

  IrradianceSpectrum ray =
      texture4d(inscatter1R_sampler_, EarthRadius + altitude, mu, mu_s, nu);
  IrradianceSpectrum mie =
      texture4d(inscatter1M_sampler_, EarthRadius + altitude, mu, mu_s, nu);
  RadianceSpectrum raymieN =
      texture4d(inscatterN_sum_sampler_, EarthRadius + altitude, mu, mu_s, nu);
  return (ray * rayleigh_phase + mie * mie_phase) + raymieN;
}

IrradianceSpectrum Bruneton::GetSkyIrradiance(Length altitude,
    Angle sun_zenith) const {
  return ::GetSkyIrradiance(
      sky_irradiance_sum_sampler_, EarthRadius + altitude, cos(sun_zenith));
}
