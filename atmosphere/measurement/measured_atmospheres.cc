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
#include "atmosphere/measurement/measured_atmospheres.h"

void MeasuredAtmospheres::AddAtmosphere(const MeasuredAtmosphere* atmosphere) {
  measurements_.push_back(
      std::unique_ptr<const MeasuredAtmosphere>(atmosphere));
}

IrradianceSpectrum MeasuredAtmospheres::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  for (unsigned int i = 0; i < measurements_.size(); ++i) {
    if (measurements_[i]->sun_zenith() == sun_zenith) {
      return measurements_[i]->GetSunIrradiance(altitude, sun_zenith);
    }
  }
  assert(false);
  return IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm);
}

RadianceSpectrum MeasuredAtmospheres::GetSkyRadiance(Length altitude,
    Angle sun_zenith, Angle view_zenith, Angle view_sun_azimuth) const {
  for (unsigned int i = 0; i < measurements_.size(); ++i) {
    if (measurements_[i]->sun_zenith() == sun_zenith) {
      return measurements_[i]->GetSkyRadiance(
          altitude, sun_zenith, view_zenith, view_sun_azimuth);
    }
  }
  assert(false);
  return RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm);
}

RadianceSpectrum MeasuredAtmospheres::GetSkyRadiance(Length altitude,
    Angle sun_zenith, Angle sun_azimuth, Angle view_zenith,
    Angle view_azimuth) const {
  for (unsigned int i = 0; i < measurements_.size(); ++i) {
    if (measurements_[i]->sun_zenith() == sun_zenith &&
        measurements_[i]->sun_azimuth() == sun_azimuth) {
      return measurements_[i]->GetSkyRadiance(
          altitude, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
    }
  }
  assert(false);
  return RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm);
}

RadianceSpectrum MeasuredAtmospheres::GetSkyRadianceMeasurement(Length altitude,
    Angle sun_zenith, Angle sun_azimuth, Angle view_zenith,
    Angle view_azimuth) const {
  for (unsigned int i = 0; i < measurements_.size(); ++i) {
    if (measurements_[i]->sun_zenith() == sun_zenith &&
        measurements_[i]->sun_azimuth() == sun_azimuth) {
      return measurements_[i]->GetSkyRadianceMeasurement(
          altitude, sun_zenith, sun_azimuth, view_zenith, view_azimuth);
    }
  }
  assert(false);
  return RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm);
}

IrradianceSpectrum MeasuredAtmospheres::GetSkyIrradiance(Length altitude,
    Angle sun_zenith) const {
  for (unsigned int i = 0; i < measurements_.size(); ++i) {
    if (measurements_[i]->sun_zenith() == sun_zenith) {
      return measurements_[i]->GetSkyIrradiance(altitude, sun_zenith);
    }
  }
  assert(false);
  return IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm);
}
