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
#ifndef ATMOSPHERE_MODEL_BRUNETON_BRUNETON_H_
#define ATMOSPHERE_MODEL_BRUNETON_BRUNETON_H_

#include "atmosphere/atmosphere.h"
#include "atmosphere/model/bruneton/core.h"
#include "math/angle.h"
#include "physics/spectrum.h"
#include "physics/units.h"

class Bruneton : public Atmosphere {
 public:
  enum ScatteringType {
    SINGLE_SCATTERING_ONLY, DOUBLE_SCATTERING_ONLY, ALL_ORDERS
  };

  // Use 3 wavelength for the Bruneton model, 15 for the Elek model.
  Bruneton(ScatteringType scattering_type, int original_number_of_wavelength);

  int GetOriginalNumberOfWavelengths() const override {
    return original_number_of_wavelength_;
  }

  IrradianceSpectrum GetSunIrradiance(Length altitude,
      Angle sun_zenith) const override;

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun_azimuth) const override;

  IrradianceSpectrum GetSkyIrradiance(Length altitude,
      Angle sun_zenith) const override;

 private:
  int original_number_of_wavelength_;
  TransmittanceTexture transmittance_sampler_;
  IrradianceTexture inscatter1R_sampler_;
  IrradianceTexture inscatter1M_sampler_;
  RadianceTexture inscatterN_sum_sampler_;
  SkyIrradianceTexture sky_irradiance_sum_sampler_;
};

#endif  // ATMOSPHERE_MODEL_BRUNETON_BRUNETON_H_
