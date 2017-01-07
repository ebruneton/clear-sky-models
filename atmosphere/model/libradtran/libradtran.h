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
#ifndef ATMOSPHERE_MODEL_LIBRADTRAN_LIBRADTRAN_H_
#define ATMOSPHERE_MODEL_LIBRADTRAN_LIBRADTRAN_H_

#include <string>

#include "atmosphere/atmosphere.h"
#include "atmosphere/hemispherical_function.h"
#include "atmosphere/measurement/measured_atmospheres.h"
#include "math/angle.h"
#include "math/binary_function.h"
#include "physics/units.h"

class LibRadtran : public Atmosphere {
 public:
  enum CacheType {
    BINARY_FUNCTION_CACHE, HEMISPHERICAL_FUNCTION_CACHE
  };

  LibRadtran(const std::string& libradtran_uvspec, CacheType cache_type);

  LibRadtran(const std::string& libradtran_uvspec, double mie_angstrom_alpha,
      double mie_angstrom_beta, double mie_phase_function_g, bool ground_albedo,
      CacheType cache_type);

  // Not implemented.
  virtual IrradianceSpectrum GetSunIrradiance(Length altitude,
      Angle sun_zenith) const;

  virtual RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun_azimuth) const;

  virtual RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle sun_azimuth, Angle view_zenith, Angle view_azimuth) const;

 private:
  static constexpr int kNumPhi = 120;
  static constexpr int kNumTheta = kNumPhi / 4;
  static constexpr Angle kDeltaPhi = 2.0 * pi / kNumPhi;

  void MaybeComputeBinaryFunctionCache(Angle sun_zenith) const;
  void MaybeComputeHemisphericalFunctionCache(Angle sun_zenith,
      Angle sun_azimuth) const;

  std::string libradtran_uvspec_;
  CacheType cache_type_;
  mutable Angle current_sun_zenith_;
  mutable Angle current_sun_azimuth_;
  mutable dimensional::BinaryFunction<kNumTheta, kNumPhi / 2, RadianceSpectrum>
      binary_function_cache_;
  mutable HemisphericalFunction<RadianceSpectrum> hemispherical_function_cache_;
};

#endif  // ATMOSPHERE_MODEL_LIBRADTRAN_LIBRADTRAN_H_
