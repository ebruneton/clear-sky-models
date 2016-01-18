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
#ifndef ATMOSPHERE_MEASUREMENT_MEASURED_ATMOSPHERE_H_
#define ATMOSPHERE_MEASUREMENT_MEASURED_ATMOSPHERE_H_

#include <string>

#include "atmosphere/atmosphere.h"
#include "math/angle.h"
#include "math/hemispherical_function.h"
#include "physics/spectrum.h"
#include "physics/units.h"

class MeasuredAtmosphere : public Atmosphere {
 public:
  MeasuredAtmosphere(const std::string& directory, const std::string& date,
      const std::string& hour, const std::string& minutes, Angle sun_zenith,
      Angle sun_azimuth, const std::string& cache_directory,
      bool compute_azimuth_from_data = true);

  inline const Angle sun_zenith() const { return sun_zenith_; }
  inline const Angle sun_azimuth() const { return sun_azimuth_; }

  IrradianceSpectrum GetSunIrradiance(Length altitude,
      Angle sun_zenith) const;

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun_azimuth) const;

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle sun_azimuth, Angle view_zenith, Angle view_azimuth) const;

  // Same as GetSkyRadiance, but returns 0 for view directions that were not
  // measured directly, instead of interpolating the neasest samples.
  RadianceSpectrum GetSkyRadianceMeasurement(Length altitude, Angle sun_zenith,
      Angle sun_azimuth, Angle view_zenith, Angle view_azimuth) const;

 private:
  Radiance GetAxiSymmetryError(Angle axis_azimuth) const;
  Angle EstimateSunAzimuth(Angle initial_azimuth) const;

  Angle sun_zenith_;
  Angle sun_azimuth_;
  HemisphericalFunction<RadianceSpectrum> measurements_;
};

#endif  // ATMOSPHERE_MEASUREMENT_MEASURED_ATMOSPHERE_H_
