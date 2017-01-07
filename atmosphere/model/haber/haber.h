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
#ifndef ATMOSPHERE_MODEL_HABER_HABER_H_
#define ATMOSPHERE_MODEL_HABER_HABER_H_

#include <vector>

#include "atmosphere/atmosphere.h"
#include "math/angle.h"
#include "math/binary_function.h"
#include "math/vector.h"
#include "physics/units.h"

class Haber : public Atmosphere {
 public:
  enum ScatteringType {
    SINGLE_SCATTERING_ONLY, DOUBLE_SCATTERING_ONLY, ALL_ORDERS
  };

  explicit Haber(ScatteringType scattering_type);

  int GetOriginalNumberOfWavelengths() const override { return 8; }

  IrradianceSpectrum GetSunIrradiance(Length altitude,
      Angle sun_zenith) const override;

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun_azimuth) const override;

 private:
  static constexpr int kNumLayers = 50;
  static constexpr Length kMaxLayerHeight = 35.0 * km;

  static constexpr int kNumPhi = 72;
  static constexpr Angle kDeltaPhi = 2.0 * pi / kNumPhi;

  static constexpr Length kMinShellRadius = 10.0 * m;
  static constexpr Length kMaxShellRadius = sqrt(
      kMaxLayerHeight * kMaxLayerHeight + 2.0 * EarthRadius * kMaxLayerHeight);
  static constexpr Number kShellRatio =
      (1.0 + PI / kNumPhi) / (1.0 - PI / kNumPhi);
  static constexpr int kNumShell =
      ceil((log(kMaxShellRadius / kMinShellRadius) / log(kShellRatio))());

  static constexpr int kNumTheta =
     ceil(((pi / 2.0 + 2.0 * asin(kMaxShellRadius / (2.0 * EarthRadius))) /
         kDeltaPhi)());

  typedef dimensional::Vector3<Length> Position;
  typedef dimensional::Vector3<Number> Direction;

  struct Cell {
    Position center;
    Volume volume;
    Length radial_width;
    int shell_index;
    int layer_index;
    PowerSpectrum iso;
    PowerSpectrum aniso;
    PowerSpectrum last;
    PowerSpectrum current;
  };

  void MaybeInit(Angle sun_zenith) const;
  void ComputeSingleScatter(Angle sun_zenith) const;
  void ComputeMultipleScatter(Angle sun_zenith, bool double_scatter) const;
  void ComputeMultipleScatter(Angle sun_zenith, bool double_scatter,
      unsigned int thread_id) const;
  void InterpolateMultipleScatter(bool double_scatter) const;
  void AccumulateMultipleScatter() const;
  void ComputeSkyDome(Angle sun_zenith, ScatteringType scattering_type) const;

  DimensionlessSpectrum GetTransmittance(const Position& p,
      const Position& q) const;
  DimensionlessSpectrum GetTransmittanceSimple(const Position& p,
      const Position& q) const;

  static Position GetRayIntersectionWithLastLayer(const Position& p,
      Angle view_zenith);
  static Length GetLayerHeight(int layer_index);
  static int GetLayerIndex(const Position &p);

  ScatteringType scattering_type_;
  Number rayleigh_density_[kNumLayers];
  Number mie_density_[kNumLayers];
  mutable Angle current_sun_zenith_;
  mutable std::vector<Cell> cells_[kNumTheta][kNumPhi / 2];
  mutable dimensional::BinaryFunction<kNumTheta, kNumPhi / 2, RadianceSpectrum>
      sky_dome_;
};

#endif  // ATMOSPHERE_MODEL_HABER_HABER_H_
