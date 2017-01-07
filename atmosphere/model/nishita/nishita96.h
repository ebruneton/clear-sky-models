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
#ifndef ATMOSPHERE_MODEL_NISHITA_NISHITA96_H_
#define ATMOSPHERE_MODEL_NISHITA_NISHITA96_H_

#include "atmosphere/atmosphere.h"
#include "atmosphere/model/nishita/nishita93.h"
#include "math/angle.h"
#include "math/ternary_function.h"
#include "math/vector.h"
#include "physics/units.h"

class Nishita96 : public Nishita93 {
 public:
  enum ScatteringType {
    SINGLE_SCATTERING_ONLY, DOUBLE_SCATTERING_ONLY, ALL_ORDERS
  };

  explicit Nishita96(ScatteringType scattering_type);

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun_azimuth) const override;

 private:
  static constexpr int kNumSteps = 33;

  typedef dimensional::TernaryFunction<kNumSteps, kNumSteps, kNumSteps,
      WavelengthFunction<0, 0, -1, 0, 0>> SingleScatteringFunction;
  typedef dimensional::Vector3<Length> Position;
  typedef dimensional::Vector3<Number> Direction;

  class SingleScatteringTable {
   public:
    void Init(Angle sun_zenith, Angle view_zenith);

    WavelengthFunction<0, 0, -1, 0, 0> GetSingleScattering(
        const Position& p) const;

    void GetPosition(const dimensional::vec3& u, Position* p) const;
    void GetU(const Position& p, dimensional::vec3* u) const;

    Angle view_zenith_;
    Number cos_view_zenith_;
    Number sin_view_zenith_;
    SingleScatteringFunction value_;
  };

  void MaybePrecomputeSingleScatteringTables(Angle sun_zenith) const;

  void PreComputeSingleScatteringTable(Angle sun_zenith, Angle view_zenith,
      SingleScatteringTable* output) const;

  void PreComputeSingleScattering(const Position& p, const Direction& d,
      Length step, const Direction& sun_dir, int i, int j,
      SingleScatteringFunction* output) const;

  void GetTransmittance(const Position& p, const Direction& d,
     DimensionlessSpectrum* transmittance) const;

  // d must be normalized.
  static bool IntersectsGround(const Position& p, const Direction& d);

  // Returns the distance from a point at radius r to the sphere of radius
  // sphere_radius in a direction whose angle with the local vertical is
  // acos(rmu / r), or 0 if there is no intersection.
  static Length DistanceToSphere(Length r, Length rmu, Length sphere_radius);

  // p must be outside the atmosphere and d must be normalized.
  static Length DistanceToTopOfAtmosphere(const Position& p,
      const Direction& d);

  // p must be inside the atmosphere and d must be normalized.
  static Length DistanceToGroundOrTopOfAtmosphere(const Position& p,
      const Direction& d);

  ScatteringType scattering_type_;
  mutable Angle current_sun_zenith_;
  mutable dimensional::Vector3<Number> sample_directions_[8];
  mutable SingleScatteringTable sample_single_scattering_[8];
};

#endif  // ATMOSPHERE_MODEL_NISHITA_NISHITA96_H_
