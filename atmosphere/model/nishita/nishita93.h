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
#ifndef ATMOSPHERE_MODEL_NISHITA_NISHITA93_H_
#define ATMOSPHERE_MODEL_NISHITA_NISHITA93_H_

#include "atmosphere/atmosphere.h"
#include "math/angle.h"
#include "math/binary_function.h"
#include "physics/units.h"

class Nishita93 : public Atmosphere {
 public:
  Nishita93();

  int GetOriginalNumberOfWavelengths() const override { return 3; }

  IrradianceSpectrum GetSunIrradiance(Length altitude,
      Angle sun_zenith) const override;

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun_azimuth) const override;

 protected:
  // Returns the optical lengths for rayleigh and mie particles between a point
  // at radius r and the top of the atmosphere in a direction whose angle with
  // the local vertical is acos(mu). This is done by using the precomputed
  // optical length lookup tables.
  void GetOpticalLengths(Length r, Number mu, Length* rayleigh_length,
      Length* mie_length) const;

  // Computes the sphere radius noted r_i in Nishita93.
  static Length GetSphereRadius(int sphere_index);

  // Computes the cylinder radius noted C_j in Nishita93.
  static Length GetCylinderRadius(int cylinder_index);

  // Returns the distance from a point at radius r to the sphere of radius
  // sphere_radius in a direction whose angle with the local vertical is
  // acos(rmu / r), or 0 if there is no intersection.
  static Length DistanceToSphere(Length r, Length rmu, Length sphere_radius);

  // Number of virtual spheres and cylinders for the optical length lookup table
  // described in Section 4.3.4 and Fig. 4 of Nishita93.
  static constexpr int kNumSphere = 64;
  static constexpr int kNumCylinder = 64;

  // Optical lengths between a point of the intersection circle between a sphere
  // and a cylinder, and the Sun. The cylinder axis being the Sun direction.
  dimensional::BinaryFunction<kNumSphere, kNumCylinder, Length>
      rayleigh_optical_length_;
  dimensional::BinaryFunction<kNumSphere, kNumCylinder, Length>
      mie_optical_length_;
  // Same, with the Sun in the opposite direction of the cylinder axis.
  dimensional::BinaryFunction<kNumSphere, kNumCylinder, Length>
      rayleigh_opposite_optical_length_;
  dimensional::BinaryFunction<kNumSphere, kNumCylinder, Length>
      mie_opposite_optical_length_;
};

#endif  // ATMOSPHERE_MODEL_NISHITA_NISHITA93_H_
