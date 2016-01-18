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
#include "atmosphere/model/nishita/nishita93.h"

#include <algorithm>

Nishita93::Nishita93() {
  // Precomputes the optical length lookup tables.
  for (int i = 0; i < kNumSphere; ++i) {
    Length r = GetSphereRadius(i);
    Length h = r - EarthRadius;
    for (int j = 0; j < kNumCylinder; ++j) {
      Length c = std::min(GetCylinderRadius(j), r);
      Length rmu = sqrt(r * r - c * c);

      Length rayleigh_length = 0.0 * m;
      Length mie_length = 0.0 * m;
      Number previous_rayleigh_density = exp(-h / RayleighScaleHeight);
      Number previous_mie_density = exp(-h / MieScaleHeight);
      Length distance_to_previous_sphere = 0.0 * m;
      for (int k = i + 1; k < kNumSphere; ++k) {
        Length r_k = GetSphereRadius(k);
        Length h_k = r_k - EarthRadius;
        Length distance_to_sphere = DistanceToSphere(r, rmu, r_k);
        Number rayleigh_density = exp(-h_k / RayleighScaleHeight);
        Number mie_density = exp(-h_k / MieScaleHeight);
        Length segment_length =
            distance_to_sphere - distance_to_previous_sphere;
        rayleigh_length += (rayleigh_density + previous_rayleigh_density) / 2 *
            segment_length;
        mie_length += (mie_density + previous_mie_density) / 2 * segment_length;
        previous_rayleigh_density = rayleigh_density;
        previous_mie_density = mie_density;
        distance_to_previous_sphere = distance_to_sphere;
      }
      rayleigh_optical_length_.Set(i, j, rayleigh_length);
      mie_optical_length_.Set(i, j, mie_length);

      rmu = -rmu;
      rayleigh_length = 0.0 * m;
      mie_length = 0.0 * m;
      previous_rayleigh_density = exp(-h / RayleighScaleHeight);
      previous_mie_density = exp(-h / MieScaleHeight);
      distance_to_previous_sphere = 0.0 * m;
      for (int k = i - 1; k > -kNumSphere; --k) {
        Length r_k = GetSphereRadius(std::abs(k));
        Length h_k = r_k - EarthRadius;
        Length distance_to_sphere = DistanceToSphere(r, rmu, r_k);
        if (distance_to_sphere == 0.0 * m) {
          continue;
        }
        Number rayleigh_density = exp(-h_k / RayleighScaleHeight);
        Number mie_density = exp(-h_k / MieScaleHeight);
        Length segment_length =
            distance_to_sphere - distance_to_previous_sphere;
        rayleigh_length += (rayleigh_density + previous_rayleigh_density) / 2 *
            segment_length;
        mie_length += (mie_density + previous_mie_density) / 2 * segment_length;
        previous_rayleigh_density = rayleigh_density;
        previous_mie_density = mie_density;
        distance_to_previous_sphere = distance_to_sphere;
      }
      rayleigh_opposite_optical_length_.Set(i, j, rayleigh_length);
      mie_opposite_optical_length_.Set(i, j, mie_length);
    }
  }
}

IrradianceSpectrum Nishita93::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  Length rayleigh_length;
  Length mie_length;
  GetOpticalLengths(EarthRadius + altitude, cos(sun_zenith), &rayleigh_length,
      &mie_length);
  DimensionlessSpectrum optical_depth = RayleighScattering() * rayleigh_length +
      MieExtinction() * mie_length;
  DimensionlessSpectrum transmittance(exp(-optical_depth));
  return transmittance * SolarSpectrum();
}

RadianceSpectrum Nishita93::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  Length r = EarthRadius + altitude;
  Number mu_s = cos(sun_zenith);
  Number mu = cos(view_zenith);
  Number nu =
      cos(view_sun_azimuth) * sin(view_zenith) * sin(sun_zenith) + mu * mu_s;
  Length rmu = r * mu;
  Length rmu_s = r * mu_s;

  WavelengthFunction<Length> rayleigh_integral(0.0 * m);
  WavelengthFunction<Length> mie_integral(0.0 * m);
  Length rayleigh_length = 0.0 * m;
  Length mie_length = 0.0 * m;
  Number previous_rayleigh_density = exp(-altitude / RayleighScaleHeight);
  Number previous_mie_density = exp(-altitude / MieScaleHeight);
  Length distance_to_previous_sphere = 0.0 * m;
  DimensionlessSpectrum previous_rayleigh_sample(0.0);
  DimensionlessSpectrum previous_mie_sample(0.0);
  for (int i = 0; i < kNumSphere; ++i) {
    Length r_i = GetSphereRadius(i);
    if (r_i <= r) {
      continue;
    }
    Length h_i = r_i - EarthRadius;
    Number rayleigh_density = exp(-h_i / RayleighScaleHeight);
    Number mie_density = exp(-h_i / MieScaleHeight);
    Length distance_to_sphere = DistanceToSphere(r, rmu, r_i);
    Length half_segment_length =
        (distance_to_sphere - distance_to_previous_sphere) * 0.5;
    rayleigh_length +=
        (rayleigh_density + previous_rayleigh_density) * half_segment_length;
    mie_length += (mie_density + previous_mie_density) * half_segment_length;

    Length rayleigh_sun_length;
    Length mie_sun_length;
    Number mu_s_i = (rmu_s + distance_to_sphere * nu) / r_i;
    GetOpticalLengths(r_i, mu_s_i, &rayleigh_sun_length, &mie_sun_length);

    DimensionlessSpectrum optical_depth =
        RayleighScattering() * (rayleigh_length + rayleigh_sun_length) +
        MieExtinction() * (mie_length + mie_sun_length);
    DimensionlessSpectrum transmittance(exp(-optical_depth));
    DimensionlessSpectrum rayleigh_sample(transmittance * rayleigh_density);
    DimensionlessSpectrum mie_sample(transmittance * mie_density);
    rayleigh_integral +=
        (rayleigh_sample + previous_rayleigh_sample) * half_segment_length;
    mie_integral += (mie_sample + previous_mie_sample) * half_segment_length;

    previous_rayleigh_density = rayleigh_density;
    previous_mie_density = mie_density;
    distance_to_previous_sphere = distance_to_sphere;
    previous_rayleigh_sample = rayleigh_sample;
    previous_mie_sample = mie_sample;
  }

  InverseSolidAngle rayleigh_phase = RayleighPhaseFunction(nu);
  InverseSolidAngle mie_phase = MiePhaseFunction(nu);
  return (rayleigh_integral * RayleighScattering() * rayleigh_phase +
      mie_integral * MieScattering() * mie_phase) * SolarSpectrum();
}

void Nishita93::GetOpticalLengths(Length r, Number mu, Length* rayleigh_length,
    Length* mie_length) const {
  constexpr Number a =
      exp(-(AtmosphereRadius - EarthRadius) / RayleighScaleHeight) - 1.0;
  Number x = (exp(-(r - EarthRadius) / RayleighScaleHeight) - 1.0) / a;
  Number y = r * sqrt(1.0 - mu * mu) / AtmosphereRadius;
  x = 0.5 / kNumSphere + (kNumSphere - 1.0) / kNumSphere * x;
  y = 0.5 / kNumCylinder + (kNumCylinder - 1.0) / kNumCylinder * y;
  if (mu >= 0.0) {
    *rayleigh_length = rayleigh_optical_length_(x(), y());
    *mie_length = mie_optical_length_(x(), y());
  } else {
    *rayleigh_length = rayleigh_opposite_optical_length_(x(), y());
    *mie_length = mie_opposite_optical_length_(x(), y());
  }
}

Length Nishita93::GetSphereRadius(int sphere_index) {
  constexpr Number a =
      exp(-(AtmosphereRadius - EarthRadius) / RayleighScaleHeight) - 1.0;
  double x = sphere_index / static_cast<double>(kNumSphere - 1);
  return EarthRadius - RayleighScaleHeight * log(a * x + 1.0);
}

Length Nishita93::GetCylinderRadius(int cylinder_index) {
  double x = cylinder_index / static_cast<double>(kNumCylinder - 1);
  return AtmosphereRadius * x;
}

Length Nishita93::DistanceToSphere(Length r, Length rmu, Length sphere_radius) {
  Area delta_sq = sphere_radius * sphere_radius - r * r + rmu * rmu;
  return delta_sq < 0.0 * m2 ? 0.0 * m :
      (r < sphere_radius ? -rmu + sqrt(delta_sq) : -rmu - sqrt(delta_sq));
}
