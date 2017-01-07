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

/*
 Adapted from Sean O'Neil "Accurate Atmospheric Scattering", GPU Gems 2, 2004.
 More precisely from the SkyFromAtmosphere.vert file from the GPU Gems 2 CD.

 Original header:
    // Atmospheric scattering vertex shader
    //
    // Author: Sean O'Neil
    //
    // Copyright (c) 2004 Sean O'Neil
    //
*/

#include "atmosphere/model/oneal/oneal.h"

#include "math/vector.h"

typedef dimensional::Vector3<Length> Position;
typedef dimensional::Vector3<Length> Vector;
typedef dimensional::Vector3<Number> Direction;

namespace {

// Returns the distance from a point at radius r to the sphere of radius
// sphere_radius in a direction whose angle with the local vertical is
// acos(rmu / r), or 0 if there is no intersection.
Length DistanceToSphere(Length r, Length rmu, Length sphere_radius) {
  Area delta_sq = sphere_radius * sphere_radius - r * r + rmu * rmu;
  return delta_sq < 0.0 * m2 ? 0.0 * m :
      (r < sphere_radius ? -rmu + sqrt(delta_sq) : -rmu - sqrt(delta_sq));
}

// Scale function based on "An Approximation to the Chapman Grazing-Incidence
// Function for Atmospheric Scattering", C. Schüler, GPU Pro 3, 2012. We use
// this to replace the original function from O'Neal, which is specific to the
// atmospheric parameters he used.
Length Scale(Length H, Number mu) {
  Number c = sqrt(EarthRadius / H);
  return H * c / (c * mu + 1.0);
}

}  // namespace

ONeal::ONeal() {
}

IrradianceSpectrum ONeal::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  Number f_rayleigh_depth = exp(-altitude / RayleighScaleHeight);
  Number f_mie_depth = exp(-altitude / MieScaleHeight);
  Number f_light_angle = cos(sun_zenith);
  Length f_rayleigh_scatter =
      f_rayleigh_depth * Scale(RayleighScaleHeight, f_light_angle);
  Length f_mie_scatter = f_mie_depth * Scale(MieScaleHeight, f_light_angle);
  DimensionlessSpectrum attenuate = exp(-(MieExtinction() * f_mie_scatter +
      RayleighScattering() * f_rayleigh_scatter));
  return attenuate * SolarSpectrum();
}

RadianceSpectrum ONeal::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  Length f_inner_radius = EarthRadius;
  Length f_camera_height = f_inner_radius + altitude;
  Position v3_camera_pos(0.0 * m, 0.0 * m, f_camera_height);
  Direction v3_light_pos(cos(view_sun_azimuth) * sin(sun_zenith),
      sin(view_sun_azimuth) * sin(sun_zenith), cos(sun_zenith));

  // Get the ray from the camera, and its length.
  Direction v3_ray(sin(view_zenith), 0.0, cos(view_zenith));
  Length f_far = DistanceToSphere(
      f_camera_height, f_camera_height * cos(view_zenith), AtmosphereRadius);

  // Calculate the ray's starting position, and calculate its scattering offset.
  Position v3_start = v3_camera_pos;
  Length f_height = length(v3_start);
  Number f_rayleigh_depth =
      exp((f_inner_radius - f_height) / RayleighScaleHeight);
  Number f_mie_depth = exp((f_inner_radius - f_height) / MieScaleHeight);
  Number f_start_angle = dot(v3_ray, v3_start) / f_height;
  Length f_rayleigh_start_offset =
      f_rayleigh_depth * Scale(RayleighScaleHeight, f_start_angle);
  Length f_mie_start_offset =
      f_mie_depth * Scale(MieScaleHeight, f_start_angle);

  // Note: the original implementation uses only 2 samples.
  const int kNumSamples = 4;

  // Initialize the scattering loop variables.
  Length f_sample_length = f_far / kNumSamples;
  Vector v3_sample_ray = v3_ray * f_sample_length;
  Position v3_sample_point = v3_start + v3_sample_ray * Number(0.5);

  // Now loop through the sample rays.
  WavelengthFunction<1, 0, 0, 0, 0> front_color(0.0 * m);
  WavelengthFunction<1, 0, 0, 0, 0> front_secondary_color(0.0 * m);
  for (int i = 0; i < kNumSamples; ++i) {
    f_height = length(v3_sample_point);
    f_rayleigh_depth = exp((f_inner_radius - f_height) / RayleighScaleHeight);
    f_mie_depth = exp((f_inner_radius - f_height) / MieScaleHeight);
    Number f_light_angle = dot(v3_light_pos, v3_sample_point) / f_height;
    Number f_camera_angle = dot(v3_ray, v3_sample_point) / f_height;
    Length f_rayleigh_scatter = f_rayleigh_start_offset + f_rayleigh_depth * (
        Scale(RayleighScaleHeight, f_light_angle) -
        Scale(RayleighScaleHeight, f_camera_angle));
    Length f_mie_scatter = f_mie_start_offset + f_mie_depth * (
        Scale(MieScaleHeight, f_light_angle) -
        Scale(MieScaleHeight, f_camera_angle));
    DimensionlessSpectrum attenuate = exp(-(MieExtinction() * f_mie_scatter +
        RayleighScattering() * f_rayleigh_scatter));
    front_color = front_color +
        attenuate * (f_rayleigh_depth * f_sample_length);
    front_secondary_color = front_secondary_color +
        attenuate * (f_mie_depth * f_sample_length);
    v3_sample_point = v3_sample_point + v3_sample_ray;
  }

  Number nu = dot(v3_ray, v3_light_pos);
  InverseSolidAngle rayleigh_phase = RayleighPhaseFunction(nu);
  InverseSolidAngle mie_phase = MiePhaseFunction(nu);
  return (front_color * RayleighScattering() * rayleigh_phase +
      front_secondary_color * MieScattering() * mie_phase) * SolarSpectrum();
}
