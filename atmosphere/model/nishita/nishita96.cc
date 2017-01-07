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
#include "atmosphere/model/nishita/nishita96.h"

#include <algorithm>
#include <sstream>
#include <string>

namespace {

using dimensional::vec3;

constexpr Length kHorizonDist =
    sqrt(AtmosphereRadius * AtmosphereRadius - EarthRadius * EarthRadius);

}  // anonymous namespace

Nishita96::Nishita96(ScatteringType scattering_type)
    : Nishita93(), scattering_type_(scattering_type) {
}

RadianceSpectrum Nishita96::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  MaybePrecomputeSingleScatteringTables(sun_zenith);
  Position p(0.0 * m, 0.0 * m, EarthRadius + altitude);
  Direction view_dir(cos(view_sun_azimuth) * sin(view_zenith),
      sin(view_sun_azimuth) * sin(view_zenith), cos(view_zenith));

  // Integral of the Rayleigh and Mie phase functions over the solid angle
  // associated with each sampling direction (the solid angle for a sampling
  // direction d_i is formed by the unit vectors v for which v is nearer to d_i
  // than to any other sampling direction d_j, j!=i).
  Number rayleigh_phase_integral[8];
  Number mie_phase_integral[8];
  for (int i = 0; i < 8; ++i) {
    rayleigh_phase_integral[i] = 0.0;
    mie_phase_integral[i] = 0.0;
  }

  constexpr int kNumOmega = 32;
  constexpr Angle dtheta = pi / kNumOmega;
  constexpr Angle dphi = pi / kNumOmega;
  for (int i = 0; i < kNumOmega; ++i) {
    Angle theta = (i + 0.5) * dtheta;
    Number cos_theta = cos(theta);
    for (int j = 0; j < 2 * kNumOmega; ++j) {
      Angle phi = (j + 0.5) * dphi;
      SolidAngle dw = sin(theta) * dtheta.to(rad) * dphi.to(rad) * sr;
      Direction w(cos(phi) * sin(theta), sin(phi) * sin(theta), cos_theta);
      Number drayleigh = RayleighPhaseFunction(dot(view_dir, w)) * dw;
      Number dmie = MiePhaseFunction(dot(view_dir, w)) * dw;
      int nearest_index = 0;
      Number max_dot_product = dot(w, sample_directions_[0]);
      for (int k = 1; k < 8; ++k) {
        Number dot_product = dot(w, sample_directions_[i]);
        if (dot_product > max_dot_product) {
          nearest_index = k;
          max_dot_product = dot_product;
        }
      }
      rayleigh_phase_integral[nearest_index] += drayleigh;
      mie_phase_integral[nearest_index] += dmie;
    }
  }

  Length r = EarthRadius + altitude;
  Number mu_s = cos(sun_zenith);
  Number mu = cos(view_zenith);
  Number nu =
      cos(view_sun_azimuth) * sin(view_zenith) * sin(sun_zenith) + mu * mu_s;
  Length rmu = r * mu;
  Length rmu_s = r * mu_s;

  WavelengthFunction<1, 0, 0, 0, 0> rayleigh_integral(0.0 * m);
  WavelengthFunction<1, 0, 0, 0, 0> mie_integral(0.0 * m);
  WavelengthFunction<0, 0, -1, 0, 0> double_scattering_integral(0.0 / sr);
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

    if (scattering_type_ != DOUBLE_SCATTERING_ONLY) {
      Length rayleigh_sun_length;
      Length mie_sun_length;
      Number mu_s_i = (rmu_s + distance_to_sphere * nu) / r_i;
      Nishita93::GetOpticalLengths(
          r_i, mu_s_i, &rayleigh_sun_length, &mie_sun_length);

      DimensionlessSpectrum optical_depth =
          RayleighScattering() * (rayleigh_length + rayleigh_sun_length) +
          MieExtinction() * (mie_length + mie_sun_length);
      DimensionlessSpectrum transmittance(exp(-optical_depth));
      DimensionlessSpectrum rayleigh_sample(transmittance * rayleigh_density);
      DimensionlessSpectrum mie_sample(transmittance * mie_density);
      rayleigh_integral +=
          (rayleigh_sample + previous_rayleigh_sample) * half_segment_length;
      mie_integral += (mie_sample + previous_mie_sample) * half_segment_length;
      previous_rayleigh_sample = rayleigh_sample;
      previous_mie_sample = mie_sample;
    }

    if (scattering_type_ != SINGLE_SCATTERING_ONLY) {
      Position p_i = p + view_dir * distance_to_sphere;
      // Transmittance from the viewer to p_i.
      DimensionlessSpectrum transmittance_i(exp(-(RayleighScattering() *
          rayleigh_length + MieExtinction() * mie_length)));
      WavelengthFunction<-1, 0, -1, 0, 0> double_scattering(0.0 / m / sr);
      for (int j = 0; j < 8; ++j) {
        double_scattering +=
            sample_single_scattering_[j].GetSingleScattering(p_i) * (
            RayleighScattering() * rayleigh_density * rayleigh_phase_integral[j]
            + MieScattering() * mie_density * mie_phase_integral[j]);
      }
      double_scattering_integral +=
          transmittance_i * double_scattering * (2.0 * half_segment_length);
    }

    previous_rayleigh_density = rayleigh_density;
    previous_mie_density = mie_density;
    distance_to_previous_sphere = distance_to_sphere;
  }

  InverseSolidAngle rayleigh_phase = RayleighPhaseFunction(nu);
  InverseSolidAngle mie_phase = MiePhaseFunction(nu);
  return (rayleigh_integral * RayleighScattering() * rayleigh_phase +
      mie_integral * MieScattering() * mie_phase +
      double_scattering_integral) * SolarSpectrum();
}

void Nishita96::SingleScatteringTable::Init(
    Angle sun_zenith, Angle view_zenith) {
  view_zenith_ = view_zenith;
  cos_view_zenith_ = cos(view_zenith);
  sin_view_zenith_ = sin(view_zenith);
  value_.Set(WavelengthFunction<0, 0, -1, 0, 0>(0.0 / sr));
}

WavelengthFunction<0, 0, -1, 0, 0>
Nishita96::SingleScatteringTable::GetSingleScattering(const Position& p) const {
  vec3 u;
  GetU(p, &u);
  return value_(u.x(), u.y(), u.z());
}

void Nishita96::SingleScatteringTable::GetPosition(const vec3& u,
    Position* p) const {
  Number z =
      view_zenith_ > -pi / 2 && view_zenith_ <= pi / 2.0 ? 1.0 - u.z : u.z;
  if (view_zenith_ == pi / 2.0 || view_zenith_ == -pi / 2.0) {
    p->x = (2.0 * z - 1.0) * kHorizonDist;
    p->y = (2.0 * u.y - 1.0) * kHorizonDist;
    p->z = lerp(EarthRadius, AtmosphereRadius, u.x());
  } else {
    Number tan_view_zenith = sin_view_zenith_ / cos_view_zenith_;
    Length altitude = z * (AtmosphereRadius - EarthRadius);
    p->x = (2.0 * u.x - 1.0) * kHorizonDist + tan_view_zenith * altitude;
    p->y = (2.0 * u.y - 1.0) * kHorizonDist;
    p->z = EarthRadius + altitude;
  }
}

void Nishita96::SingleScatteringTable::GetU(const Position& p, vec3* u) const {
  double z;
  if (view_zenith_ == pi / 2.0 || view_zenith_ == -pi / 2.0) {
    u->x = ((p.z - EarthRadius) / (AtmosphereRadius - EarthRadius))();
    u->y = ((p.y / kHorizonDist + 1.0) / 2.0)();
    z = ((p.x / kHorizonDist + 1.0) / 2.0)();
  } else {
    Number tan_view_zenith = sin_view_zenith_ / cos_view_zenith_;
    Length altitude = p.z - EarthRadius;
    u->x =
        (((p.x - tan_view_zenith * altitude) / kHorizonDist + 1.0) / 2.0)();
    u->y = ((p.y / kHorizonDist + 1.0) / 2.0)();
    z = (altitude / (AtmosphereRadius - EarthRadius))();
  }
  u->z = view_zenith_ > -pi / 2 && view_zenith_ <= pi / 2.0 ? 1.0 - z : z;
  u->x = 0.5 / kNumSteps + (kNumSteps - 1.0) / kNumSteps * u->x;
  u->y = 0.5 / kNumSteps + (kNumSteps - 1.0) / kNumSteps * u->y;
  u->z = 0.5 / kNumSteps + (kNumSteps - 1.0) / kNumSteps * u->z;
}

void Nishita96::MaybePrecomputeSingleScatteringTables(Angle sun_zenith) const {
  if (current_sun_zenith_ == sun_zenith) {
    return;
  }
  current_sun_zenith_ = sun_zenith;
  for (int i = 0; i < 8; ++i) {
    Angle view_zenith;
    switch (i) {
      case 0:
        view_zenith = 0.0 * rad;
        break;
      case 1:
        view_zenith = pi;
        break;
      case 2:
        view_zenith = pi / 2.0;
        break;
      case 3:
        view_zenith = -pi / 2.0;
        break;
      case 4:
        view_zenith = sun_zenith;
        break;
      case 5:
        view_zenith = sun_zenith + pi;
        break;
      case 6:
        view_zenith = sun_zenith - pi / 2.0;
        break;
      case 7:
        view_zenith = sun_zenith + pi / 2.0;
        break;
    }
    sample_directions_[i] = Direction(sin(view_zenith), 0.0, cos(view_zenith));
    sample_single_scattering_[i].Init(sun_zenith, view_zenith);

    std::stringstream filename;
    filename << "output/cache/nishita/single_scatter_" << sun_zenith.to(deg)
        << "_" << i << ".dat";
    std::ifstream f;
    f.open(filename.str());
    if (f.good()) {
      f.close();
      if (i == 0) {
        std::cout << "Loading single scattering tables: "
                  << sun_zenith.to(deg) << std::endl;
      }
      sample_single_scattering_[i].value_.Load(filename.str());
    } else {
      if (i == 0) {
        std::cout << "Precompute single scattering tables: "
                  << sun_zenith.to(deg) << std::endl;
      }
      PreComputeSingleScatteringTable(
          sun_zenith, view_zenith, &sample_single_scattering_[i]);
      sample_single_scattering_[i].value_.Save(filename.str());
    }
  }
}

void Nishita96::PreComputeSingleScatteringTable(Angle sun_zenith,
    Angle view_zenith, SingleScatteringTable* output) const {
  Direction sun_dir(sin(sun_zenith), 0.0, cos(sun_zenith));
  for (int i = 0; i < kNumSteps; ++i) {
    for (int j = 0; j < kNumSteps; ++j) {
      Position p0;
      Position p1;
      output->GetPosition(vec3(i, j, 0) / Number(kNumSteps - 1.0), &p0);
      output->GetPosition(vec3(i, j, 1) / Number(kNumSteps - 1.0), &p1);
      Length step = length(p1 - p0);
      Direction d = (p1 - p0) / step;
      PreComputeSingleScattering(p0, d, step, sun_dir, i, j, &(output->value_));
    }
  }
}

void Nishita96::PreComputeSingleScattering(const Position& p,
    const Direction& d, Length step, const Direction& sun_dir, int i, int j,
    SingleScatteringFunction* output) const {
  Length start_t;
  Length end_t;
  if (dot(p, p) >= AtmosphereRadius * AtmosphereRadius) {
    start_t = DistanceToTopOfAtmosphere(p, d);
    if (start_t < 0.0 * m) {
      return;
    }
    end_t = DistanceToGroundOrTopOfAtmosphere(p + d * start_t, d) + start_t;
  } else {
    start_t = -DistanceToGroundOrTopOfAtmosphere(p, -d);
    end_t = DistanceToGroundOrTopOfAtmosphere(p, d);
  }
  assert(end_t >= start_t);

  // Integral of the inscattering from the current point to the top of the
  // atmosphere or to the ground, in view direction -d (without the constant
  // phase function, solar spectrum and RayleighScattering() or MieScattering()
  // factors, added later). This integral is computed from point to point in
  // direction d, where the integral at a new point is the integral at the
  // previous point times the transmittance between these two points, plus the
  // integral over the segment between the two points.
  WavelengthFunction<1, 0, 0, 0, 0> rayleigh_integral(0.0 * m);
  WavelengthFunction<1, 0, 0, 0, 0> mie_integral(0.0 * m);

  Length altitude = length(p + d * start_t) - EarthRadius;
  Number previous_rayleigh_density = exp(-altitude / RayleighScaleHeight);
  Number previous_mie_density = exp(-altitude / MieScaleHeight);
  DimensionlessSpectrum previous_sun_transmittance;
  GetTransmittance(p + d * start_t, sun_dir, &previous_sun_transmittance);

  // Background radiance, i.e. either 0 if the start point is the top of the
  // atmosphere, or the ground radiance (Sun irradiance times the transmittance
  // times the albedo times the Lambert BRDF). Does not include the solar
  // spectrum factor, added later.
  WavelengthFunction<0, 0, -1, 0, 0> background(0.0 / sr);
  if (altitude < (AtmosphereRadius - EarthRadius) * 0.5) {
    background = previous_sun_transmittance *
        (GroundAlbedo() * (1.0 / (PI * sr)) *
        std::max(0.0, dot(normalize(p + d * start_t), sun_dir)()));
  }

  Number nu = dot(-d, sun_dir);
  const auto rayleigh_factor = RayleighScattering() * RayleighPhaseFunction(nu);
  const auto mie_factor = MieScattering() * MiePhaseFunction(nu);

  int start_k = static_cast<int>(ceil((start_t / step)()));
  int end_k = static_cast<int>(ceil((end_t / step)()));
  assert(end_k >= start_k);
  for (int k = start_k; k <= end_k; ++k) {
    Position p_k = p + d * (k * step);
    Length altitude_k = length(p_k) - EarthRadius;
    Number rayleigh_density = exp(-altitude_k / RayleighScaleHeight);
    Number mie_density = exp(-altitude_k / MieScaleHeight);
    Length half_segment_length =
        (k == start_k ? k * step - start_t : step) * 0.5;
    Length rayleigh_length =
        (rayleigh_density + previous_rayleigh_density) * half_segment_length;
    Length mie_length =
        (mie_density + previous_mie_density) * half_segment_length;
    DimensionlessSpectrum segment_optical_depth =
        RayleighScattering() * rayleigh_length + MieExtinction() * mie_length;
    DimensionlessSpectrum segment_transmittance(exp(-segment_optical_depth));

    DimensionlessSpectrum sun_transmittance;
    GetTransmittance(p_k, sun_dir, &sun_transmittance);

    // Integral of the inscattering over the current segment (without the
    // constant phase function, solar spectrum and RayleighScattering() or
    // MieScattering() factors, added later). Integral computed using
    // trapezoidal integration.
    WavelengthFunction<1, 0, 0, 0, 0> rayleigh_segment_integral =
        (sun_transmittance * rayleigh_density + previous_sun_transmittance *
         segment_transmittance * previous_rayleigh_density) *
         half_segment_length;
    WavelengthFunction<1, 0, 0, 0, 0> mie_segment_integral =
        (sun_transmittance * mie_density + previous_sun_transmittance *
         segment_transmittance * previous_mie_density) * half_segment_length;

    // Incrementally update the inscattering integral value as described above.
    rayleigh_integral =
        rayleigh_segment_integral + segment_transmittance * rayleigh_integral;
    mie_integral = mie_segment_integral + segment_transmittance * mie_integral;
    // Incrementally update the background radiance.
    background = background * segment_transmittance;

    if (k >= 0 && k < kNumSteps) {
      output->Set(i, j, k, rayleigh_integral * rayleigh_factor +
          mie_integral * mie_factor + background);
    }
    previous_rayleigh_density = rayleigh_density;
    previous_mie_density = mie_density;
    previous_sun_transmittance = sun_transmittance;
  }
}

// p must be inside the atmosphere and d must be normalized.
void Nishita96::GetTransmittance(const Position& p, const Direction& d,
     DimensionlessSpectrum* transmittance) const {
  if (IntersectsGround(p, d)) {
    *transmittance = DimensionlessSpectrum(0.0);
  } else {
    Length rayleigh_length;
    Length mie_length;
    Length r = length(p);
    Number mu = dot(p, d) / r;
    GetOpticalLengths(r, mu, &rayleigh_length, &mie_length);
    DimensionlessSpectrum optical_depth =
        RayleighScattering() * rayleigh_length + MieExtinction() * mie_length;
    *transmittance = exp(-optical_depth);
  }
}

bool Nishita96::IntersectsGround(const Position& p, const Direction& d) {
  Length b = dot(p, d);
  Area c = dot(p, p) - EarthRadius * EarthRadius;
  Area delta_sq = b * b - c;
  return delta_sq >= 0.0 * m2 && -b - sqrt(delta_sq) > 0.0 * m;
}

Length Nishita96::DistanceToSphere(Length r, Length rmu, Length sphere_radius) {
  Area delta_sq = sphere_radius * sphere_radius - r * r + rmu * rmu;
  return delta_sq < 0.0 * m2 ? 0.0 * m :
      (r < sphere_radius ? -rmu + sqrt(delta_sq) : -rmu - sqrt(delta_sq));
}

Length Nishita96::DistanceToTopOfAtmosphere(const Position& p,
    const Direction& d) {
  Length b = dot(p, d);
  Area c = dot(p, p) - AtmosphereRadius * AtmosphereRadius;
  Area delta_sq = b * b - c;
  return delta_sq >= 0.0 * m2 ? -b - sqrt(delta_sq) : -1.0 * m;
}

Length Nishita96::DistanceToGroundOrTopOfAtmosphere(const Position& p,
    const Direction& d) {
  Length b = dot(p, d);
  Area c = dot(p, p) - AtmosphereRadius * AtmosphereRadius;
  Area delta_sq = b * b - c;
  assert(delta_sq >= 0.0 * m2);
  Length distance_to_top_of_atmosphere = -b + sqrt(delta_sq);

  delta_sq += EarthRadius * EarthRadius - AtmosphereRadius * AtmosphereRadius;
  if (delta_sq >= 0.0 * m2) {
    Length distance_to_ground = -b - sqrt(delta_sq);
    if (distance_to_ground >= 0.0 * m) {
      return std::min(distance_to_top_of_atmosphere, distance_to_ground);
    }
  }
  return distance_to_top_of_atmosphere;
}
