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
#include "atmosphere/model/haber/haber.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>

#include "util/progress_bar.h"

namespace {

int Lookup(const std::map<int, int>& m, int key) {
  if (m.find(key) == m.end()) {
    return -1;
  }
  return m.find(key)->second;
}

std::string GetCacheFileName(Angle sun_zenith,
    Haber::ScatteringType scattering_type) {
  std::stringstream name;
  name << "output/cache/haber/sky_dome_" << sun_zenith.to(deg) << "_";
  switch (scattering_type) {
    case Haber::SINGLE_SCATTERING_ONLY:
      name << "single.dat";
      break;
    case Haber::DOUBLE_SCATTERING_ONLY:
      name << "double.dat";
      break;
    default:
      name << "all.dat";
      break;
  }
  return name.str();
}

}  // anonymous namespace

constexpr Angle Haber::kDeltaPhi;
constexpr Number Haber::kShellRatio;
constexpr Length Haber::kMinShellRadius;
constexpr Length Haber::kMaxShellRadius;
constexpr Length Haber::kMaxLayerHeight;

Haber::Haber(ScatteringType scattering_type)
    : scattering_type_(scattering_type) {
  for (int i = 0; i < kNumLayers; ++i) {
    Length layer_center = (GetLayerHeight(i) + GetLayerHeight(i + 1)) * 0.5;
    rayleigh_density_[i] = exp(-layer_center / RayleighScaleHeight);
    mie_density_[i] = exp(-layer_center / MieScaleHeight);
  }
  int num_cells = 0;
  for (int i = 0; i < kNumTheta; ++i) {
    Angle theta_min = i * kDeltaPhi;
    Angle theta_max = (i + 1) * kDeltaPhi;
    Angle theta = (i + 0.5) * kDeltaPhi;
    for (int j = 0; j < kNumPhi / 2; ++j) {
      Angle phi = (j + 0.5) * kDeltaPhi;
      for (int k = 0; k < kNumShell; ++k) {
        Length r_min = kMinShellRadius * pow(kShellRatio, k);
        Length r_max = kMinShellRadius * pow(kShellRatio, k + 1);
        Length r = (r_min + r_max) * 0.5;
        Cell cell;
        cell.center.x = r * cos(phi) * sin(theta);
        cell.center.y = r * sin(phi) * sin(theta);
        cell.center.z = r * cos(theta) + EarthRadius;
        if (length(cell.center) <= EarthRadius) {
          continue;
        }
        if (length(cell.center) >= EarthRadius + kMaxLayerHeight) {
          break;
        }
        cell.volume = (r_max * r_max * r_max - r_min * r_min * r_min) / 3.0 *
            (2.0 * PI / kNumPhi) * (cos(theta_min) - cos(theta_max));
        cell.radial_width = r_max - r_min;
        cell.shell_index = k;
        cell.layer_index = GetLayerIndex(cell.center);
        assert(cell.layer_index >= 0 && cell.layer_index < kNumLayers);
        cells_[i][j].push_back(cell);
        ++num_cells;
      }
    }
  }
  for (int i = 0; i < kNumTheta; ++i) {
    for (int j = 0; j < kNumPhi / 2; ++j) {
      assert(cells_[i][j].size() == cells_[i][0].size());
    }
  }
}

IrradianceSpectrum Haber::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  Position p = Position(0.0 * m, 0.0 * m, EarthRadius + 1.0 * m);
  Position q = GetRayIntersectionWithLastLayer(p, sun_zenith);
  return GetTransmittance(p, q) * SolarSpectrum();
}

RadianceSpectrum Haber::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  assert(altitude == 0.0 * m);
  MaybeInit(sun_zenith);
  if (view_sun_azimuth < 0.0 * deg) {
    view_sun_azimuth = -view_sun_azimuth;
  }
  while (view_sun_azimuth > 2.0 * pi) {
    view_sun_azimuth = view_sun_azimuth - 2.0 * pi;
  }
  if (view_sun_azimuth > pi) {
    view_sun_azimuth = 2.0 * pi - view_sun_azimuth;
  }
  Number u = view_zenith / (kNumTheta * kDeltaPhi);
  Number v = view_sun_azimuth / (kNumPhi / 2 * kDeltaPhi);
  return sky_dome_(u(), v());
}

void Haber::MaybeInit(Angle sun_zenith) const {
  if (sun_zenith == current_sun_zenith_) {
    return;
  }
  current_sun_zenith_ = sun_zenith;

  const std::string name = GetCacheFileName(sun_zenith, scattering_type_);
  std::ifstream f;
  f.open(name);
  if (f.good()) {
    f.close();
  } else {
    ComputeSingleScatter(sun_zenith);
    ComputeSkyDome(sun_zenith, SINGLE_SCATTERING_ONLY);
    sky_dome_.Save(GetCacheFileName(sun_zenith, SINGLE_SCATTERING_ONLY));

    constexpr int kNumScatteringOrders = 4;
    for (int i = 2; i <= kNumScatteringOrders; ++i) {
      std::cout << "Precomputing (sun zenith angle " << sun_zenith.to(deg)
          << "), step " << i - 1 << "/" << kNumScatteringOrders - 1 << "..."
          << std::endl;
      ComputeMultipleScatter(sun_zenith, i == 2);
      InterpolateMultipleScatter(i == 2);
      AccumulateMultipleScatter();
      if (i == 2) {
        ComputeSkyDome(sun_zenith, DOUBLE_SCATTERING_ONLY);
        sky_dome_.Save(GetCacheFileName(sun_zenith, DOUBLE_SCATTERING_ONLY));
      }
    }
    ComputeSkyDome(sun_zenith, ALL_ORDERS);
    sky_dome_.Save(GetCacheFileName(sun_zenith, ALL_ORDERS));
  }

  sky_dome_.Load(name);
}

void Haber::ComputeSingleScatter(Angle sun_zenith) const {
  for (int i = 0; i < kNumTheta; ++i) {
    for (int j = 0; j < kNumPhi / 2; ++j) {
      std::vector<Cell>& cells = cells_[i][j];
      for (unsigned int k = 0; k < cells.size(); ++k) {
        Cell& cell = cells[k];
        Position q =
            GetRayIntersectionWithLastLayer(cell.center, sun_zenith);
        auto I =
            GetTransmittance(cell.center, q) * SolarSpectrum() * cell.volume;
        cell.iso =
            I * RayleighScattering() * rayleigh_density_[cell.layer_index];
        cell.aniso = I * MieScattering() * mie_density_[cell.layer_index];
        cell.current = PowerSpectrum(0.0 * watt / nm);
      }
    }
  }
}

void Haber::ComputeMultipleScatter(Angle sun_zenith,
    bool double_scatter) const {
  int total_progress = 0;
  for (int i_tgt = 0; i_tgt < kNumTheta; ++i_tgt) {
    for (size_t k_tgt = 0; k_tgt < cells_[i_tgt][0].size(); ++k_tgt) {
      if (cells_[i_tgt][0][k_tgt].shell_index % 4 != 0) {
        continue;
      }
      for (int i_src = 0; i_src < kNumTheta; ++i_src) {
        for (size_t k_src = 0; k_src < cells_[i_src][0].size(); ++k_src) {
          if (cells_[i_src][0][k_src].shell_index % 4 != 0) {
            continue;
          }
          total_progress += 1;
        }
      }
    }
  }
  ProgressBar progress_bar(total_progress);

  RunJobs([&](int i_tgt) {
    Direction sun_dir = Direction(sin(sun_zenith), 0.0, cos(sun_zenith));
    DimensionlessSpectrum transmittances[kNumPhi];
    for (size_t k_tgt = 0; k_tgt < cells_[i_tgt][0].size(); ++k_tgt) {
      Cell& target0 = cells_[i_tgt][0][k_tgt];
      if (target0.shell_index % 4 != 0) {
        continue;
      }
      auto target0_scatter = (
          RayleighScattering() * rayleigh_density_[target0.layer_index] +
          MieScattering() * mie_density_[target0.layer_index]) * target0.volume;

      for (int i_src = 0; i_src < kNumTheta; ++i_src) {
        for (size_t k_src = 0; k_src < cells_[i_src][0].size(); ++k_src) {
          if (cells_[i_src][0][k_src].shell_index % 4 != 0) {
            continue;
          }
          // Precompute and cache the transmittances from target0 to all cells
          // in ring (i_src,*,k_src).
          for (int j_src = 0; j_src < kNumPhi; ++j_src) {
            int j = j_src < kNumPhi / 2 ? j_src : kNumPhi - 1 - j_src;
            Cell& src = cells_[i_src][j][k_src];
            if (i_src == i_tgt && j_src == 0 && k_src == k_tgt) {
              transmittances[j_src] = DimensionlessSpectrum(1.0);
            } else {
              Position src_center = src.center;
              if (j_src >= kNumPhi / 2) {
                src_center.y = -src_center.y;
              }
              transmittances[j_src] =
                  GetTransmittance(src_center, target0.center);
            }
          }

          for (int j_tgt = 0; j_tgt < kNumPhi / 2; ++j_tgt) {
            Cell& target = cells_[i_tgt][j_tgt][k_tgt];
            for (int j_src = 0; j_src < kNumPhi; ++j_src) {
              int j = j_src < kNumPhi / 2 ? j_src : kNumPhi - 1 - j_src;
              Cell& src = cells_[i_src][j][k_src];
              if (i_src == i_tgt && j_src == j_tgt && k_src == k_tgt) {
                continue;
              }
              Position src_center = src.center;
              if (j_src >= kNumPhi / 2) {
                src_center.y = -src_center.y;
              }
              Position segment = src_center - target.center;
              Area segment_length_square = dot(segment, segment);
              Direction d = segment / sqrt(segment_length_square);
              InverseSolidAngle iso_phase = 1.0 / (4.0 * PI * sr);
              InverseSolidAngle mie_phase = MiePhaseFunction(dot(d, sun_dir));

              RadianceSpectrum src_radiance;
              if (double_scatter) {
                src_radiance =
                    (src.iso * iso_phase + src.aniso * mie_phase) *
                    (src.radial_width / src.volume);
              } else {
                src_radiance =
                    src.last * iso_phase * (src.radial_width / src.volume);
              }
              auto src_solid_angle = src.volume /
                  (src.radial_width * 4.0 * PI * segment_length_square) * sr;
              int transmittance_index = (j_src - j_tgt + kNumPhi) % kNumPhi;
              target.current += transmittances[transmittance_index] *
                  src_radiance * target0_scatter * src_solid_angle;
            }
          }
          progress_bar.Increment(1);
        }
      }
    }
  }, kNumTheta);
}

void Haber::InterpolateMultipleScatter(bool double_scatter) const {
  PowerSpectrum sum_over_all_cells = PowerSpectrum(0.0 * watt / nm);
  PowerSpectrum sum_over_active_cells = PowerSpectrum(0.0 * watt / nm);
  for (int i = 0; i < kNumTheta; ++i) {
    for (int j = 0; j < kNumPhi / 2; ++j) {
      std::vector<Cell>& cells = cells_[i][j];
      std::map<int, int> index_from_shell_index;
      for (unsigned int k = 0; k < cells.size(); ++k) {
        Cell& cell = cells[k];
        index_from_shell_index[cell.shell_index] = k;
        if (double_scatter) {
          if (cell.shell_index % 4 == 0) {
            sum_over_active_cells += cell.iso + cell.aniso;
          }
          sum_over_all_cells += cell.iso + cell.aniso;
        } else {
          if (cell.shell_index % 4 == 0) {
            sum_over_active_cells += cell.last;
          }
          sum_over_all_cells += cell.last;
        }
      }
      for (unsigned int k = 0; k < cells.size(); ++k) {
        Cell& cell = cells[k];
        if (cell.shell_index % 4 != 0) {
          int shell_index0 = 4 * (cell.shell_index / 4);
          int shell_index1 = shell_index0 + 4;
          int index0 = Lookup(index_from_shell_index, shell_index0);
          int index1 = Lookup(index_from_shell_index, shell_index1);
          if (index0 >= 0 && index1 >= 0) {
            double u = (cell.shell_index % 4) / 4.0;
            cell.current =
                cells[index0].current * (1.0 - u) + cells[index1].current * u;
          } else if (index0 >= 0) {
            cell.current = cells[index0].current;
          } else if (index1 >= 0) {
            cell.current = cells[index1].current;
          } else {
            cell.current = PowerSpectrum(0.0 * watt / nm);
          }
        }
      }
    }
  }
  Number scale = Integral(sum_over_all_cells) / Integral(sum_over_active_cells);
  for (int i = 0; i < kNumTheta; ++i) {
    for (int j = 0; j < kNumPhi / 2; ++j) {
      std::vector<Cell>& cells = cells_[i][j];
      for (unsigned int k = 0; k < cells.size(); ++k) {
        Cell& cell = cells[k];
        cell.current = cell.current * scale;
      }
    }
  }
}

void Haber::AccumulateMultipleScatter() const {
  for (int i = 0; i < kNumTheta; ++i) {
    for (int j = 0; j < kNumPhi / 2; ++j) {
      std::vector<Cell>& cells = cells_[i][j];
      for (unsigned int k = 0; k < cells.size(); ++k) {
        Cell& cell = cells[k];
        cell.iso = cell.iso + cell.current;
        cell.last = cell.current;
        cell.current = PowerSpectrum(0.0 * watt / nm);
      }
    }
  }
}

void Haber::ComputeSkyDome(Angle sun_zenith,
    ScatteringType scattering_type) const {
  Position viewer = Position(0.0 * m, 0.0 * m, EarthRadius);
  for (int i = 0; i < kNumTheta; ++i) {
    Angle view_zenith = (i + 0.5) * kDeltaPhi;
    for (int j = 0; j < kNumPhi / 2; ++j) {
      Angle view_azimuth = (j + 0.5) * kDeltaPhi;
      Angle view_sun = GetViewSunAngle(sun_zenith, view_zenith, view_azimuth);
      InverseSolidAngle iso_phase = 1.0 / (4.0 * PI * sr);
      InverseSolidAngle mie_phase = MiePhaseFunction(view_sun);
      RadianceSpectrum radiance =
          RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm);
      std::vector<Cell>& cells = cells_[i][j];
      for (unsigned int k = 0; k < cells.size(); ++k) {
        Cell& cell = cells[k];
        auto inverse_area = cell.radial_width / cell.volume;
        if (scattering_type == DOUBLE_SCATTERING_ONLY) {
          radiance += GetTransmittance(viewer, cell.center) *
              (cell.last * iso_phase) * inverse_area;
        } else {
          radiance += GetTransmittance(viewer, cell.center) *
              (cell.iso * iso_phase + cell.aniso * mie_phase) * inverse_area;
        }
      }
      sky_dome_.Set(i, j, radiance);
    }
  }
}

DimensionlessSpectrum Haber::GetTransmittance(const Position& p,
    const Position& q) const {
  // Compute where the line (p,q) is the closest to the Earth center.
  Number t = dot(p, p - q) / dot(p - q, p - q);
  if (t > 0.0 && t < 1.0) {
    // If it is between p and q, compute the transmittance in two steps, unless
    // the nearest distance to the center is less than the Earth radius.
    Position r = p + (q - p) * t;
    if (length(r) <= EarthRadius) {
      return DimensionlessSpectrum(0.0);
    } else {
      return GetTransmittanceSimple(p, r) * GetTransmittanceSimple(r, q);
    }
  } else {
    return GetTransmittanceSimple(p, q);
  }
}

DimensionlessSpectrum Haber::GetTransmittanceSimple(const Position& p,
    const Position& q) const {
  // We assume here that the segment p,q traverses each layer only once. This is
  // ensured by GetTransmittance above, which splits the segment in two if
  // necessary.
  int start_layer_index = GetLayerIndex(p);
  int end_layer_index = GetLayerIndex(q);
  Length pq = length(q - p);
  Length b;
  Area c;
  if (start_layer_index > end_layer_index) {
    std::swap(start_layer_index, end_layer_index);
    b = dot(q, p - q) / pq;
    c = dot(q, q);
  } else {
    b = dot(p, q - p) / pq;
    c = dot(p, p);
  }
  Length rayleigh_length = 0.0 * m;
  Length mie_length = 0.0 * m;
  Length distance_to_previous_layer = 0.0 * m;
  for (int i = start_layer_index; i < end_layer_index; ++i) {
    Length r_i = EarthRadius + GetLayerHeight(i + 1);
    Length distance_to_layer = -b + sqrt(b * b - c + r_i * r_i);
    Length layer_length = distance_to_layer - distance_to_previous_layer;
    rayleigh_length += rayleigh_density_[i] * layer_length;
    mie_length += mie_density_[i] * layer_length;
    distance_to_previous_layer = distance_to_layer;
  }
  rayleigh_length +=
      rayleigh_density_[end_layer_index] * (pq - distance_to_previous_layer);
  mie_length +=
      mie_density_[end_layer_index] * (pq - distance_to_previous_layer);

  DimensionlessSpectrum optical_depth = RayleighScattering() * rayleigh_length +
      MieExtinction() * mie_length;
  return exp(-optical_depth);
}

Haber::Position Haber::GetRayIntersectionWithLastLayer(const Position& p,
    Angle view_zenith) {
  constexpr Length max_layer_radius = EarthRadius + kMaxLayerHeight - 1.0 * m;
  Length b = p.x * sin(view_zenith) + p.z * cos(view_zenith);
  Area c = dot(p, p) - max_layer_radius * max_layer_radius;
  Length d = -b + sqrt(b * b - c);
  return Position(p.x + d * sin(view_zenith), p.y, p.z + d * cos(view_zenith));
}

Length Haber::GetLayerHeight(int layer_index) {
  constexpr Number M =
      (1.0 - exp(-kMaxLayerHeight / RayleighScaleHeight)) / kNumLayers;
  return -RayleighScaleHeight * log(1.0 - layer_index * M);
}

int Haber::GetLayerIndex(const Position &p) {
  constexpr Number M =
      (1.0 - exp(-kMaxLayerHeight/ RayleighScaleHeight)) / kNumLayers;
  Length height = length(p) - EarthRadius;
  return floor((1.0 - exp(-height / RayleighScaleHeight)) / M)();
}
