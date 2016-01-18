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
#include "atmosphere/measurement/measured_atmosphere.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace {

std::string ToString(double x) {
  std::stringstream out;
  out << x;
  if (out.str().find('.') == std::string::npos) {
    out << ".0";
  }
  return out.str();
}

std::string GetFileNameAndSampleLocation(const std::string& directory,
    const std::string& date, const std::string& hour,
    const std::string& minutes, int sample_index, int* x, int* y) {
  // Sample order used in Kider measurements, starts at North low on the
  // horizon, does a full turn, then a full turn in the opposite direction at a
  // higher elevation, and so on until the zenith is reached.
  static const int coords[2 * 81] = {
    // First turn.
    8, 4, 8, 5, 8, 6, 8, 7, 8, 8, 7, 8, 6, 8, 5, 8, 4, 8, 3, 8, 2, 8, 1, 8,
    0, 8, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1, 0, 0, 1, 0, 2, 0, 3, 0,
    4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 8, 1, 8, 2, 8, 3,
    // Second turn.
    7, 3, 7, 2, 7, 1, 6, 1, 5, 1, 4, 1, 3, 1, 2, 1, 1, 1, 1, 2, 1, 3, 1, 4,
    1, 5, 1, 6, 1, 7, 2, 7, 3, 7, 4, 7, 5, 7, 6, 7, 7, 7, 7, 6, 7, 5, 7, 4,
    // Third turn.
    6, 4, 6, 5, 6, 6, 5, 6, 4, 6, 3, 6, 2, 6, 2, 5, 2, 4, 2, 3, 2, 2,
    3, 2, 4, 2, 5, 2, 6, 2, 6, 3,
    // Fourth turn.
    5, 3, 4, 3, 3, 3, 3, 4, 3, 5, 4, 5, 5, 5, 5, 4,
    // Zenith.
    4, 4
  };
  *x = coords[2 * sample_index];
  *y = coords[2 * sample_index + 1];
  Angle view_zenith;
  Angle view_azimuth;
  HemisphericalFunction<RadianceSpectrum>::GetSampleDirection(
      *x, *y, &view_zenith, &view_azimuth);
  if (view_azimuth.to(rad) < 0.0) {
    view_azimuth =  view_azimuth + 2.0 * pi;
  }
  std::stringstream out;
  out << directory << "/" << date << "___" << hour << "." << minutes << ".00/"
      << date << "___" << hour << "." << minutes << ".00_-_" << sample_index
      << "____" << ToString(view_azimuth.to(deg)) << "___"
      << ToString(90.0 - view_zenith.to(deg)) << ".asd.rad.txt";
  return out.str();
}

RadianceSpectrum LoadKiderRadianceFile(const std::string& filename) {
  std::vector<Wavelength>  wavelengths;
  std::vector<SpectralRadiance> radiances;
  std::ifstream file(filename, std::ifstream::in);
  if (!file.good()) {
    std::cerr << "Cannot open file '" << filename << "'" << std::endl;
    ::exit(-1);
    return RadianceSpectrum(wavelengths, radiances);
  }
  std::string first_line_ignored;
  std::getline(file, first_line_ignored);
  while (file.good()) {
    int wavelength;
    double radiance;
    file >> wavelength >> radiance;
    wavelengths.push_back(wavelength * nm);
    radiances.push_back(radiance * watt_per_square_meter_per_sr_per_nm);
  }
  file.close();
  return RadianceSpectrum(wavelengths, radiances);
}

}  // anonymous namespace

MeasuredAtmosphere::MeasuredAtmosphere(const std::string& directory,
    const std::string& date, const std::string& hour,
    const std::string& minutes, Angle sun_zenith, Angle sun_azimuth,
    const std::string& cache_directory, bool compute_azimuth_from_data) {
  const std::string cache_file_name =
      cache_directory + "/" + date + hour + minutes + ".dat";
  if (cache_directory.length() > 0) {
    std::ifstream f;
    f.open(cache_file_name);
    if (f.good()) {
      f.close();
      measurements_.Load(cache_file_name);
      std::ifstream file(cache_file_name + ".params",
          std::ifstream::binary | std::ifstream::in);
      file.read(reinterpret_cast<char*>(&sun_zenith_), sizeof(sun_zenith_));
      file.read(reinterpret_cast<char*>(&sun_azimuth_), sizeof(sun_azimuth_));
      file.close();
      return;
    }
  }

  std::cout << "Loading " << hour << "h" << minutes << "..." << std::endl;
  for (int i = 0; i < 81; ++i) {
    int x, y;
    std::string filename =
        GetFileNameAndSampleLocation(directory, date, hour, minutes, i, &x, &y);
    // Azimuth in Kider data file names is measured counterclockwise from North,
    // while the usual convention is to measure it clockwise from North. To use
    // the normal azimuth convention in GetSkyRadiance, we need to invert the
    // azimuth sign, which translates to inverting the y coordinate of the
    // samples in the unit square / unit disk (where the y axis corresponds to
    // East).
    measurements_.Set(x, 8 - y, LoadKiderRadianceFile(filename));
  }
  // Measurement sample 0 (north) seems systematically wrong. Fix it by
  // interpolating the two nearest measurements at the same elevation.
  measurements_.Set(
      8, 4, (measurements_.Get(8, 3) + measurements_.Get(8, 5)) * 0.5);

  sun_zenith_ = sun_zenith;
  if (compute_azimuth_from_data) {
    // The theoretical sun zenith angle from the measurement location and date
    // corresponds to the data, but the azimuth angle does not, and the
    // difference between the theoretical and the actual values is not constant
    // (maybe a drift in the measuring instrument?). Thus, instead of using the
    // theoretical azimuth, we estimate it by finding the axis of symmetry of
    // the interpolated measurements.
    std::cout << "Computing sun direction..." << std::endl;
    sun_azimuth_ = EstimateSunAzimuth(sun_azimuth);
  } else {
    sun_azimuth_ = sun_azimuth;
  }

  // Cache the results on disk.
  if (cache_directory.length() > 0) {
    measurements_.Save(cache_file_name);
    std::ofstream file(cache_file_name + ".params", std::ofstream::binary);
    file.write(
        reinterpret_cast<const char*>(&sun_zenith_), sizeof(sun_zenith_));
    file.write(
        reinterpret_cast<const char*>(&sun_azimuth_), sizeof(sun_azimuth_));
    file.close();
  }
}

IrradianceSpectrum MeasuredAtmosphere::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  assert(altitude == 0.0 * m);
  assert(sun_zenith == sun_zenith_);
  return IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm);
}

RadianceSpectrum MeasuredAtmosphere::GetSkyRadiance(Length altitude,
    Angle sun_zenith, Angle view_zenith, Angle view_sun_azimuth) const {
  return GetSkyRadiance(altitude, sun_zenith, sun_azimuth_, view_zenith,
      sun_azimuth_ + view_sun_azimuth);
}

RadianceSpectrum MeasuredAtmosphere::GetSkyRadiance(Length altitude,
    Angle sun_zenith, Angle sun_azimuth, Angle view_zenith,
    Angle view_azimuth) const {
  assert(altitude == 0.0 * m);
  assert(sun_zenith == sun_zenith_);
  assert(sun_azimuth == sun_azimuth_);
  return measurements_(view_zenith, view_azimuth);
}

RadianceSpectrum MeasuredAtmosphere::GetSkyRadianceMeasurement(Length altitude,
    Angle sun_zenith, Angle sun_azimuth, Angle view_zenith,
    Angle view_azimuth) const {
  assert(altitude == 0.0 * m);
  assert(sun_zenith == sun_zenith_);
  assert(sun_azimuth == sun_azimuth_);

  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      Angle zenith_angle;
      Angle azimuth;
      HemisphericalFunction<RadianceSpectrum>::GetSampleDirection(
          i, j, &zenith_angle, &azimuth);
      Number cos_angle = cos(view_zenith) * cos(zenith_angle) +
          sin(view_zenith) * sin(zenith_angle) * cos(view_azimuth - azimuth);
      if (acos(cos_angle).to(deg) < 0.1) {
        return measurements_.Get(i, j);
      }
    }
  }
  return RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm);
}

Radiance MeasuredAtmosphere::GetAxiSymmetryError(Angle axis_azimuth) const {
  Radiance error = 0.0 * watt_per_square_meter_per_sr;
  for (int i = 0; i < 16; ++i) {
    Angle view_zenith = acos(Number((i + 0.5) / 16));
    for (int j = 0; j < 64; ++j) {
      Angle view_azimuth = j / 64.0 * 2.0 * pi;
      Radiance err = Integral(measurements_(view_zenith, view_azimuth) -
          measurements_(view_zenith, axis_azimuth * 2.0 - view_azimuth));
      error += err < 0.0 * watt_per_square_meter_per_sr ? -err : err;
    }
  }
  return error;
}

Angle MeasuredAtmosphere::EstimateSunAzimuth(Angle initial_azimuth) const {
  Angle min_sun_azimuth = initial_azimuth - 20.0 * deg;
  Angle max_sun_azimuth = initial_azimuth + 20.0 * deg;
  Radiance min_error = GetAxiSymmetryError(min_sun_azimuth);
  Angle best_sun_azimuth = min_sun_azimuth;
  Angle sun_azimuth = min_sun_azimuth;
  do {
    sun_azimuth = sun_azimuth + 0.2 * deg;
    Radiance error = GetAxiSymmetryError(sun_azimuth);
    if (error < min_error) {
      min_error = error;
      best_sun_azimuth = sun_azimuth;
    }
  } while (sun_azimuth < max_sun_azimuth);
  return best_sun_azimuth;
}
