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

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "atmosphere/measurement/measured_atmosphere.h"
#include "atmosphere/measurement/measured_atmospheres.h"
#include "atmosphere/model/bruneton/bruneton.h"
#include "atmosphere/model/haber/haber.h"
#include "atmosphere/model/hosek/hosek.h"
#include "atmosphere/model/libradtran/libradtran.h"
#include "atmosphere/model/nishita/nishita93.h"
#include "atmosphere/model/nishita/nishita96.h"
#include "atmosphere/model/oneal/oneal.h"
#include "atmosphere/model/polradtran/polradtran.h"
#include "atmosphere/model/preetham/preetham.h"
#include "atmosphere/comparisons.h"
#include "atmosphere/sun_direction.h"
#include "math/angle.h"
#include "physics/units.h"

namespace {

std::string ToString(int i) {
  std::stringstream out;
  if (i < 10) {
    out << "0";
  }
  out << i;
  return out.str();
}

std::string ToString(double x) {
  std::stringstream out;
  out << x;
  if (out.str().find('.') == std::string::npos) {
    out << ".0";
  }
  return out.str();
}

constexpr int kNumMeasurements = 17;
constexpr int kNumModels = 10;
constexpr int kNumViewSamples = 5;

const std::string kModels[kNumModels] = {
  "nishita93", "nishita96", "preetham", "oneal", "haber", "bruneton", "elek",
  "hosek", "libradtran", "measurements"
};

const std::string kCaptions[kNumModels] = {
  "Nishita93", "Nishita96", "Preetham", "O'Neal", "Haber", "Bruneton", "Elek",
  "Hosek", "libRadtran", "Measurements"
};

const std::string kLineStyle[kNumModels] = {
  "lt 1 lc rgbcolor \"#666666\" pt 9",   // nishita93
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 8",   // nishita96
  "lt 1 lc rgbcolor \"#999999\" pt 13",  // preetham
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 10",  // oneal
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 5",   // haber
  "lt 1 lc rgbcolor \"#666666\" pt 5",   // bruneton
  "lt 1 lc rgbcolor \"#666666\" pt 4",   // elek
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 12",  // hosek
  "lt 1 lc rgbcolor \"#0000ff\" pt 4",   // libradtran
  "lt 1 lc rgbcolor \"#ff0000\" lw 3"    // measurements
};

const std::string kSingleScatteringLineStyle[kNumModels] = {
  "",                                    // nishita93
  "lt 1 lc rgbcolor \"#666666\" pt 8",   // nishita96
  "",                                    // preetham
  "",                                    // oneal
  "lt 1 lc rgbcolor \"#666666\" pt 4",   // haber
  "lt 1 lc rgbcolor \"#666666\" pt 5",   // bruneton
  "",                                    // hosek
  "",                                    // libradtran
  ""                                     // measurements
};

const std::string kDoubleScatteringLineStyle[kNumModels] = {
  "",                                    // nishita93
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 8",   // nishita96
  "",                                    // preetham
  "",                                    // oneal
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 4",   // haber
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 5",   // bruneton
  "",                                    // hosek
  "",                                    // libradtran
  ""                                     // measurements
};

const std::string kMeasurements = "measurements";

const double kViewZenithSamples[kNumViewSamples] = {
  90.0 - 12.1151, 90.0 - 53.3665, 90.0 - 53.3665, 90.0 - 12.1151, 90.0 - 71.9187
};

const double kViewAzimuthSamples[kNumViewSamples] = {
  360.0 - 326.25, 360.0 - 67.5, 360.0 - 225.0, 360.0 - 225.0, 360.0 - 225.0
};

const int kIndices[] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
};

void SaveLibRadtranRmseTable(const std::string& libradtran_uvspec,
    const MeasuredAtmospheres& measurements,
    const std::vector<Angle>& sun_zenith, const std::vector<Angle>& sun_azimuth,
    Wavelength min_wavelength, Wavelength max_wavelength) {
  std::cout << "Computing libRadtran RMSE table..." << std::endl;
  auto index = [](double x) { return static_cast<int>(round(x * 1000)); };
  auto key = [index](double alpha, double beta, double g) {
    return std::make_pair(std::make_pair(index(alpha), index(beta)), index(g));
  };
  std::map<std::pair<std::pair<int, int>, int>, double> already_computed;
  std::ifstream input(
      Comparisons::GetOutputDirectory() + "libradtran_rmse.txt");
  if (input) {
    while (!input.eof()) {
      double alpha, beta, g, rmse;
      input >> alpha >> beta >> g >> rmse;
      already_computed.insert(std::make_pair(key(alpha, beta, g), rmse));
    }
    input.close();
  }

  SpectralRadiance min_rmse = 0.0 * watt_per_square_meter_per_sr_per_nm;
  double min_alpha, min_beta, min_g;
  std::ofstream output(
      Comparisons::GetOutputDirectory() + "libradtran_rmse.txt",
      std::ofstream::app);
  for (double alpha = 0.0; alpha <= 2.0; alpha += 0.2) {
    for (double beta = 0.02; beta <= 0.2; beta += 0.02) {
      for (double g = 0.5; g <= 0.9; g += 0.1) {
        SpectralRadiance rmse;
        auto it = already_computed.find(key(alpha, beta, g));
        if (it != already_computed.end()) {
          rmse = (it->second * 1e-3) * watt_per_square_meter_per_sr_per_nm;
        } else {
          std::cout << alpha << " " << beta << " " << g << std::endl;
          LibRadtran lib_radtran(
              libradtran_uvspec, alpha, beta, g, true /* ground_albedo */,
              LibRadtran::HEMISPHERICAL_FUNCTION_CACHE);
          Comparisons comparisons(
              "", lib_radtran, measurements, min_wavelength, max_wavelength);
          rmse = comparisons.ComputeTotalRmse(sun_zenith, sun_azimuth);
          output << alpha << " " << beta << " " << g << " "
                 << rmse.to(1e-3 * watt_per_square_meter_per_sr_per_nm)
                 << std::endl;
        }
        if (min_rmse == 0.0 * watt_per_square_meter_per_sr_per_nm ||
            rmse < min_rmse) {
          min_rmse = rmse;
          min_alpha = alpha;
          min_beta = beta;
          min_g = g;
        }
      }
    }
  }
  output.close();
  std::cout << "Smallest RMSE: "
            << min_rmse.to(1e-3 * watt_per_square_meter_per_sr_per_nm)
            << " at alpha=" << min_alpha << " beta=" << min_beta
            << " g=" << min_g << std::endl;
}

void SaveZenithLuminanceRmseTable(const MeasuredAtmospheres& measurements,
    const std::vector<Angle>& sun_zenith, Wavelength min_wavelength,
    Wavelength max_wavelength) {
  std::ifstream input(
      Comparisons::GetOutputDirectory() + "zenith_luminance_rmse.txt");
  if (input) {
    input.close();
  } else {
    std::cout << "Computing zenith luminance RMSE table..." << std::endl;
    std::ofstream output(
        Comparisons::GetOutputDirectory() + "zenith_luminance_rmse.txt");
    Comparisons comparisons("", measurements, measurements, min_wavelength,
        max_wavelength);
    for (double turbidity = 2.0; turbidity <= 4.0; turbidity += 0.01) {
      auto error_square_sum =
          0.0 * cd_per_square_meter * cd_per_square_meter;
      for (unsigned int i = 0; i < sun_zenith.size(); ++i) {
        Luminance measured = comparisons.ComputeZenithLuminance(sun_zenith[i]);
        // See Eq. 1 in "Zenith luminance and sky luminance distributions for
        // daylighting calculations", Karayel et al, Energy and Buildings, 1984.
        Angle gamma = pi / 2.0 - sun_zenith[i];
        Luminance model =
            ((1.376 * turbidity - 1.81) * tan(gamma) + 0.38) *
                kcd_per_square_meter;
        error_square_sum += (measured - model) * (measured - model);
      }
      Luminance rmse = sqrt(error_square_sum / sun_zenith.size());
      output << turbidity << " "
             << rmse.to(cd_per_square_meter)
             << std::endl;
    }
    output.close();
  }
}

void SavePreethamRmseTable(const MeasuredAtmospheres& measurements,
    const std::vector<Angle>& sun_zenith, const std::vector<Angle>& sun_azimuth,
    Wavelength min_wavelength, Wavelength max_wavelength) {
  std::ifstream input(Comparisons::GetOutputDirectory() + "preetham_rmse.txt");
  if (input) {
    input.close();
  } else {
    std::cout << "Computing Preetham RMSE table..." << std::endl;
    std::ofstream output(
        Comparisons::GetOutputDirectory() + "preetham_rmse.txt");
    for (double turbidity = 2.0; turbidity <= 4.0; turbidity += 0.01) {
      SpectralRadiance rmse = Comparisons("", Preetham(turbidity), measurements,
          min_wavelength, max_wavelength).ComputeTotalRmse(
              sun_zenith, sun_azimuth);
      output << turbidity << " "
             << rmse.to(1e-3 * watt_per_square_meter_per_sr_per_nm)
             << std::endl;
    }
    output.close();
  }
}

void SaveHosekRmseTable(const MeasuredAtmospheres& measurements,
    const std::vector<Angle>& sun_zenith, const std::vector<Angle>& sun_azimuth,
    Wavelength min_wavelength, Wavelength max_wavelength) {
  std::ifstream input(Comparisons::GetOutputDirectory() + "hosek_rmse.txt");
  if (input) {
    input.close();
  } else {
    std::cout << "Computing Hosek RMSE table..." << std::endl;
    std::ofstream output(Comparisons::GetOutputDirectory() + "hosek_rmse.txt");
    for (double turbidity = 2.0; turbidity <= 4.0; turbidity += 0.01) {
      SpectralRadiance rmse = Comparisons("", Hosek(turbidity), measurements,
          min_wavelength, max_wavelength).ComputeTotalRmse(
              sun_zenith, sun_azimuth);
      output << turbidity << " "
             << rmse.to(1e-3 * watt_per_square_meter_per_sr_per_nm)
             << std::endl;
    }
    output.close();
  }
}

void SaveRadiances(const Comparisons& comparisons,
    const std::vector<std::string>& name, const std::vector<Angle>& sun_zenith,
    const std::vector<Angle>& sun_azimuth) {
  std::cout << "Computing radiance..." << std::endl;
  for (int i : kIndices) {
    if (i == 9) {
      for (int x = 0; x < 9; ++x) {
        for (int y = 0; y < 9; ++y) {
          Angle view_zenith;
          Angle view_azimuth;
          HemisphericalFunction<Number>::GetSampleDirection(
              x, y, &view_zenith, &view_azimuth);
          comparisons.PlotRadiance(name[i], sun_zenith[i], sun_azimuth[i],
              view_zenith, view_azimuth);
        }
      }
    } else {
      for (int j = 0; j < kNumViewSamples; ++j) {
        Angle view_zenith = kViewZenithSamples[j] * deg;
        Angle view_azimuth = kViewAzimuthSamples[j] * deg;
        comparisons.PlotRadiance(
            name[i], sun_zenith[i], sun_azimuth[i], view_zenith, view_azimuth);
      }
    }
  }
}

void SaveComparisons(const Comparisons& comparisons,
    const std::vector<std::string>& name, const std::vector<Angle>& sun_zenith,
    const std::vector<Angle>& sun_azimuth, Location measurement_location,
    Time measurement_date, bool save_radiances, bool is_measurement) {
  if (!is_measurement) {
    std::cout << "Computing sun illuminance attenuation..." << std::endl;
    comparisons.PlotSunIlluminanceAttenuation();

    std::cout << "Rendering sky images..." << std::endl;
    measurement_date.hours = 6 + 4;
    measurement_date.minutes = 0;
    SunCoordinates sun_direction;
    GetSunDirection(measurement_date, measurement_location, &sun_direction);
    comparisons.RenderSkyImage(
        "sunrise", sun_direction.zenith * deg, -30.0 * deg);

    measurement_date.hours = 9 + 4;
    measurement_date.minutes = 30;
    GetSunDirection(measurement_date, measurement_location, &sun_direction);
    comparisons.RenderSkyImage(
        "morning", sun_direction.zenith * deg, 30.0 * deg);
  }

  std::vector<Luminance> zenith_luminance;
  std::vector<Irradiance> sun_irradiance;
  std::vector<Irradiance> sky_irradiance;
  int count = 0;
  auto error_square_sum = 0.0 * watt_per_square_meter_per_sr_per_nm *
      watt_per_square_meter_per_sr_per_nm;
  auto error_square_sum_with_approximate_spectrum = error_square_sum;
  for (int i : kIndices) {
    if (save_radiances) {
      std::cout << "Computing radiance..." << std::endl;
      if (i == 9) {
        for (int x = 0; x < 9; ++x) {
          for (int y = 0; y < 9; ++y) {
            Angle view_zenith;
            Angle view_azimuth;
            HemisphericalFunction<Number>::GetSampleDirection(
                x, y, &view_zenith, &view_azimuth);
            comparisons.PlotRadiance(name[i], sun_zenith[i], sun_azimuth[i],
                view_zenith, view_azimuth);
          }
        }
      } else {
        for (int j = 0; j < kNumViewSamples; ++j) {
          Angle view_zenith = kViewZenithSamples[j] * deg;
          Angle view_azimuth = kViewAzimuthSamples[j] * deg;
          comparisons.PlotRadiance(name[i], sun_zenith[i], sun_azimuth[i],
              view_zenith, view_azimuth);
        }
      }
    }
    std::cout << "Computing luminance and image..." << std::endl;
    comparisons.RenderLuminanceAndImage(name[i], sun_zenith[i], sun_azimuth[i]);
    std::cout << "Computing luminance profile..." << std::endl;
    comparisons.PlotLuminanceProfile(name[i], sun_zenith[i], sun_azimuth[i]);
    if (!is_measurement) {
      std::cout << "Computing relative error..." << std::endl;
      SpectralRadiance rmse_with_approximate_spectrum;
      SpectralRadiance rmse = comparisons.PlotRelativeError(name[i],
          sun_zenith[i], sun_azimuth[i], &rmse_with_approximate_spectrum);
      count += 1;
      error_square_sum += rmse * rmse;
      error_square_sum_with_approximate_spectrum +=
          rmse_with_approximate_spectrum * rmse_with_approximate_spectrum;
    }
    std::cout << "Computing zenith luminance..." << std::endl;
    zenith_luminance.push_back(
        comparisons.ComputeZenithLuminance(sun_zenith[i]));
    std::cout << "Computing irradiance..." << std::endl;
    Irradiance sun, sky;
    comparisons.ComputeIrradiance(sun_zenith[i], &sun, &sky);
    sun_irradiance.push_back(sun);
    sky_irradiance.push_back(sky);
  }
  std::cout << "Computing sunrise luminance and image..." << std::endl;
  for (int i = 0; i < 3; ++i) {
    int minutes = (6 + i) * 60;
    measurement_date.hours = minutes / 60 + 4;
    measurement_date.minutes = minutes % 60;
    SunCoordinates sun_direction;
    GetSunDirection(measurement_date, measurement_location, &sun_direction);
    std::string name = ToString(minutes / 60) + "h" + ToString(minutes % 60);
    if (is_measurement) {
      std::ofstream sza(
          Comparisons::GetOutputDirectory() + "sza_" + name + ".txt");
      sza << round(sun_direction.zenith) << std::endl;
      sza.close();
    } else {
      comparisons.RenderLuminanceAndImage(name, sun_direction.zenith * deg,
          sun_direction.azimuth * deg);
    }
  }

  comparisons.PlotDayZenithLuminance(sun_zenith, zenith_luminance);
  comparisons.PlotDayIrradiance(sun_zenith, sun_irradiance, sky_irradiance);
  if (!is_measurement) {
    SpectralRadiance rmse = sqrt(error_square_sum / count);
    double rounded_rmse =
        round(rmse.to(1e-4 * watt_per_square_meter_per_sr_per_nm)) / 10.0;
    std::ofstream file_percent(Comparisons::GetOutputDirectory() +
        "error_total_" + comparisons.name() + ".txt");
    file_percent << rounded_rmse << std::endl;
    file_percent.close();
  }
  if (!is_measurement) {
    SpectralRadiance rmse_with_approximate_spectrum =
        sqrt(error_square_sum_with_approximate_spectrum / count);
    double rounded_rmse = round(rmse_with_approximate_spectrum.to(
        1e-4 * watt_per_square_meter_per_sr_per_nm)) / 10.0;
    std::ofstream file_percent(Comparisons::GetOutputDirectory() +
        "error_total_with_approximate_spectrum_" + comparisons.name() + ".txt");
    file_percent << rounded_rmse << std::endl;
    file_percent.close();
  }
}

void SaveSolarSpectrum() {
  IrradianceSpectrum solar = SolarSpectrum();
  std::ofstream file(Comparisons::GetOutputDirectory() + "solar_spectrum.txt");
  for (unsigned int i = 0; i < solar.size(); ++i) {
    file << solar.GetWavelength(i).to(nm) << " "
        << solar[i].to(watt_per_square_meter_per_nm) << std::endl;
  }
  file.close();
}

double CornetteShanks(double g, double theta_radians) {
  double k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
  double mu = cos(theta_radians);
  return k * (1.0 + mu * mu) / pow(1.0 + g * g - 2.0 * g * mu, 1.5);
}

// The asymmetry factor, integral of the phase function times cos(theta) over
// all solid angles, is NOT equal to g (unlike for the Henyey-Greenstein phase
// function).
double CornetteShanksAsymmetryFactor(double g) {
  const int N = 500;
  double result = 0.0;
  for (unsigned int i = 1; i < N; ++i) {
    double theta0 = static_cast<double>(i - 1) / N * PI;
    double theta1 = static_cast<double>(i) / N * PI;
    double value0 = CornetteShanks(g, theta0) * sin(theta0) * cos(theta0);
    double value1 = CornetteShanks(g, theta1) * sin(theta1) * cos(theta1);
    result += (value1 + value0) * (theta1 - theta0) * PI;
  }
  return result;
}

// Returns the parameter g such that CornetteShanks(g) has the given asymmetry
// factor.
double CornetteShanksG(double asymmetry_factor) {
  assert(asymmetry_factor > 0.0);
  // Naive search: go through all g between 0 and 1 and returns the best match.
  const int N = 100;
  double best_g = 0.0;
  double asymmetry_factor_error_for_best_g = asymmetry_factor;
  for (int i = 1; i < N; ++i) {
    double g = static_cast<double>(i) / N;
    double asymmetry_factor_error_for_g =
        std::abs(CornetteShanksAsymmetryFactor(g) - asymmetry_factor);
    if (asymmetry_factor_error_for_g < asymmetry_factor_error_for_best_g) {
      best_g = g;
      asymmetry_factor_error_for_best_g = asymmetry_factor_error_for_g;
    }
  }
  return best_g;
}

void SaveMiePhaseFunction() {
  std::ifstream input("input/130527_130527_Egbert.pfn");
  std::string line;
  for (int i = 0; i < 3; ++i) {
    std::getline(input, line);
  }
  std::vector<double> values[2];
  for (int i = 0; i < 2; ++i) {
    std::getline(input, line);
    std::stringstream line_stream(line);
    std::string value;
    while (getline(line_stream, value, ',')) {
      std::stringstream value_stream(value);
      double double_value;
      value_stream >> double_value;
      values[i].push_back(i == 0 ? double_value : double_value / (4.0 * PI));
    }
  }
  input.close();
  double asymmetry_factor = 0.0;
  for (unsigned int i = 4; i < values[0].size(); ++i) {
    double theta0 = values[0][i - 1] / 180.0 * PI;
    double theta1 = values[0][i] / 180.0 * PI;
    double value0 = values[1][i - 1] * sin(theta0) * cos(theta0);
    double value1 = values[1][i] * sin(theta1) * cos(theta1);
    asymmetry_factor += (value1 + value0) * (theta1 - theta0) * PI;
    if (values[0][i] == 180.0) {
      break;
    }
  }
  double g = CornetteShanksG(asymmetry_factor);
  std::ofstream file(Comparisons::GetOutputDirectory() + "phasefunction.txt");
  for (unsigned int i = 3; i < values[0].size(); ++i) {
    double theta = values[0][i] / 180.0 * PI;
    double mie = CornetteShanks(g, theta);
    file << values[0][i] << " " << values[1][i] << " " << mie << std::endl;
    if (values[0][i] == 180.0) {
      break;
    }
  }
  file.close();
}

void SaveSubPlot(std::ofstream* file, Angle sun_zenith, Angle sun_azimuth,
    Angle view_zenith, Angle view_azimuth, bool profile) {
  (*file) << "set size square 0.47\n"
      << (profile ? "set origin 0.15,0.5\n" : "set origin 0.59,0.50\n")
      << "set object circle at 0.5,0.5 radius screen 0.103 behind "
      "fc rgb \"white\" fs solid\n"
      "set polar\n"
      "unset raxis\n"
      "set angle degrees\n"
      "set grid xtics noytics polar 45 linestyle 1\n"
      "set xtics axis (\"\" 10, \"\" 20, \"\" 30, \"\" 40, \"\" 50, "
      "\"\" 60, \"\" 70, \"\" 80, \"\" 90) scale 0\n"
      "set xrange [-100:100]\n"
      "set yrange [-100:100]\n"
      "set label \"E\" at 80,0 center\n"
      "set label \"S\" at 0,75 center\n"
      "set label \"W\" at -75,0 center\n"
      "set label \"N\" at 0,-75 center\n"
      "set pointsize 2\n";
  if (profile) {
    (*file) << "plot '-' with points lc rgbcolor \"#f0f000\" pt 7, "
        "'-' with lines lt 1 lc rgbcolor \"#000000\" lw 3, "
        "90 with lines lt 1 lc rgbcolor \"#aaaaaa\"\n"
        << (sun_azimuth.to(deg) - 90) << " " << sun_zenith.to(deg) << std::endl
        << "e\n"
        << (view_azimuth.to(deg) - 90) << " -90\n"
        << (view_azimuth.to(deg) - 90) << " +90\n"
        << "e\n";
  } else {
    (*file) << "plot '-' with points lc rgbcolor \"#f0f000\" pt 7, "
        "'-' with points lc rgbcolor \"#000000\" lw 3 pt 4, "
        "90 with lines lt 1 lc rgbcolor \"#aaaaaa\"\n"
        << (sun_azimuth.to(deg) - 90) << " " << sun_zenith.to(deg) << std::endl
        << "e\n" << (view_azimuth.to(deg) - 90) << " " << view_zenith.to(deg)
        << std::endl << "e\n";
  }
}

void SavePlot(const std::vector<MeasuredAtmosphere*>& measurements,
    const std::vector<std::string>& names,
    const std::vector<Angle>& sun_zenith,
    const std::vector<Angle>& sun_azimuth) {
  std::ofstream file(Comparisons::GetOutputDirectory() + "main.plot");
  file << "load \"" << Comparisons::GetOutputDirectory()
       << "error_caption.plot\"\n"
       << "load \"" << Comparisons::GetOutputDirectory()
       << "scale_caption.plot\"\n"
       << "reset\n"
       << "set terminal postscript eps size 6.8cm,4cm \"NimbusSanL-Regu\"\n"
       << "set tics in nomirror scale 0.2\n"
       << "set style line 1 lc rgbcolor \"#eeeeee\"\n"
       << "set grid noxtics ytics linestyle 1\n\n";

  file << "set output \"" << Comparisons::GetOutputDirectory()
       << "solar_spectrum.eps\"\n"
       << "set xlabel \"Wavelength (nm)\"\n"
       << "set ylabel \"Spectral Irradiance (W/m2/nm)\"\n"
       << "plot [360:830] \"input/astm-g173.txt\" with lines t \"ASTM-G173\", "
       << "\"" << Comparisons::GetOutputDirectory() << "solar_spectrum.txt\" "
       << "t \"Our Solar Spectrum\" with lines\n\n";

  file << "set output \"" << Comparisons::GetOutputDirectory()
       << "phase_function.eps\"\n"
       << "set xlabel \"Scattering angle (degrees)\"\n"
       << "set ylabel \"Value (unitless)\"\n"
       << "set logscale y\n"
       << "plot \"" << Comparisons::GetOutputDirectory()
       << "phasefunction.txt\" using 1:2 with lines "
       << "t \"Measured\", \"\" using 1:3 with lines t \"Cornette-Shanks\", "
       << "\"\" using 1:($3/$2) with lines t \"Ratio\"\n"
       << "unset logscale y\n\n";

  file << "set output \"" << Comparisons::GetOutputDirectory()
       << "sun_illuminance_attenuation.eps\"\n"
       << "set xlabel \"Sun zenith angle (degrees)\"\n"
       << "set ylabel \"Sun illuminance attenuation\"\n"
       << "plot [0:90][0:1] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (kModels[i] == kMeasurements) {
      continue;
    } else if (i > 0) {
      file << ", ";
    }
    file << "\"" << Comparisons::GetOutputDirectory()
         << "sun_illuminance_attenuation_" << kModels[i] << ".txt\" "
         << "t \"" << kCaptions[i] << "\" with lines";
  }
  file << std::endl << std::endl;

  file << "set terminal postscript eps size 8.5cm,5.0cm \"NimbusSanL-Regu\"\n";
  for (int i : kIndices) {
    file << "set output \"" << Comparisons::GetOutputDirectory()
         << "luminance_profile_" << names[i] << ".eps\"\n"
         << "reset\n"
         << "set title \"" << names[i] << " / "
         << round(sun_zenith[i].to(deg)) << "^o\" offset 17,-3\n"
         << "set object rectangle from screen 0,0 to screen 1,1 behind "
             "fc rgb \"white\" fs solid noborder\n"
             "set multiplot\n"
             "set tics in nomirror scale 0.5\n"
             "set style line 1 lc rgbcolor \"#eeeeee\"\n"
             "set grid noxtics ytics linestyle 1\n"
             "set key outside center bottom horizontal Left reverse samplen 3 "
             "width -3 maxrows 2\n";
    file << "set pointsize 0.75\n";
    file << "plot [-90:90][0:25000] ";
    for (int j = 0; j < kNumModels; ++j) {
      if (j > 0) {
        file << ", ";
      }
      file << "\"" << Comparisons::GetOutputDirectory()
           << "luminance_profile_" << names[i] << "_" << kModels[j]
           << ".txt\" t \"" << kCaptions[j] << "\"";
      if (kModels[j] == kMeasurements) {
        file << " with points lc rgbcolor \"#ff0000\" lw 3 pt 5";
      } else {
        file << " with linespoints " << kLineStyle[j];
      }
    }
    file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
        "unset border\nunset key\nunset object\n unset title\n";
    SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i], 0 * deg, 135 * deg, true);
    file << "unset multiplot\n\n";
  }

  int i = 9;
  for (int x = 0; x < 9; ++x) {
    for (int y = 0; y < 9; ++y) {
      Angle view_zenith;
      Angle view_azimuth;
      HemisphericalFunction<Number>::GetSampleDirection(
          x, y, &view_zenith, &view_azimuth);
      file << "set output \"" << Comparisons::GetOutputDirectory()
           << "radiance_" << names[i] << "_" << x << "_" << y << ".eps\"\n";
      file << "reset\n"
          "set multiplot\n"
          "set object rectangle from screen 0,0 to screen 1,1 behind "
          "fc rgb \"white\" fs solid noborder\n"
          "set tics in nomirror scale 0.5\n"
          "set style line 1 lc rgbcolor \"#eeeeee\"\n"
          "set pointsize 1.25\n"
          "set grid noxtics ytics linestyle 1\n"
          "set key outside center bottom horizontal Left reverse samplen 3 "
          "width -3 maxrows 2\n";
      file << "plot [360:720][0:0.28] ";
      for (int k = 0; k < kNumModels; ++k) {
        if (k > 0) {
          file << ", ";
        }
        if (kModels[k] == kMeasurements) {
          file << "\"" << measurements[i]->GetSourceFileName(x, y) << "\""
               << " every ::1 t \"" << kCaptions[k]
               << "\" with lines " << kLineStyle[k];
        } else {
          file << "\"" << Comparisons::GetOutputDirectory()
               << "radiance_" << names[i] << "_" << view_zenith.to(deg) << "_"
               << view_azimuth.to(deg) << "_" << kModels[k]
               << ".txt\" t \"" << kCaptions[k] << "\""
               << " with linespoints " << kLineStyle[k];
        }
      }
      file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
          "unset border\nunset key\nunset object\n";
      SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i],
          view_zenith, view_azimuth, false);
      file << "unset multiplot" << std::endl << std::endl;
    }
  }

  for (int i : kIndices) {
    for (int j = 0; j < kNumViewSamples; ++j) {
      file << "set output \"" << Comparisons::GetOutputDirectory()
           << "albedo_impact_" << names[i] << "_" << kViewZenithSamples[j]
           << "_" << kViewAzimuthSamples[j] << ".eps\"" << std::endl;
      file << "reset\n"
          "set multiplot\n"
          "set object rectangle from screen 0,0 to screen 1,1 behind "
          "fc rgb \"white\" fs solid noborder\n"
          "set tics in nomirror scale 0.5\n"
          "set style line 1 lc rgbcolor \"#eeeeee\"\n"
          "set pointsize 1.25\n"
          "set grid noxtics ytics linestyle 1\n"
          "set key outside center bottom horizontal Left reverse maxrows 2\n";
      file << "plot [360:720][0:0.28] ";
      file << "\"" << Comparisons::GetOutputDirectory()
           << "radiance_" << names[i] << "_" << kViewZenithSamples[j] << "_"
           << kViewAzimuthSamples[j] << "_libradtran.txt\" "
           << "t \"libRadtran with ground albedo\" with lines, ";
      file << "\"" << Comparisons::GetOutputDirectory()
           << "radiance_" << names[i] << "_" << kViewZenithSamples[j] << "_"
           << kViewAzimuthSamples[j] << "_libradtran_no_ground_albedo.txt\" "
           << "t \"libRadtran without ground albedo\" with lines, ";
      file << "\"" << measurements[i]->GetSourceFileName(
               kViewZenithSamples[j] * deg, kViewAzimuthSamples[j] * deg)
           << "\" every ::1 t \"" << kCaptions[kNumModels - 1]
           << "\" with lines " << kLineStyle[kNumModels - 1];
      file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
          "unset border\nunset key\nunset object\n";
      SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i],
          kViewZenithSamples[j] * deg, kViewAzimuthSamples[j] * deg, false);
      file << "unset multiplot" << std::endl << std::endl;
    }
  }

  for (int i : kIndices) {
    for (int j = 0; j < kNumViewSamples; ++j) {
      file << "set output \"" << Comparisons::GetOutputDirectory()
           << "ms_radiance_" << names[i] << "_" << kViewZenithSamples[j] << "_"
           << kViewAzimuthSamples[j] << ".eps\"" << std::endl;
      file << "reset\n"
          "set multiplot\n"
          "set object rectangle from screen 0,0 to screen 1,1 behind "
          "fc rgb \"white\" fs solid noborder\n"
          "set tics in nomirror scale 0.5\n"
          "set style line 1 lc rgbcolor \"#eeeeee\"\n"
          "set pointsize 1.25\n"
          "set grid noxtics ytics linestyle 1\n"
          "set key outside center bottom horizontal Left reverse maxrows 2\n";
      file << "plot [360:720][0:0.28] ";
      bool is_first_iteration = true;
      for (int k = 0; k < kNumModels; ++k) {
        if (kModels[k] != "bruneton" && kModels[k] != "haber" &&
            kModels[k] != "nishita96") {
          continue;
        }
        if (!is_first_iteration) {
          file << ", ";
        }
        file << "\"" << Comparisons::GetOutputDirectory()
             << "radiance_" << names[i] << "_" << kViewZenithSamples[j] << "_"
             << kViewAzimuthSamples[j] << "_" << kModels[k] << "_ss.txt\" t \""
             << kCaptions[k] << " SS\"" << " with linespoints "
             << kSingleScatteringLineStyle[k];
        is_first_iteration = false;
      }
      for (int k = 0; k < kNumModels; ++k) {
        if (kModels[k] != "bruneton" && kModels[k] != "haber" &&
            kModels[k] != "nishita96") {
          continue;
        }
        file << ", \"" << Comparisons::GetOutputDirectory()
             << "radiance_" << names[i] << "_" << kViewZenithSamples[j] << "_"
             << kViewAzimuthSamples[j] << "_" << kModels[k] << "_ds.txt\" t \""
             << kCaptions[k] << " DS\"" << " with linespoints "
             << kDoubleScatteringLineStyle[k];
      }
      file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
          "unset border\nunset key\nunset object\n";
      SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i],
          kViewZenithSamples[j] * deg, kViewAzimuthSamples[j] * deg, false);
      file << "unset multiplot" << std::endl << std::endl;
    }
  }

  file << "reset\n";

  file << "set output \"" << Comparisons::GetOutputDirectory()
       << "day_zenith_luminance.eps\"\n"
          "set object rectangle from screen 0,0 to screen 1,1 behind "
          "fc rgb \"white\" fs solid noborder\n"
          "set tics in nomirror scale 0.5\n"
          "set style line 1 lc rgbcolor \"#eeeeee\"\n"
          "set grid noxtics ytics linestyle 1\n"
          "set key outside center bottom horizontal Left reverse samplen 3 "
          "width -3 maxrows 2\n"
          "set ylabel \"Zenith Luminance (cd/m^2)\"\n"
          "set xlabel \"Time\"\n"
          "set xtics (";
  for (int i = 0; i < kNumMeasurements; i += 2) {
    if (i > 0) {
      file << ", ";
    }
    file << "\"" << names[i] << "\"" << " " << i;
  }
  file << ")" << std::endl;
  file << "plot [][0:] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (i > 0) {
      file << ", ";
    }
    file << "\"" << Comparisons::GetOutputDirectory()
         << "day_zenith_luminance_" << kModels[i] << ".txt\" "
         << "t \"" << kCaptions[i] << "\" with "
         << (kModels[i] == kMeasurements ? "lines " : "linespoints ")
         << kLineStyle[i];
  }
  file << std::endl << std::endl;

  file << "set output \"" << Comparisons::GetOutputDirectory()
      << "day_irradiance_sky.eps\"\n";
  file << "unset xlabel\n unset ylabel\n";
  file << "plot [][0:] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (i > 0) {
      file << ", ";
    }
    file << "\"" << Comparisons::GetOutputDirectory()
         << "day_irradiance_" << kModels[i] << ".txt\" "
         << "using 1:3 t \"" << kCaptions[i] << "\" with "
         << (kModels[i] == kMeasurements ? "lines " : "linespoints ")
         << kLineStyle[i];
  }
  file << std::endl << std::endl;

  file << "set output \"" << Comparisons::GetOutputDirectory()
       << "day_irradiance_total.eps\"\n";
  file << "unset xlabel\n unset ylabel\n";
  file << "plot [][0:] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (kModels[i] == kMeasurements || kModels[i] == "libradtran") {
      continue;
    } else if (i > 0) {
      file << ", ";
    }
    file << "\"" << Comparisons::GetOutputDirectory()
         << "day_irradiance_" << kModels[i] << ".txt\" "
         << "using 1:2 t \"" << kCaptions[i]
         << "\" with linespoints " << kLineStyle[i];
  }
  file << std::endl << std::endl;

  file << "set terminal postscript eps size 10cm,5.0cm\n";
  file << "set output \"" << Comparisons::GetOutputDirectory()
       << "turbidity_rmse.eps\"\n";
  file << "reset\n"
      "set xrange [2:3]\n"
      "set xtics nomirror\n"
      "set yrange [0:]\n"
      "set ytics nomirror\n"
      "set y2tics\n"
      "set xlabel \"Turbidity\"\n"
      "set ylabel \"Zenith luminance (cd/m^2)\"\n"
      "set y2label \"Spectral radiance (mW/m^2/sr/nm)\"\n"
      "set key center top Left reverse\n"
      "plot \"";
  file << Comparisons::GetOutputDirectory()
       << "zenith_luminance_rmse.txt\" with lines axes x1y1 "
       << "t \"Karayek et al. RMSE (luminance)\", \""
       << Comparisons::GetOutputDirectory()
       << "preetham_rmse.txt\" with lines axes x1y2 "
       << "t \"Preetham RMSE (spectral radiance)\", \""
       << Comparisons::GetOutputDirectory()
       << "hosek_rmse.txt\" using 1:($2*2) with lines axes x1y2 "
       << "t \"Hosek RMSE * 2 (spectral radiance)\"";
  file << std::endl << std::endl;

  file << Comparisons::SaveErrorCaption() << std::endl;
  file << Comparisons::SaveScaleCaption() << std::endl;
  file << Comparisons::SaveGroundAlbedo() << std::endl;
  file.close();
}

}  // anonymous namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <libRatran uvspec path> <libRadtran data path>" << std::endl;
    return -1;
  }
  const std::string libradtran_uvspec(argv[1]);
  const std::string libradtran_data(argv[2]);

  // The measurements were made at Cornell University, Frank H.T. Rhodes Hall,
  // in 2013/05/27.
  std::cout << "Measurements..." << std::endl;
  Location measurement_location;
  measurement_location.latitude = 42.44337;
  measurement_location.longitude = -76.48163;
  Time measurement_time;
  measurement_time.year = 2013;
  measurement_time.month = 5;
  measurement_time.day = 27;

  std::vector<std::string> name;
  std::vector<Angle> sun_zenith;
  std::vector<Angle> sun_azimuth;
  std::vector<MeasuredAtmosphere*> measurements;
  MeasuredAtmospheres measured;
  for (int i = 0; i < kNumMeasurements; ++i) {
    int minutes = 9 * 60 + 30 + i * 15;
    measurement_time.hours = minutes / 60 + 4;
    measurement_time.minutes = minutes % 60;
    measurement_time.seconds = 0;
    SunCoordinates sun_direction;
    GetSunDirection(measurement_time, measurement_location, &sun_direction);
    MeasuredAtmosphere* atmosphere = new MeasuredAtmosphere(
        "input", "2013-05-27", ToString(minutes / 60), ToString(minutes % 60),
        sun_direction.zenith * deg, sun_direction.azimuth * deg,
        "output/cache/input");
    name.push_back(ToString(minutes / 60) + "h" + ToString(minutes % 60));
    sun_zenith.push_back(atmosphere->sun_zenith());
    sun_azimuth.push_back(atmosphere->sun_azimuth());
    measurements.push_back(atmosphere);
    measured.AddAtmosphere(atmosphere);
    std::ofstream sza(
        Comparisons::GetOutputDirectory() + "sza_" + name.back() + ".txt");
    sza << round(sun_direction.zenith) << std::endl;
    sza.close();
  }

  // The Hosek model does not support wavelengths larger than 720 nm. This is
  // not an issue for luminance comparisons since the cie_y_bar_function is very
  // small in the 720-830nm range. However, for radiance and irradiance
  // comparisons with other models, it is important to integrate the spectrums
  // over the same range for all models. Thus we limit the radiance and
  // irradiance comparisons to the 360-720nm range.
  const Wavelength min_wavelength = 360.0 * nm;
  const Wavelength max_wavelength = 720.0 * nm;

  SaveLibRadtranRmseTable(libradtran_uvspec, measured, sun_zenith, sun_azimuth,
      min_wavelength, max_wavelength);
  SaveZenithLuminanceRmseTable(measured, sun_zenith, min_wavelength,
      max_wavelength);
  SavePreethamRmseTable(measured, sun_zenith, sun_azimuth, min_wavelength,
      max_wavelength);
  SaveHosekRmseTable(measured, sun_zenith, sun_azimuth, min_wavelength,
      max_wavelength);

  std::cout << std::endl << "Nishita93 model..." << std::endl;
  SaveComparisons(Comparisons("nishita93", Nishita93(), measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "Nishita96 model..." << std::endl;
  SaveComparisons(Comparisons("nishita96", Nishita96(Nishita96::ALL_ORDERS),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "Preetham model..." << std::endl;
  SaveComparisons(Comparisons("preetham", Preetham(Turbidity), measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "O'Neal model..." << std::endl;
  SaveComparisons(Comparisons("oneal", ONeal(), measured, min_wavelength,
      max_wavelength), name, sun_zenith, sun_azimuth, measurement_location,
      measurement_time, true, false);

  std::cout << std::endl << "Haber model..." << std::endl;
  SaveComparisons(Comparisons("haber", Haber(Haber::ALL_ORDERS), measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "Bruneton model..." << std::endl;
  SaveComparisons(Comparisons("bruneton", Bruneton(Bruneton::ALL_ORDERS, 3),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "Elek model..." << std::endl;
  SaveComparisons(Comparisons("elek", Bruneton(Bruneton::ALL_ORDERS, 15),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "Hosek model..." << std::endl;
  SaveComparisons(Comparisons("hosek", Hosek(Turbidity), measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, false);

  std::cout << std::endl << "libRadtran model..." << std::endl;
  SaveRadiances(Comparisons("libradtran",
      LibRadtran(libradtran_uvspec, LibRadtran::HEMISPHERICAL_FUNCTION_CACHE),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);

  SaveRadiances(Comparisons("libradtran_no_ground_albedo",
      LibRadtran(libradtran_uvspec, MieAngstromAlpha, MieAngstromBeta,
          MiePhaseFunctionG, false /* ground_albedo */,
          LibRadtran::HEMISPHERICAL_FUNCTION_CACHE),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);

  SaveComparisons(Comparisons("libradtran",
      LibRadtran(libradtran_uvspec, LibRadtran::BINARY_FUNCTION_CACHE),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, false, false);

  SaveComparisons(Comparisons("measurements", measured, measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      measurement_location, measurement_time, true, true);

  std::cout << std::endl << "Single scattering comparisons..." << std::endl;
  SaveRadiances(Comparisons("nishita96_ss",
      Nishita96(Nishita96::SINGLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("haber_ss",
      Haber(Haber::SINGLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("bruneton_ss",
      Bruneton(Bruneton::SINGLE_SCATTERING_ONLY, 3),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);

  std::cout << std::endl << "Double scattering comparisons..." << std::endl;
  SaveRadiances(Comparisons("nishita96_ds",
      Nishita96(Nishita96::DOUBLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("haber_ds",
      Haber(Haber::DOUBLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("bruneton_ds",
      Bruneton(Bruneton::DOUBLE_SCATTERING_ONLY, 3),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);

  std::cout << std::endl << "Polarization effect evaluation..." << std::endl;
  {
    PolRadtran polarized(libradtran_uvspec, libradtran_data, true);
    Comparisons polradtran_vector("polradtran_vector", polarized, measured,
        min_wavelength, max_wavelength);
    polradtran_vector.RenderLuminanceAndImage("sza30", 30.0 * deg, 90.0 * deg);

    PolRadtran non_polarized(libradtran_uvspec, libradtran_data, false);
    Comparisons polradtran_scalar("polradtran_scalar", non_polarized, measured,
        min_wavelength, max_wavelength);
    polradtran_scalar.RenderLuminanceAndImage("sza30", 30.0 * deg, 90.0 * deg);
    polradtran_scalar.PlotRelativeError("sza30", polarized, 30.0 * deg,
        90.0 * deg);
  }

  SaveSolarSpectrum();
  SaveMiePhaseFunction();

  SavePlot(measurements, name, sun_zenith, sun_azimuth);

  return 0;
}
