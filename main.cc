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
constexpr int kNumModels = 8;
constexpr int kNumViewSamples = 5;

const std::string kModels[kNumModels] = {
  "bruneton", "haber", "hosek", "libradtran", "nishita93", "nishita96",
  "preetham", "measurements"
};

const std::string kCaptions[kNumModels] = {
  "Bruneton", "Haber", "Hosek", "libRadtran", "Nishita93", "Nishita96",
  "Preetham", "Measurements"
};

const std::string kLineStyle[kNumModels] = {
  "lt 1 lc rgbcolor \"#666666\" pt 5", "lt 1 lc rgbcolor \"#aaaaaa\" pt 5",
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 12", "lt 1 lc rgbcolor \"#0000ff\" pt 4",
  "lt 1 lc rgbcolor \"#666666\" pt 9", "lt 1 lc rgbcolor \"#aaaaaa\" pt 8",
  "lt 1 lc rgbcolor \"#999999\" pt 13", "lt 1 lc rgbcolor \"#ff0000\" lw 3"
};

const std::string kSingleScatteringLineStyle[kNumModels] = {
  "lt 1 lc rgbcolor \"#666666\" pt 5", "lt 1 lc rgbcolor \"#666666\" pt 4",
  "", "", "", "lt 1 lc rgbcolor \"#666666\" pt 8", "", ""
};

const std::string kDoubleScatteringLineStyle[kNumModels] = {
  "lt 1 lc rgbcolor \"#aaaaaa\" pt 5", "lt 1 lc rgbcolor \"#aaaaaa\" pt 4",
  "", "", "", "lt 1 lc rgbcolor \"#aaaaaa\" pt 8", "", ""
};

const std::string kMeasurements = "measurements";

const double kViewZenithSamples[kNumMeasurements] = {
  90.0 - 12.1151, 90.0 - 53.3665, 90.0 - 53.3665, 90.0 - 12.1151, 90.0 - 71.9187
};

const double kViewAzimuthSamples[kNumMeasurements] = {
  360.0 - 326.25, 360.0 - 67.5, 360.0 - 225.0, 360.0 - 225.0, 360.0 - 225.0
};

const int kViewIdSamples[kNumViewSamples] = { 29, 59, 66, 20, 74 };

const int kIndices[] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
};

void SaveRadiances(const Comparisons& comparisons,
    const std::vector<std::string>& name, const std::vector<Angle>& sun_zenith,
    const std::vector<Angle>& sun_azimuth) {
  std::cout << "Computing radiance..." << std::endl;
  for (int i : kIndices) {
    for (int j = 0; j < kNumViewSamples; ++j) {
      Angle view_zenith = kViewZenithSamples[j] * deg;
      Angle view_azimuth = kViewAzimuthSamples[j] * deg;
      comparisons.PlotRadiance(
          name[i], sun_zenith[i], sun_azimuth[i], view_zenith, view_azimuth);
    }
  }
}

void SaveComparisons(const Comparisons& comparisons,
    const std::vector<std::string>& name, const std::vector<Angle>& sun_zenith,
    const std::vector<Angle>& sun_azimuth, bool save_radiances,
    bool is_measurement) {
  if (!is_measurement) {
    std::cout << "Computing sun illuminance attenuation..." << std::endl;
    comparisons.PlotSunIlluminanceAttenuation();
  }
  std::vector<Radiance> zenith_luminance;
  std::vector<Irradiance> sun_irradiance;
  std::vector<Irradiance> sky_irradiance;
  for (int i : kIndices) {
    if (save_radiances) {
      std::cout << "Computing radiance..." << std::endl;
      for (int j = 0; j < kNumViewSamples; ++j) {
        Angle view_zenith = kViewZenithSamples[j] * deg;
        Angle view_azimuth = kViewAzimuthSamples[j] * deg;
        comparisons.PlotRadiance(
            name[i], sun_zenith[i], sun_azimuth[i], view_zenith, view_azimuth);
      }
    }
    std::cout << "Computing luminance and image..." << std::endl;
    comparisons.RenderLuminanceAndImage(name[i], sun_zenith[i], sun_azimuth[i]);
    std::cout << "Computing luminance profile..." << std::endl;
    comparisons.PlotLuminanceProfile(name[i], sun_zenith[i], sun_azimuth[i]);
    if (!is_measurement) {
      std::cout << "Computing relative error..." << std::endl;
      comparisons.PlotRelativeError(name[i], sun_zenith[i], sun_azimuth[i]);
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
  comparisons.PlotDayZenithLuminance(sun_zenith, zenith_luminance);
  comparisons.PlotDayIrradiance(sun_zenith, sun_irradiance, sky_irradiance);
}

void SaveSolarSpectrum() {
  IrradianceSpectrum solar = SolarSpectrum();
  std::ofstream file("output/comparisons/solar_spectrum.txt");
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
  std::ofstream file("output/comparisons/phasefunction.txt");
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

void SaveTable(const std::string& title, const std::string& name,
    const std::vector<std::string>& names, bool include_measurements) {
  std::ofstream file("output/tables/" + name + ".html");
  file << "<html>";
  if (name == "relative_error") {
    file << "<style type=\"text/css\">";
    file << ".cell { position:relative; }";
    file << ".caption { position:absolute; bottom:0; right:0; font-size:75%; }";
    file << "</style>";
  }
  file << "<body><h1>" << title << "</h1>";
  if (name == "relative_error") {
    file << "<img src=\"../figures/error_caption.png\" "
        << "style=\"vertical-align:middle;\">[%]<p>";
  } else if (name == "relative_luminance") {
    file << "<img src=\"../figures/scale_caption.png\" "
        << "style=\"vertical-align:middle;\">[%]<p>";
  } else if (name == "absolute_luminance") {
    file << "<img src=\"../figures/scale_caption.png\" "
        << "style=\"vertical-align:middle;\">[cd/m2]<p>";
  }
  file << "<table><tr align=\"center\"><td></td>";
  for (int i = 0; i < kNumModels; ++i) {
    if (kModels[i] != kMeasurements || include_measurements) {
      file << "<td>" << kCaptions[i] << "</td>";
    }
  }
  file << "</tr>" << std::endl;
  for (int i = 0; i < kNumMeasurements; ++i) {
    file << "<tr><td>" << names[i] << "</td>" << std::endl;
    for (int j = 0; j < kNumModels; ++j) {
      if (!include_measurements && kModels[j] == kMeasurements) {
        continue;
      }
      file << "<td class=\"cell\"><img src=\"../figures/" << name << "_"
           << names[i] << "_" << kModels[j] << ".png\" width=128 height=128>";
      if (name == "relative_error") {
        std::ifstream in(
            "output/comparisons/error_" + names[i] + "_" + kModels[j] + ".txt");
        double avg, rmse;
        in >> avg >> rmse;
        in.close();
        double cvrmse = static_cast<int>(round((rmse / avg) * 100.0));
        file << "<div class=\"caption\">" << cvrmse << "%</div>";
      }
      file << "</td>" << std::endl;
    }
    file << "</tr>" << std::endl;
  }
  file << "</table>";
  if (name == "relative_error") {
    file << "<p>percentage values = Root Mean Square Error / Average";
  }
  file << "</body></html>" << std::endl;
  file.close();
}

void SaveRadianceTable(const std::string& title, const std::string& name,
    const std::vector<std::string>& names) {
  std::ofstream file("output/tables/" + name + ".html");
  file << "<html><body><h1>" << title << "</h1><table>\n";
  for (int i = 0; i < kNumMeasurements; ++i) {
    file << "<tr><td>" << names[i] << "</td>" << std::endl;
    for (int j = 0; j < kNumViewSamples; ++j) {
      Angle view_zenith = kViewZenithSamples[j] * deg;
      Angle view_azimuth = kViewAzimuthSamples[j] * deg;
      file << "<td><a href=\"../figures/" << name << "_" << names[i] << "_"
          << view_zenith.to(deg) << "_" << view_azimuth.to(deg) << ".png\">"
          << "<img src=\"../figures/" << name << "_" << names[i] << "_"
          << view_zenith.to(deg) << "_" << view_azimuth.to(deg)
          << ".png\" width=200></a></td>" << std::endl;
    }
    file << "</tr>" << std::endl;
  }
  file << "</table>";
  file << "</body></html>" << std::endl;
  file.close();
}

void SaveProfileTable(const std::vector<std::string>& names) {
  std::ofstream file("output/tables/profile.html");
  file << "<html><body><h1>Luminance profile</h1><table>\n";
  for (int i = 0; i < kNumMeasurements; ++i) {
    file << "<tr><td>" << names[i] << "</td>" << std::endl;
    file << "<td><a href=\"../figures/luminance_profile_" << names[i]
        << ".png\"><img src=\"../figures/luminance_profile_" << names[i]
        << ".png\"></a></td></tr>" << std::endl;
  }
  file << "</table>";
  file << "</body></html>" << std::endl;
  file.close();
}

std::string PolarizationCells(const std::string filename) {
  return "<td class=\"cell\"><img src=\"../figures/" + filename +
      "_sza30_polradtran_scalar.png\" width=128 height=128></td>" +
      "<td class=\"cell\"><img src=\"../figures/" + filename +
      "_sza30_polradtran_vector.png\" width=128 height=128></td>";
}

void SavePolarizationTable() {
  std::ofstream file("output/tables/polarization.html");
  file << "<html><body><h1>Effect of polarization</h1><p>"
      "<table><tr align=\"center\"><td></td><td>Without</td><td>With</td></tr>"
      "<tr><td>Rendering</td>";
  file << PolarizationCells("image");
  file << "</tr><tr><td>Rendering (non spectral)</td>";
  file << PolarizationCells("image_rgb");
  file << "</tr><tr><td>Relative error</td>";
  file << PolarizationCells("relative_error");
  file << "</tr><tr><td>Absolute luminance</td>";
  file << PolarizationCells("absolute_luminance");
  file << "</tr><tr><td>Relative luminance</td>";
  file << PolarizationCells("relative_luminance");
  file << "</tr><tr><td>Chromaticity</td>";
  file << PolarizationCells("chromaticity");
  file << "</tr></table></body></html>";
  file.close();
}

void SaveSubPlot(std::ofstream* file, Angle sun_zenith, Angle sun_azimuth,
    Angle view_zenith, Angle view_azimuth, bool profile) {
  (*file) << "set size square 0.47\n"
      << (profile ? "set origin 0.15,0.5\n" : "set origin 0.59,0.50\n")
      << "set object circle at 0.5,0.5 radius screen 0.115 behind "
      "fc rgb \"white\" fs solid\n"
      "set polar\n"
      "set angle degrees\n"
      "set grid polar 45 linestyle 1\n"
      "set xtics axis (\"\" 10, \"\" 20, \"\" 30, \"\" 40, \"\" 50, "
      "\"\" 60, \"\" 70, \"\" 80, \"\" 90) scale 0\n"
      "set xrange [-90:90]\n"
      "set yrange [-90:90]\n"
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

// Note: Fig10 a,b,c,d in Kider2014 correspond to, respectively:
// "output/figures/radiance_11h45_77.8849_33.75.png"
// "output/figures/radiance_11h45_36.6335_292.5.png"
// "output/figures/radiance_11h45_36.6335_135.png"
// "output/figures/radiance_11h45_77.8849_135.png"
void SavePlot(const std::vector<std::string>& names,
    const std::vector<Angle>& sun_zenith,
    const std::vector<Angle>& sun_azimuth,
    bool png_output) {
  std::ofstream file(png_output ? "main.plot" : "paper.plot");
  std::string out = png_output ? "output/" : "paper/";
  std::string ext = png_output ? ".png" : ".eps";
  if (png_output) {
    file << "load \"output/comparisons/error_caption.plot\"" << std::endl;
    file << "load \"output/comparisons/scale_caption.plot\"" << std::endl;
  } else {
    file << "load \"output/comparisons/error_caption_paper.plot\"" << std::endl;
    file << "load \"output/comparisons/scale_caption_paper.plot\"" << std::endl;
  }
  file << "reset" << std::endl;
  if (png_output) {
    file << "set terminal png size 680,400 enhanced" << std::endl;
  } else {
    file << "set terminal postscript eps size 6.8cm,4cm" << std::endl;
  }
  file << "set tics in nomirror scale 0.2" << std::endl;
  file << "set style line 1 lc rgbcolor \"#eeeeee\"" << std::endl;
  file << "set grid noxtics ytics linestyle 1" << std::endl << std::endl;

  file << "set output \"" << out << "figures/solar_spectrum" << ext << "\"\n";
  file << "set xlabel \"Wavelength (nm)\"" << std::endl;
  file << "set ylabel \"Spectral Irradiance (W/m2/nm)\"" << std::endl;
  file << "plot [360:830] \"input/astm-g173.txt\" with lines t \"ASTM-G173\", ";
  file << "\"output/comparisons/solar_spectrum.txt\" "
       << "t \"Our Solar Spectrum\" with lines" << std::endl << std::endl;

  file << "set output \"" << out << "figures/phase_function" << ext << "\"\n";
  file << "set xlabel \"Scattering angle (degrees)\"" << std::endl;
  file << "set ylabel \"Value (unitless)\"" << std::endl;
  file << "set logscale y" << std::endl;
  file << "plot \"output/comparisons/phasefunction.txt\" using 1:2 with lines "
       << "t \"Measured\", \"\" using 1:3 with lines t \"Cornette-Shanks\", "
       << "\"\" using 1:($3/$2) with lines t \"Ratio\"\n";
  file << "unset logscale y" << std::endl << std::endl;

  file << "set output \"" << out<< "figures/sun_illuminance_attenuation"
       << ext << "\"" << std::endl;
  file << "set xlabel \"Sun zenith angle (degrees)\"" << std::endl;
  file << "set ylabel \"Sun illuminance attenuation\"" << std::endl;
  file << "plot [0:90][0:1] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (kModels[i] == kMeasurements) {
      continue;
    } else if (i > 0) {
      file << ", ";
    }
    file << "\"output/comparisons/sun_illuminance_attenuation_" << kModels[i]
         << ".txt\" " << "t \"" << kCaptions[i] << "\" with lines";
  }
  file << std::endl << std::endl;

  if (png_output) {
    file << "set terminal png size 850,500 enhanced\n";
  } else {
    file << "set terminal postscript eps size 8.5cm,5.0cm\n";
  }
  for (int i : kIndices) {
    if (!png_output && (i % 4) != 1) {
      continue;
    }
    file << "set output \"" << out << "figures/luminance_profile_" << names[i]
         << ext << "\"" << std::endl;
      file << "reset\n"
          "set object rectangle from screen 0,0 to screen 1,1 behind "
          "fc rgb \"white\" fs solid noborder\n"
          "set multiplot\n"
          "set tics in nomirror scale 0.5\n"
          "set style line 1 lc rgbcolor \"#eeeeee\"\n"
          "set grid noxtics ytics linestyle 1\n"
          "set key outside center bottom horizontal Left reverse samplen 3 "
          "width -3 maxrows 2\n";
    if (png_output) {
      file << "set xlabel \"View zenith angle (degrees)\"\n"
          "set ylabel \"Sky luminance (cd/m^2)\"\n";
    } else {
      file << "set pointsize 0.5" << std::endl;
    }
    file << "plot [-90:90][0:25000] ";
    for (int j = 0; j < kNumModels; ++j) {
      if (j > 0) {
        file << ", ";
      }
      file << "\"output/comparisons/luminance_profile_" << names[i] << "_"
           << kModels[j] << ".txt\" t \"" << kCaptions[j] << "\"";
      if (kModels[j] == kMeasurements) {
        file << " with points lc rgbcolor \"#ff0000\" lw 3 pt 5";
      } else {
        file << " with linespoints " << kLineStyle[j];
      }
    }
    file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
        "unset border\nunset key\nunset object\n";
    SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i], 0 * deg, 135 * deg, true);
    file << "unset multiplot" << std::endl << std::endl;
  }

  for (int i : kIndices) {
    if (!png_output && i != 9) {
      continue;
    }
    for (int j = 0; j < kNumViewSamples; ++j) {
      file << "set output \"" << out << "figures/radiance_" << names[i] << "_"
           << kViewZenithSamples[j] << "_" << kViewAzimuthSamples[j]
           << ext << "\"" << std::endl;
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
      if (png_output) {
        file << "set xlabel \"Wavelength (nm)\"\n"
            "set ylabel \"Sky spectral radiance (W/m^2/sr/nm)\"\n";
      }
      file << "plot [360:720][0:0.28] ";
      for (int k = 0; k < kNumModels; ++k) {
        if (k > 0) {
          file << ", ";
        }
        if (kModels[k] == kMeasurements) {
          file << "\"input/2013-05-27___" << names[i].substr(0, 2) << "."
               << names[i].substr(3) << ".00/2013-05-27___"
               << names[i].substr(0, 2) << "." << names[i].substr(3) << ".00_-_"
               << kViewIdSamples[j] << "____"
               << ToString(360.0 - kViewAzimuthSamples[j]) << "___"
               << ToString(90.0 - kViewZenithSamples[j])
               << ".asd.rad.txt\" every ::1 t \"" << kCaptions[k]
               << "\" with lines " << kLineStyle[k];
        } else {
          file << "\"output/comparisons/radiance_" << names[i] << "_"
               << kViewZenithSamples[j] << "_" << kViewAzimuthSamples[j] << "_"
               << kModels[k] << ".txt\" t \"" << kCaptions[k] << "\""
               << " with linespoints " << kLineStyle[k];
        }
      }
      file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
          "unset border\nunset key\nunset object\n";
      SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i],
          kViewZenithSamples[j] * deg, kViewAzimuthSamples[j] * deg, false);
      file << "unset multiplot" << std::endl << std::endl;
    }
  }

  for (int i : kIndices) {
    if (!png_output && i != 9) {
      continue;
    }
    for (int j = 0; j < kNumViewSamples; ++j) {
      file << "set output \"" << out << "figures/ms_radiance_" << names[i]
           << "_" << kViewZenithSamples[j] << "_" << kViewAzimuthSamples[j]
           << ext << "\"" << std::endl;
      file << "reset\n"
          "set multiplot\n"
          "set object rectangle from screen 0,0 to screen 1,1 behind "
          "fc rgb \"white\" fs solid noborder\n"
          "set tics in nomirror scale 0.5\n"
          "set style line 1 lc rgbcolor \"#eeeeee\"\n"
          "set pointsize 1.25\n"
          "set grid noxtics ytics linestyle 1\n"
          "set key outside center bottom horizontal Left reverse maxrows 2\n";
      if (png_output) {
        file << "set xlabel \"Wavelength (nm)\"\n"
            "set ylabel \"Sky spectral radiance (W/m^2/sr/nm)\"\n";
      }
      file << "plot [360:720][0:0.28] ";
      for (int k = 0; k < kNumModels; ++k) {
        if (kModels[k] != "bruneton" && kModels[k] != "haber" &&
            kModels[k] != "nishita96") {
          continue;
        }
        if (k > 0) {
          file << ", ";
        }
        file << "\"output/comparisons/radiance_" << names[i] << "_"
             << kViewZenithSamples[j] << "_" << kViewAzimuthSamples[j] << "_"
             << kModels[k] << "_ss.txt\" t \"" << kCaptions[k] << " SS\""
             << " with linespoints " << kSingleScatteringLineStyle[k] << ", ";
        file << "\"output/comparisons/radiance_" << names[i] << "_"
             << kViewZenithSamples[j] << "_" << kViewAzimuthSamples[j] << "_"
             << kModels[k] << "_ds.txt\" t \"" << kCaptions[k] << " DS\""
             << " with linespoints " << kDoubleScatteringLineStyle[k];
      }
      file << "\nunset xlabel\nunset ylabel\nunset grid\nunset tics\n"
          "unset border\nunset key\nunset object\n";
      SaveSubPlot(&file, sun_zenith[i], sun_azimuth[i],
          kViewZenithSamples[j] * deg, kViewAzimuthSamples[j] * deg, false);
      file << "unset multiplot" << std::endl << std::endl;
    }
  }

  file << "reset\n";

  file << "set output \"" << out << "figures/day_zenith_luminance"
      << ext << "\"\n"
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
    file << "\"output/comparisons/day_zenith_luminance_" << kModels[i]
         << ".txt\" " << "t \"" << kCaptions[i] << "\" with "
         << (kModels[i] == kMeasurements ? "lines " : "linespoints ")
         << kLineStyle[i];
  }
  file << std::endl << std::endl;

  file << "set output \"" << out << "figures/day_irradiance_sky"
      << ext << "\"\n";
  if (png_output) {
    file << "set xlabel \"Time\"\n"
        "set ylabel \"Sky Irradiance [360-720nm] (W/m^2)\"\n";
  } else {
    file << "unset xlabel\n unset ylabel\n";
  }
  file << "plot [][0:] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (i > 0) {
      file << ", ";
    }
    file << "\"output/comparisons/day_irradiance_" << kModels[i]
         << ".txt\" " << "using 1:3 t \"" << kCaptions[i] << "\" with "
         << (kModels[i] == kMeasurements ? "lines " : "linespoints ")
         << kLineStyle[i];
  }
  file << std::endl << std::endl;

  file << "set output \"" << out << "figures/day_irradiance_total"
      << ext << "\"\n";
  if (png_output) {
      file << "set xlabel \"Time\"\n"
          "set ylabel \"Sun and Sky Irradiance [360-720nm] (W/m^2)\"\n";
  } else {
    file << "unset xlabel\n unset ylabel\n";
  }
  file << "plot [][0:] ";
  for (int i = 0; i < kNumModels; ++i) {
    if (kModels[i] == kMeasurements || kModels[i] == "libradtran") {
      continue;
    } else if (i > 0) {
      file << ", ";
    }
    file << "\"output/comparisons/day_irradiance_" << kModels[i]
         << ".txt\" " << "using 1:2 t \"" << kCaptions[i]
         << "\" with linespoints " << kLineStyle[i];
  }
  file << std::endl << std::endl;

  file.close();
}

}  // anonymous namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <libRatran path>" << std::endl;
    return -1;
  }
  const std::string libRadtran_path(argv[1]);

  Number mie_optical_depth = MieExtinction()(550.0 * nm) * MieScaleHeight;
  Number rayleigh_optical_depth =
      RayleighScattering()(550.0 * nm) * RayleighScaleHeight;
  Number turbidity =
      (mie_optical_depth + rayleigh_optical_depth) / rayleigh_optical_depth;
  std::cout << "turbidity: " << turbidity() << std::endl;

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
    measured.AddAtmosphere(atmosphere);
  }

  // The Hosek model does not support wavelengths larger than 720 nm. This is
  // not an issue for luminance comparisons since the cie_y_bar_function is very
  // small in the 720-830nm range. However, for radiance and irradiance
  // comparisons with other models, it is important to integrate the spectrums
  // over the same range for all models. Thus we limit the radiance and
  // irradiance comparisons to the 360-720nm range.
  const Wavelength min_wavelength = 360.0 * nm;
  const Wavelength max_wavelength = 720.0 * nm;

  std::cout << std::endl << "Bruneton model..." << std::endl;
  SaveComparisons(Comparisons("bruneton", Bruneton(Bruneton::ALL_ORDERS),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      true, false);

  std::cout << std::endl << "Haber model..." << std::endl;
  SaveComparisons(Comparisons("haber", Haber(Haber::ALL_ORDERS), measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth, true,
      false);

  std::cout << std::endl << "Hosek model..." << std::endl;
  SaveComparisons(Comparisons("hosek", Hosek(std::max(2.0, turbidity())),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      true, false);

  std::cout << std::endl << "libRadtran model..." << std::endl;
  SaveRadiances(Comparisons("libradtran", LibRadtran(libRadtran_path, false),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveComparisons(Comparisons("libradtran", LibRadtran(libRadtran_path, true),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      false, false);

  std::cout << std::endl << "Nishita93 model..." << std::endl;
  SaveComparisons(Comparisons("nishita93", Nishita93(), measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth, true,
      false);

  std::cout << std::endl << "Nishita96 model..." << std::endl;
  SaveComparisons(Comparisons("nishita96", Nishita96(Nishita96::ALL_ORDERS),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth,
      true, false);

  std::cout << std::endl << "Preetham model..." << std::endl;
  SaveComparisons(Comparisons("preetham",
      Preetham(std::max(2.0, turbidity()), false), measured, min_wavelength,
      max_wavelength), name, sun_zenith, sun_azimuth, true, false);

  SaveComparisons(Comparisons("measurements", measured, measured,
      min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth, true,
      true);

  std::cout << std::endl << "Single scattering comparisons..." << std::endl;
  SaveRadiances(Comparisons("bruneton_ss",
      Bruneton(Bruneton::SINGLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("haber_ss",
      Haber(Haber::SINGLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("nishita96_ss",
      Nishita96(Nishita96::SINGLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);

  std::cout << std::endl << "Double scattering comparisons..." << std::endl;
  SaveRadiances(Comparisons("bruneton_ds",
      Bruneton(Bruneton::DOUBLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("haber_ds",
      Haber(Haber::DOUBLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);
  SaveRadiances(Comparisons("nishita96_ds",
      Nishita96(Nishita96::DOUBLE_SCATTERING_ONLY),
      measured, min_wavelength, max_wavelength), name, sun_zenith, sun_azimuth);

  std::cout << std::endl << "Polarization effect evaluation..." << std::endl;
  {
    PolRadtran polarized(libRadtran_path, true);
    Comparisons polradtran_vector("polradtran_vector", polarized, measured,
        min_wavelength, max_wavelength);
    polradtran_vector.RenderLuminanceAndImage("sza30", 30.0 * deg, 90.0 * deg);

    PolRadtran non_polarized(libRadtran_path, false);
    Comparisons polradtran_scalar("polradtran_scalar", non_polarized, measured,
        min_wavelength, max_wavelength);
    polradtran_scalar.RenderLuminanceAndImage("sza30", 30.0 * deg, 90.0 * deg);
    polradtran_scalar.PlotRelativeError("sza30", polarized, 30.0 * deg,
        90.0 * deg);
  }

  SaveSolarSpectrum();
  SaveMiePhaseFunction();

  SaveTable("Absolute luminance", "absolute_luminance", name, true);
  SaveTable("Relative luminance", "relative_luminance", name, true);
  SaveTable("Chromaticity", "chromaticity", name, true);
  SaveTable("Relative error", "relative_error", name, false);
  SaveTable("Rendering", "image", name, true);
  SaveTable("Rendering (non spectral)", "image_rgb", name, true);
  SaveTable("Sepctral vs RGB rendering", "image_rgb_diff", name, true);

  SaveRadianceTable("Spectral Radiance samples", "radiance", name);
  SaveRadianceTable(
      "Validation: Single and double scattering", "ms_radiance", name);
  SaveProfileTable(name);
  SavePolarizationTable();

  SavePlot(name, sun_zenith, sun_azimuth, true);
  SavePlot(name, sun_zenith, sun_azimuth, false);
  Comparisons::SaveErrorCaption(true);
  Comparisons::SaveErrorCaption(false);
  Comparisons::SaveScaleCaption(true);
  Comparisons::SaveScaleCaption(false);

  return 0;
}
