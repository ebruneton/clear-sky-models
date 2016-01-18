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
#include "atmosphere/atmosphere.h"

#include <string>

#include "test/test_case.h"

class UniformSky : public Atmosphere {
 public:
  IrradianceSpectrum GetSunIrradiance(Length altitude,
      Angle sun_zenith) const {
    return IrradianceSpectrum(0.0 * watt_per_square_meter_per_nm);
  }

  RadianceSpectrum GetSkyRadiance(Length altitude, Angle sun_zenith,
      Angle view_zenith, Angle view_sun) const {
    return RadianceSpectrum(1.0 * watt_per_square_meter_per_sr_per_nm);
  }
};

class TestAtmosphere : public TestCase {
 public:
  template<typename T>
  TestAtmosphere(const std::string& name, T test)
      : TestCase("TestAtmosphere " + name, static_cast<Test>(test)) {}

  void TestSolarSpectrum() {
    ExpectNear(734.0, Integral(SolarSpectrum()).to(watt_per_square_meter), 1.0);
  }

  void TestPhaseFunctions() {
    // Check that the phase functions integrated over 4.pi sr give 1.
    Number rayleigh_integral = 0.0;
    Number mie_integral = 0.0;
    const unsigned int N = 100;
    for (unsigned int i = 0; i < N; ++i) {
      Angle theta = (i + 0.5) * pi / N;
      SolidAngle domega = sin(theta) * (PI / N) * (2.0 * PI) * sr;
      rayleigh_integral =
          rayleigh_integral + RayleighPhaseFunction(theta) * domega;
      mie_integral =
          mie_integral + MiePhaseFunction(theta) * domega;
    }
    ExpectNear(1.0, rayleigh_integral(), 1e-3);
    ExpectNear(1.0, mie_integral(), 1e-3);
  }

  void TestSkyIrradiance() {
    // Check that the irradiance of a sky with uniform radiance 1 is pi.
    UniformSky sky;
    IrradianceSpectrum spectrum = sky.GetSkyIrradiance(0.0 * m, 0.0 * rad);
    for (unsigned int i = 0; i < spectrum.size(); ++i) {
      ExpectNear(PI, spectrum[i].to(watt_per_square_meter_per_nm), 1e-3);
    }
  }
};

namespace {

TestAtmosphere solarspectrum(
    "solarspectrum", &TestAtmosphere::TestSolarSpectrum);
TestAtmosphere phasefunctions(
    "phasefunctions", &TestAtmosphere::TestPhaseFunctions);
TestAtmosphere skyirradiance(
    "skyirradiance", &TestAtmosphere::TestSkyIrradiance);

}  // anonymous namespace
