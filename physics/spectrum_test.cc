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
#include "physics/spectrum.h"

#include <string>

#include "test/test_case.h"

class TestSpectrum : public TestCase {
 public:
  template<typename T>
  TestSpectrum(const std::string& name, T test)
      : TestCase("TestSpectrum " + name, static_cast<Test>(test)) {}

  void TestWavelengths() {
    WavelengthFunction<double> f;
    ExpectEquals(40, f.size());
    ExpectEquals(360.0, f.GetWavelength(0).to(nm));
    ExpectEquals(454.0, f.GetWavelength(8).to(nm));
    ExpectEquals(548.0, f.GetWavelength(16).to(nm));
    ExpectEquals(830.0, f.GetWavelength(40).to(nm));
  }

  void TestNewSpectrum() {
    std::vector<Wavelength> wavelengths;
    std::vector<SpectralIrradiance> spectral_irradiances;
    wavelengths.push_back(454.0 * nm);
    wavelengths.push_back(548.0 * nm);
    spectral_irradiances.push_back(10.0 * watt_per_square_meter_per_nm);
    spectral_irradiances.push_back(20.0 * watt_per_square_meter_per_nm);
    IrradianceSpectrum spectrum(wavelengths, spectral_irradiances);
    ExpectEquals(10.0 * watt_per_square_meter_per_nm, spectrum[0]);
    ExpectEquals(10.0 * watt_per_square_meter_per_nm, spectrum[8]);
    ExpectEquals(15.0 * watt_per_square_meter_per_nm, spectrum[12]);
    ExpectEquals(20.0 * watt_per_square_meter_per_nm, spectrum[16]);
    ExpectEquals(20.0 * watt_per_square_meter_per_nm, spectrum[39]);
  }

  void TestNewUniformSpectrum() {
    std::vector<SpectralIrradiance> values;
    values.push_back(0.0 * watt_per_square_meter_per_nm);
    values.push_back(10.0 * watt_per_square_meter_per_nm);
    values.push_back(30.0 * watt_per_square_meter_per_nm);
    values.push_back(0.0 * watt_per_square_meter_per_nm);
    IrradianceSpectrum spectrum(300 * nm, 600.0 * nm, values);
    ExpectNear(6.0, spectrum[0].to(watt_per_square_meter_per_nm), 1e-6);
    ExpectNear(15.6, spectrum[16].to(watt_per_square_meter_per_nm), 1e-6);
    ExpectEquals(0.0 * watt_per_square_meter_per_nm, spectrum[39]);
    for (unsigned int i = 0; i < spectrum.size(); ++i) {
      ExpectEquals(spectrum[i], spectrum(spectrum.GetWavelength(i)));
    }
  }

  void TestOperators() {
    WavelengthFunction<double> f;
    WavelengthFunction<Length> g;
    for (unsigned int i = 0; i < f.size(); ++i) {
      f[i] = i;
      g[i] = i * m;
    }
    f += f;
    WavelengthFunction<double> h = -f;
    WavelengthFunction<double> u = f + h;
    WavelengthFunction<double> v = u - f;
    WavelengthFunction<double> w = v * 2.0;
    WavelengthFunction<Area> x = g * (2.0 * m);
    WavelengthFunction<Area> y = g * g;
    DimensionlessSpectrum z = y / x;
    for (unsigned int i = 0; i < f.size(); ++i) {
      ExpectEquals(2.0 * i, f[i]);
      ExpectEquals(-2.0 * i, h[i]);
      ExpectEquals(0.0, u[i]);
      ExpectEquals(-2.0 * i, v[i]);
      ExpectEquals(-4.0 * i, w[i]);
      ExpectEquals(2.0 * i, x[i].to(m2));
      ExpectEquals(i * i, y[i].to(m2));
      if (i > 0) {
        ExpectEquals(i / 2.0, z[i]());
      }
    }
  }

  void TestFunctions() {
    WavelengthFunction<double> f, g;
    for (unsigned int i = 0; i < f.size(); ++i) {
      f[i] = i;
      g[i] = f.size() - 1 - i;
    }
    WavelengthFunction<double> u = exp(f * 1e-2);
    WavelengthFunction<double> v = min(f, g);
    for (unsigned int i = 0; i < f.size(); ++i) {
      ExpectEquals(exp(i * 1e-2), u[i]);
      ExpectEquals(std::min(i, f.size() - 1 - i), v[i]);
    }
  }

  void TestIntegral() {
    IrradianceSpectrum spectrum;
    for (unsigned int i = 0; i < spectrum.size(); ++i) {
      spectrum[i] = (i < 8 ? 1.0 : 0.0) * watt_per_square_meter_per_nm;
    }
    ExpectEquals(94.0 * watt_per_square_meter, Integral(spectrum));
  }
};

namespace {

TestSpectrum wavelengths("wavelengths", &TestSpectrum::TestWavelengths);
TestSpectrum newspectrum("newspectrum", &TestSpectrum::TestNewSpectrum);
TestSpectrum newuniformspectrum(
    "newuniformspectrum", &TestSpectrum::TestNewUniformSpectrum);
TestSpectrum operators("operators", &TestSpectrum::TestOperators);
TestSpectrum functions("functions", &TestSpectrum::TestFunctions);
TestSpectrum integral("integral", &TestSpectrum::TestIntegral);

}  // anonymous namespace
