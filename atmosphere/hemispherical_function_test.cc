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
#include "atmosphere/hemispherical_function.h"

#include <string>

#include "test/test_case.h"

class TestHemisphericalFunction : public dimensional::TestCase {
 public:
  template<typename T>
  TestHemisphericalFunction(const std::string& name, T test)
      : TestCase("TestHemisphericalFunction " + name, static_cast<Test>(test)) {
  }

  void TestGetSet() {
    HemisphericalFunction<double> f;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        f.Set(i, j, i * 10 + j);
      }
    }
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        ExpectEquals(i * 10 + j, f.Get(i, j));
      }
    }
  }

  void TestGetSampleDirection() {
    Angle view_zenith;
    Angle view_azimuth;
    HemisphericalFunction<double>::GetSampleDirection(
        4, 4, &view_zenith, &view_azimuth);
    ExpectNear(0.0, view_zenith.to(deg), 1e-9);

    HemisphericalFunction<double>::GetSampleDirection(
        8, 4, &view_zenith, &view_azimuth);
    ExpectNear(90.0 - 12.1151,  view_zenith.to(deg), 1e-4);
    ExpectNear(0.0, view_azimuth.to(deg), 1e-4);

    HemisphericalFunction<double>::GetSampleDirection(
        4, 7, &view_zenith, &view_azimuth);
    ExpectNear(90.0 - 33.749,  view_zenith.to(deg), 1e-4);
    ExpectNear(90.0, view_azimuth.to(deg), 1e-4);

    HemisphericalFunction<double>::GetSampleDirection(
        2, 4, &view_zenith, &view_azimuth);
    ExpectNear(90.0 - 53.3665,  view_zenith.to(deg), 1e-4);
    ExpectNear(180.0, view_azimuth.to(deg), 1e-4);

    HemisphericalFunction<double>::GetSampleDirection(
        4, 3, &view_zenith, &view_azimuth);
    ExpectNear(90.0 - 71.9187,  view_zenith.to(deg), 1e-4);
    ExpectNear(270.0, view_azimuth.to(deg), 1e-4);
  }

  void TestInterpolation() {
    HemisphericalFunction<double> f;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        f.Set(i, j, i * 10 + j);
      }
    }
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        Angle view_zenith;
        Angle view_azimuth;
        HemisphericalFunction<double>::GetSampleDirection(
            i, j, &view_zenith, &view_azimuth);
        ExpectNear(i * 10 + j, f(view_zenith, view_azimuth), 1e-8);
        ExpectNear(i * 10 + j, f(view_zenith, view_azimuth - 2 * pi), 1e-8);
        ExpectNear(i * 10 + j, f(view_zenith, view_azimuth + 2 * pi), 1e-8);
      }
    }
    ExpectNear(84.5, f((90.0 - 12.1151) * deg, (90.0 / 16.0) * deg), 1e-1);
    ExpectNear(
         69.0, f((90.0 - (53.3665 + 33.749) * 0.5) * deg, 0.0 * deg), 1e-1);
  }

  void TestLoadSave() {
    HemisphericalFunction<double> f;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        f.Set(i, j, i * 10 + j);
      }
    }
    f.Save("/tmp/clear_sky_models_atmosphere_hemispherical_fuction_test.dat");

    f = HemisphericalFunction<double>();

    f.Load("/tmp/clear_sky_models_atmosphere_hemispherical_fuction_test.dat");
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 9; ++j) {
        ExpectEquals(i * 10 + j, f.Get(i, j));
      }
    }
  }
};

namespace {

TestHemisphericalFunction getset(
    "getset", &TestHemisphericalFunction::TestGetSet);
TestHemisphericalFunction getsampledirection(
    "getsampledirection", &TestHemisphericalFunction::TestGetSampleDirection);
TestHemisphericalFunction interpolation(
    "interpolation", &TestHemisphericalFunction::TestInterpolation);
TestHemisphericalFunction loadsave(
    "loadsave", &TestHemisphericalFunction::TestLoadSave);

}  // anonymous namespace

