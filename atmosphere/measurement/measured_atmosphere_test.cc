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

#include <string>

#include "test/test_case.h"

class TestMeasuredAtmosphere : public dimensional::TestCase {
 public:
  template<typename T>
  TestMeasuredAtmosphere(const std::string& name, T test)
      : TestCase("TestMeasuredAtmosphere " + name, static_cast<Test>(test)) {}

  void TestLoadKiderData() {
    MeasuredAtmosphere atmosphere(
        "atmosphere/measurement/testdata", "2013-05-27", "11", "45",
        1.0 * deg, 2.0* deg, "", false /* compute_azimuth_from_data */);
    SpectralRadiance sr1 = atmosphere.GetSkyRadiance(0 * m, 1 * deg, 2 * deg,
        (90.0 - 12.1151) * deg, (360.0 - 326.25) * deg)(400.0 * nm);
    SpectralRadiance sr2 = atmosphere.GetSkyRadiance(0 * m, 1 * deg, 2 * deg,
        (90.0 - 53.3665) * deg, (360.0 - 67.5) * deg)(400.0 * nm);
    SpectralRadiance sr3 = atmosphere.GetSkyRadiance(0 * m, 1 * deg, 2 * deg,
        (90.0 - 53.3665) * deg, (360.0 - 225.0) * deg)(400.0 * nm);
    SpectralRadiance sr4 = atmosphere.GetSkyRadiance(0 * m, 1 * deg, 2 * deg,
        (90.0 - 12.1151) * deg, (360.0 - 225.0) * deg)(400.0 * nm);
    ExpectNear(29, sr1.to(watt_per_square_meter_per_sr_per_nm), 1e-3);
    ExpectNear(59, sr2.to(watt_per_square_meter_per_sr_per_nm), 1e-3);
    ExpectNear(66, sr3.to(watt_per_square_meter_per_sr_per_nm), 1e-3);
    ExpectNear(20, sr4.to(watt_per_square_meter_per_sr_per_nm), 1e-3);
  }
};

namespace {

TestMeasuredAtmosphere loadkiderdata(
    "loadkiderdata", &TestMeasuredAtmosphere::TestLoadKiderData);

}  // anonymous namespace
