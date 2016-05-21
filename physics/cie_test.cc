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
#include "physics/cie.h"

#include <string>

#include "test/test_case.h"

class TestCie : public TestCase {
 public:
  template<typename T>
  TestCie(const std::string& name, T test)
      : TestCase("TestCie " + name, static_cast<Test>(test)) {}

  void TestColorFunctions() {
    const double kEps = 1e-12;
    ExpectNear(0.0001299, cie_x_bar_function()(360.0 * nm)(), kEps);
    ExpectNear(0.000003917, cie_y_bar_function()(360.0 * nm)(), kEps);
    ExpectNear(0.0006061, cie_z_bar_function()(360.0 * nm)(), kEps);
    ExpectNear(1.0567, cie_x_bar_function()(595.0 * nm)(), kEps);
    ExpectNear(0.6949, cie_y_bar_function()(595.0 * nm)(), kEps);
    ExpectNear(0.001, cie_z_bar_function()(595.0 * nm)(), kEps);
  }

  void TestSFunctions() {
    const double kEps = 1e-12;
    ExpectNear(61.5, S0_function()(360.0 * nm)(), kEps);
    ExpectNear(38.0, S1_function()(360.0 * nm)(), kEps);
    ExpectNear(5.3, S2_function()(360.0 * nm)(), kEps);
    ExpectNear(89.8, S0_function()(595.0 * nm)(), kEps);
    ExpectNear(-4.65, S1_function()(595.0 * nm)(), kEps);
    ExpectNear(2.65, S2_function()(595.0 * nm)(), kEps);
  }
};

namespace {

TestCie colorfunctions("colorfunctions", &TestCie::TestColorFunctions);
TestCie sfunctions("sfunctions", &TestCie::TestSFunctions);

}  // anonymous namespace

