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
#include "math/function.h"

#include <string>

#include "test/test_case.h"

class TestFunction : public TestCase {
 public:
  template<typename T>
  TestFunction(const std::string& name, T test)
      : TestCase("TestFunction " + name, static_cast<Test>(test)) {}

  void TestConstant() {
    Function<10, double> f(1.0);
    ExpectEquals(10, f.size());
    for (unsigned int i = 0; i < f.size(); ++i) {
      ExpectEquals(1.0, f[i]);
    }
    ExpectEquals(1.0, f(-1.0));
    ExpectEquals(1.0, f(0.0));
    ExpectEquals(1.0, f(0.333));
    ExpectEquals(1.0, f(1.0));
    ExpectEquals(1.0, f(2.0));
  }

  void TestInterpolation() {
    Function<10, double> f;
    for (unsigned int i = 0.0; i < f.size(); ++i) {
      f[i] = i;
    }
    Function<10, double> g = f;
    ExpectNear(0.0, g(0.0), 1e-9);
    ExpectNear(0.0, g(0.05), 1e-9);
    ExpectNear(0.5, g(0.1), 1e-9);
    ExpectNear(9.0, g(0.95), 1e-9);
    ExpectNear(9.0, g(1.0), 1e-9);
  }
};

namespace {

TestFunction constant("constant", &TestFunction::TestConstant);
TestFunction interpolation("interpolation", &TestFunction::TestInterpolation);

}  // anonymous namespace
