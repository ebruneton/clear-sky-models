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
#include "math/ternary_function.h"

#include <string>

#include "test/test_case.h"

class TestTernaryFunction : public TestCase {
 public:
  template<typename T>
  TestTernaryFunction(const std::string& name, T test)
      : TestCase("TestTernaryFunction " + name, static_cast<Test>(test)) {}

  void TestConstant() {
    TernaryFunction<10, 20, 5, double> f(1.0);
    ExpectEquals(10, f.size_x());
    ExpectEquals(20, f.size_y());
    ExpectEquals(5, f.size_z());
    for (unsigned int k = 0; k < f.size_z(); ++k) {
      for (unsigned int j = 0; j < f.size_y(); ++j) {
        for (unsigned int i = 0; i < f.size_x(); ++i) {
          ExpectEquals(1.0, f.Get(i, j, k));
        }
      }
    }
    ExpectEquals(1.0, f(-1.0, -1.0, 0.5));
    ExpectEquals(1.0, f(0.0, 0.0, 0.0));
    ExpectEquals(1.0, f(0.333, 0.666, 0.999));
    ExpectEquals(1.0, f(1.0, 1.0, 1.0));
    ExpectEquals(1.0, f(2.0, 2.0, 2.0));
    f += f;
    ExpectEquals(2.0, f(-1.0, -1.0, 0.5));
    ExpectEquals(2.0, f(0.0, 0.0, 0.0));
    ExpectEquals(2.0, f(0.333, 0.666, 0.999));
    ExpectEquals(2.0, f(1.0, 1.0, 1.0));
    ExpectEquals(2.0, f(2.0, 2.0, 2.0));
  }

  void TestInterpolation() {
    TernaryFunction<4, 8, 16, double> f;
    for (unsigned int k = 0; k < f.size_z(); ++k) {
      for (unsigned int j = 0; j < f.size_y(); ++j) {
        for (unsigned int i = 0; i < f.size_x(); ++i) {
          f.Set(i, j, k, i + 2.0 * j + 4.0 * k);
        }
      }
    }
    ExpectNear(0.0, f(0.5 / 4.0, 0.5 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(1.0, f(1.5 / 4.0, 0.5 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(2.0, f(0.5 / 4.0, 1.5 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(3.0, f(1.5 / 4.0, 1.5 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(4.0, f(0.5 / 4.0, 0.5 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(5.0, f(1.5 / 4.0, 0.5 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(6.0, f(0.5 / 4.0, 1.5 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(7.0, f(1.5 / 4.0, 1.5 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(0.75, f(1.25 / 4.0, 0.5 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(1.25, f(1.25 / 4.0, 0.75 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(2.75, f(1.25 / 4.0, 1.5 / 8.0, 0.5 / 16.0), 1e-9);
    ExpectNear(4.75, f(1.25 / 4.0, 0.5 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(5.25, f(1.25 / 4.0, 0.75 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(6.75, f(1.25 / 4.0, 1.5 / 8.0, 1.5 / 16.0), 1e-9);
    ExpectNear(3.25, f(1.25 / 4.0, 0.75 / 8.0, 1.0 / 16.0), 1e-9);
  }

  void TestLoadSave() {
    TernaryFunction<4, 8, 16, double> f;
    for (unsigned int k = 0; k < f.size_z(); ++k) {
      for (unsigned int j = 0; j < f.size_y(); ++j) {
        for (unsigned int i = 0; i < f.size_x(); ++i) {
          f.Set(i, j, k, i + 2.0 * j + 4.0 * k);
        }
      }
    }
    f.Save("output/Debug/math/ternary_fuction_test.dat");

    f = TernaryFunction<4, 8, 16, double>(0.0);

    f.Load("output/Debug/math/ternary_fuction_test.dat");
    for (unsigned int k = 0; k < f.size_z(); ++k) {
      for (unsigned int j = 0; j < f.size_y(); ++j) {
        for (unsigned int i = 0; i < f.size_x(); ++i) {
          ExpectEquals(i + 2.0 * j + 4.0 * k, f.Get(i, j, k));
        }
      }
    }
  }
};

namespace {

TestTernaryFunction constant("constant", &TestTernaryFunction::TestConstant);
TestTernaryFunction interpolation(
    "interpolation", &TestTernaryFunction::TestInterpolation);
TestTernaryFunction loadsave("loadsave", &TestTernaryFunction::TestLoadSave);

}  // anonymous namespace
