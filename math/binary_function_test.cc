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
#include "math/binary_function.h"

#include <string>

#include "test/test_case.h"

class TestBinaryFunction : public TestCase {
 public:
  template<typename T>
  TestBinaryFunction(const std::string& name, T test)
      : TestCase("TestBinaryFunction " + name, static_cast<Test>(test)) {}

  void TestConstant() {
    BinaryFunction<10, 20, double> f(1.0);
    ExpectEquals(10, f.size_x());
    ExpectEquals(20, f.size_y());
    for (unsigned int j = 0; j < f.size_y(); ++j) {
      for (unsigned int i = 0; i < f.size_x(); ++i) {
        ExpectEquals(1.0, f.Get(i, j));
      }
    }
    ExpectEquals(1.0, f(-1.0, -1.0));
    ExpectEquals(1.0, f(0.0, 0.0));
    ExpectEquals(1.0, f(0.333, 0.666));
    ExpectEquals(1.0, f(1.0, 1.0));
    ExpectEquals(1.0, f(2.0, 2.0));
    f += f;
    ExpectEquals(2.0, f(-1.0, -1.0));
    ExpectEquals(2.0, f(0.0, 0.0));
    ExpectEquals(2.0, f(0.333, 0.666));
    ExpectEquals(2.0, f(1.0, 1.0));
    ExpectEquals(2.0, f(2.0, 2.0));
  }

  void TestInterpolation() {
    BinaryFunction<4, 8, double> f;
    for (unsigned int j = 0; j < f.size_y(); ++j) {
      for (unsigned int i = 0; i < f.size_x(); ++i) {
        f.Set(i, j, i + 2.0 * j);
      }
    }
    ExpectNear(0.0, f(0.5 / 4.0, 0.5 / 8.0), 1e-9);
    ExpectNear(1.0, f(1.5 / 4.0, 0.5 / 8.0), 1e-9);
    ExpectNear(2.0, f(0.5 / 4.0, 1.5 / 8.0), 1e-9);
    ExpectNear(3.0, f(1.5 / 4.0, 1.5 / 8.0), 1e-9);
    ExpectNear(0.75, f(1.25 / 4.0, 0.5 / 8.0), 1e-9);
    ExpectNear(1.25, f(1.25 / 4.0, 0.75 / 8.0), 1e-9);
    ExpectNear(2.75, f(1.25 / 4.0, 1.5 / 8.0), 1e-9);
  }

  void TestLoadSave() {
    BinaryFunction<4, 8, double> f;
    for (unsigned int j = 0; j < f.size_y(); ++j) {
      for (unsigned int i = 0; i < f.size_x(); ++i) {
        f.Set(i, j, i + 2.0 * j);
      }
    }
    f.Save("output/Debug/math/binary_fuction_test.dat");

    f = BinaryFunction<4, 8, double>(0.0);

    f.Load("output/Debug/math/binary_fuction_test.dat");
    for (unsigned int j = 0; j < f.size_y(); ++j) {
      for (unsigned int i = 0; i < f.size_x(); ++i) {
        ExpectEquals(i + 2.0 * j, f.Get(i, j));
      }
    }
  }
};

namespace {

TestBinaryFunction constant("constant", &TestBinaryFunction::TestConstant);
TestBinaryFunction interpolation(
    "interpolation", &TestBinaryFunction::TestInterpolation);
TestBinaryFunction loadsave("loadsave", &TestBinaryFunction::TestLoadSave);

}  // anonymous namespace
