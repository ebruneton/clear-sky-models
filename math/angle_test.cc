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
#include "math/angle.h"

#include <string>

#include "test/test_case.h"

class TestAngle : public TestCase {
 public:
  template<typename T>
  TestAngle(const std::string& name, T test)
      : TestCase("TestAngle " + name, static_cast<Test>(test)) {}

  void TestOperators() {
    Angle a = Angle::Unit();
    Angle b = 2.0 * a;
    Angle c = a + b;
    Angle d = c - b;
    Angle e = d * 3.0;
    Angle f = e / 2.0;
    Number g = -e / b;
    ExpectEquals(2.0, b.to(a));
    ExpectEquals(3.0, c.to(a));
    ExpectEquals(1.0, d.to(a));
    ExpectEquals(3.0, e.to(a));
    ExpectEquals(1.5, f.to(a));
    ExpectEquals(-1.5, g());
    ExpectTrue(e == c);
    ExpectFalse(e == a);
    ExpectTrue(b < c);
    ExpectFalse(b < b);
    ExpectFalse(c < b);
    ExpectTrue(b <= c);
    ExpectTrue(b <= b);
    ExpectFalse(c <= b);
    ExpectTrue(c > b);
    ExpectFalse(c > c);
    ExpectFalse(b > c);
    ExpectTrue(c >= b);
    ExpectTrue(c >= c);
    ExpectFalse(b >= c);
  }

  void TestFunctions() {
    ExpectEquals(std::acos(1.0), acos(Number(1.0)).to(rad));
    ExpectEquals(std::asin(1.0), asin(Number(1.0)).to(rad));
    ExpectEquals(std::atan(2.0), atan(Number(2.0)).to(rad));
    ExpectEquals(std::atan2(1.0, 2.0), atan2(Number(1.0), Number(2.0)).to(rad));
    ExpectEquals(Number(std::cos(1.0)), cos(1.0 * rad));
    ExpectEquals(Number(std::sin(1.0)), sin(1.0 * rad));
    ExpectEquals(Number(std::tan(1.0)), tan(1.0 * rad));
  }
};

namespace {

TestAngle operators("operators", &TestAngle::TestOperators);
TestAngle functions("functions", &TestAngle::TestFunctions);

}  // anonymous namespace
