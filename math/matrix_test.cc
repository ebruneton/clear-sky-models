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
#include "math/matrix.h"

#include <string>

#include "test/test_case.h"

class TestMatrix : public TestCase {
 public:
  template<typename T>
  TestMatrix(const std::string& name, T test)
      : TestCase("TestMatrix " + name, static_cast<Test>(test)) {}

  void TestVectorProduct() {
    Scalar<1, 2, 3, 4, 5> u = Scalar<1, 2, 3, 4, 5>::Unit();
    Scalar<10, 20, 30, 40, 50> t = Scalar<10, 20, 30, 40, 50>::Unit();
    Matrix3<Scalar<1, 2, 3, 4, 5>> m(
        1.0 * u, 2.0 * u, 3.0 * u,
        4.0 * u, 5.0 * u, 6.0 * u,
        7.0 * u, 8.0 * u, 9.0 * u);
    Vector3<Scalar<10, 20, 30, 40, 50>> x(1.0 * t, 2.0 * t, 3.0 * t);
    Vector3<Scalar<11, 22, 33, 44, 55>> y = m * x;

    const double kEps = 1e-12;
    ExpectNear(14.0, y.x.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(32.0, y.y.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(50.0, y.z.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
  }

  void TestMatrixProduct() {
    Scalar<1, 2, 3, 4, 5> u = Scalar<1, 2, 3, 4, 5>::Unit();
    Scalar<10, 20, 30, 40, 50> t = Scalar<10, 20, 30, 40, 50>::Unit();
    Matrix3<Scalar<1, 2, 3, 4, 5>> m(
        1.0 * u, 2.0 * u, 3.0 * u,
        4.0 * u, 5.0 * u, 6.0 * u,
        7.0 * u, 8.0 * u, 9.0 * u);
    Matrix3<Scalar<10, 20, 30, 40, 50>> n(
        10.0 * t, 20.0 * t, 30.0 * t,
        40.0 * t, 50.0 * t, 60.0 * t,
        70.0 * t, 80.0 * t, 90.0 * t);
    Matrix3<Scalar<11, 22, 33, 44, 55>> p = m * n;

    const double kEps = 1e-12;
    ExpectNear(300.0, p.m00.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(360.0, p.m01.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(420.0, p.m02.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(660.0, p.m10.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(810.0, p.m11.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(960.0, p.m12.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(1020.0, p.m20.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(1260.0, p.m21.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
    ExpectNear(1500.0, p.m22.to(Scalar<11, 22, 33, 44, 55>::Unit()), kEps);
  }

  void TestMatrixInverse() {
    Scalar<1, 2, 3, 4, 5> u = Scalar<1, 2, 3, 4, 5>::Unit();
    Matrix3<Scalar<1, 2, 3, 4, 5>> m(
        +1.0 * u, +1.0 * u, 0.0 * u,
        -2.0 * u, -1.0 * u, 0.0 * u,
        +3.0 * u, -2.0 * u, 1.0 * u);
    Matrix3<Scalar<-1, -2, -3, -4, -5>> p = inverse(m);

    const double kEps = 1e-12;
    ExpectNear(-1.0, p.m00.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(-1.0, p.m01.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(0.0, p.m02.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(2.0, p.m10.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(1.0, p.m11.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(0.0, p.m12.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(7.0, p.m20.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(5.0, p.m21.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
    ExpectNear(1.0, p.m22.to(Scalar<-1, -2, -3, -4, -5>::Unit()), kEps);
  }
};

namespace {

TestMatrix vectorproduct("vectorproduct", &TestMatrix::TestVectorProduct);
TestMatrix matrixproduct("matrixproduct", &TestMatrix::TestMatrixProduct);
TestMatrix matrixinverse("matrixinverse", &TestMatrix::TestMatrixInverse);

}  // anonymous namespace
