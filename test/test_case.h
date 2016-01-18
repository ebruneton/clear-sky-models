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
#ifndef TEST_TEST_CASE_H_
#define TEST_TEST_CASE_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "math/scalar.h"

class TestCase {
 public:
  typedef void (TestCase::*Test)();

  TestCase(const std::string& name, Test test)
      : name_(name), test_(test), pass_(true) {
    if (!tests_) {
      tests_ = new std::vector<TestCase*>();
    }
    tests_->push_back(this);
  }

  virtual void SetUp() {}
  virtual void TearDown() {}

  void Run() {
    SetUp();
    (this->*test_)();
    TearDown();
  }

  static void RunAllTests() {
    if (tests_) {
      int pass = 0;
      int fail = 0;
      for (unsigned int i = 0; i < tests_->size(); ++i) {
        TestCase* test = (*tests_)[i];
        test->Run();
        if (test->pass_) {
          std::cout << "[  OK  ] " << test->name_ << std::endl;
          ++pass;
        } else {
          std::cout << "[ FAIL ] " << test->name_ << std::endl;
          ++fail;
        }
      }
      std::cout << (pass + fail) << " test(s) run." << std::endl;
      std::cout << pass  << " test(s) pass." << std::endl;
      std::cout << fail  << " test(s) fail." << std::endl;
    }
  }

  void ExpectTrue(bool actual) {
    if (!actual) {
      std::cout << "True expected but got False" << std::endl;
      pass_ = false;
    }
  }

  void ExpectFalse(bool actual) {
    if (actual) {
      std::cout << "False expected but got True" << std::endl;
      pass_ = false;
    }
  }

  void ExpectEquals(double expected, double actual) {
    if (!(actual == expected)) {
      std::cout << expected << " expected but got " << actual << std::endl;
      pass_ = false;
    }
  }

  template<int U1, int U2, int U3, int U4>
  void ExpectEquals(Scalar<U1, U2, U3, U4> expected,
      Scalar<U1, U2, U3, U4> actual) {
    if (!(actual == expected)) {
      std::cout << expected.to(Scalar<U1, U2, U3, U4>::Unit())
                << " expected but got "
                << actual.to(Scalar<U1, U2, U3, U4>::Unit()) << std::endl;
      pass_ = false;
    }
  }

  void ExpectNear(double expected, double actual, double tolerance) {
    if (actual - expected > tolerance || expected - actual > tolerance) {
      std::cout << expected << " +/- " << tolerance << " expected but got "
                << actual << std::endl;
      pass_ = false;
    }
  }

 private:
  std::string name_;
  Test test_;
  bool pass_;
  static std::vector<TestCase*>* tests_;
};

#endif  // TEST_TEST_CASE_H_
