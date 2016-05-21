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
#ifndef PHYSICS_SPECTRUM_H_
#define PHYSICS_SPECTRUM_H_

#include <algorithm>
#include <cassert>
#include <vector>

#include "math/scalar.h"
#include "physics/units.h"

namespace spectrum {

constexpr unsigned int NUM_WAVELENGTH = 40;
constexpr Wavelength MIN_WAVELENGTH = 360.0 * nm;
constexpr Wavelength MAX_WAVELENGTH = 830.0 * nm;

}  // namespace spectrum

// A function of wavelength from spectrum::MIN_WAVELENGTH to
// spectrum::MAX_WAVELENGTH, to values of type T, represented
// by its values at spectrum::NUM_WAVELENGTH uniformly distributed wavelengths,
// and linearly interpolated in between.
template<typename T>
class WavelengthFunction {
 public:
  WavelengthFunction() {}

  explicit WavelengthFunction(const T& constant_value) {
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      value_[i] = constant_value;
    }
  }

  WavelengthFunction(const WavelengthFunction& rhs) {
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      value_[i] = rhs.value_[i];
    }
  }

  // Creates a new spectrum from the given values, sampled at the given sampling
  // points. The given vectors must have the same size n and must be sorted in
  // increasing order of wavelength. The input function is supposed constant
  // outside the sampling interval (i.e. for lambda < sampling_points[0] and for
  // lambda > sampling_points[n - 1]) and is linearly interpolated between
  // samples.
  WavelengthFunction(const std::vector<Wavelength>& sampling_points,
      const std::vector<T>& sampled_values) {
    Init(sampling_points, sampled_values);
  }

  // Creates a new spectrum from the given values, uniformly sampled from
  // 'min_wavelength' and 'max_wavelength' wavelengths (inclusive). The input
  // function is supposed constant outside the sampling interval (i.e. for
  // lambda < min_wavelength and for lambda > max_wavelength) and is linearly
  // interpolated between samples.
  WavelengthFunction(Wavelength min_wavelength, Wavelength max_wavelength,
      const std::vector<T> values) {
    std::vector<Wavelength> sampling_points;
    for (unsigned int i = 0; i < values.size(); ++i) {
      sampling_points.push_back(
          lerp(min_wavelength, max_wavelength, i / (values.size() - 1.0)));
    }
    Init(sampling_points, values);
  }

  WavelengthFunction& operator=(const WavelengthFunction& rhs) {
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      value_[i] = rhs.value_[i];
    }
    return *this;
  }

  WavelengthFunction& operator+=(const WavelengthFunction& rhs) {
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      value_[i] += rhs.value_[i];
    }
    return *this;
  }

  WavelengthFunction operator+(const WavelengthFunction& rhs) const {
    WavelengthFunction result;
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      result.value_[i] = value_[i] + rhs.value_[i];
    }
    return result;
  }

  WavelengthFunction operator-() const {
    WavelengthFunction result;
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      result.value_[i] = -value_[i];
    }
    return result;
  }

  WavelengthFunction operator-(const WavelengthFunction& rhs) const {
    WavelengthFunction result;
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      result.value_[i] = value_[i] - rhs.value_[i];
    }
    return result;
  }

  WavelengthFunction operator*(double rhs) const {
    WavelengthFunction result;
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      result.value_[i] = value_[i] * rhs;
    }
    return result;
  }

  inline unsigned int size() const {
    return spectrum::NUM_WAVELENGTH;
  }

  inline const T operator[](int index) const {
    return value_[index];
  }

  inline T& operator[](int index) {
    return value_[index];
  }

  inline Wavelength GetWavelength(int wavelength_index) const {
    double u = static_cast<double>(wavelength_index) / spectrum::NUM_WAVELENGTH;
    return lerp(spectrum::MIN_WAVELENGTH, spectrum::MAX_WAVELENGTH, u);
  }

  const T operator()(Wavelength lambda) const {
    Number x = (lambda - spectrum::MIN_WAVELENGTH) /
        (spectrum::MAX_WAVELENGTH - spectrum::MIN_WAVELENGTH);
    double u = x() * spectrum::NUM_WAVELENGTH;
    int i = std::floor(u);
    u -= i;
    int i0 = std::max(0,
        std::min(static_cast<int>(spectrum::NUM_WAVELENGTH) - 1, i));
    int i1 = std::max(0,
        std::min(static_cast<int>(spectrum::NUM_WAVELENGTH) - 1, i + 1));
    return value_[i0] * (1.0 - u) + value_[i1] * u;
  }

 private:
  void Init(const std::vector<Wavelength>& sampling_points,
      const std::vector<T>& sampled_values) {
    assert(sampling_points.size() == sampled_values.size());
    for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
      Wavelength lambda = GetWavelength(i);
      // Find the two nearest sampling points around 'lambda' using binary
      // search and return the corresponding sampled values, linearly
      // interpolated.
      int min_index = 0;
      int max_index = sampling_points.size() - 1;
      if (lambda <= sampling_points[min_index]) {
        value_[i] = sampled_values[min_index];
        continue;
      }
      if (lambda >= sampling_points[max_index]) {
        value_[i] = sampled_values[max_index];
        continue;
      }
      while (max_index - min_index > 1) {
        int mid_index = (min_index + max_index) / 2;
        if (lambda < sampling_points[mid_index]) {
          max_index = mid_index;
        } else {
          min_index = mid_index;
        }
      }
      Number u = (lambda - sampling_points[min_index]) /
          (sampling_points[max_index] - sampling_points[min_index]);
      value_[i] = sampled_values[min_index] * (1.0 - u()) +
          sampled_values[max_index] * u();
    }
  }

  T value_[spectrum::NUM_WAVELENGTH];
};

// A function from wavelength to spectral irradiance values.
typedef WavelengthFunction<Dimensionless> DimensionlessSpectrum;

// A function from wavelength to spectral power values.
typedef WavelengthFunction<SpectralPower> PowerSpectrum;

// A function from wavelength to spectral irradiance values.
typedef WavelengthFunction<SpectralIrradiance> IrradianceSpectrum;

// A function from wavelength to spectral radiance values.
typedef WavelengthFunction<SpectralRadiance> RadianceSpectrum;

// A function from wavelength to scattering coefficient values.
typedef WavelengthFunction<ScatteringCoefficient> ScatteringSpectrum;

// A function from wavelength to phase function values.
typedef WavelengthFunction<InverseSolidAngle> PhaseFunctionSpectrum;

template<int L1, int WL1, int SA1, int P1, int LP1,
    int L2, int WL2, int SA2, int P2, int LP2>
WavelengthFunction<Scalar<L1 + L2, WL1 + WL2, SA1 + SA2, P1 + P2, LP1 + LP2>>
operator*(
    WavelengthFunction<Scalar<L1, WL1, SA1, P1, LP1>> lhs,
    Scalar<L2, WL2, SA2, P2, LP2> rhs) {
  WavelengthFunction<Scalar<L1 + L2, WL1 + WL2, SA1 + SA2, P1 + P2, LP1 + LP2>>
      result;
  for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
    result[i] = lhs[i] * rhs;
  }
  return result;
}

template<int L1, int WL1, int SA1, int P1, int LP1,
    int L2, int WL2, int SA2, int P2, int LP2>
WavelengthFunction<Scalar<L1 + L2, WL1 + WL2, SA1 + SA2, P1 + P2, LP1 + LP2>>
operator*(
    WavelengthFunction<Scalar<L1, WL1, SA1, P1, LP1>> lhs,
    WavelengthFunction<Scalar<L2, WL2, SA2, P2, LP2>> rhs) {
  WavelengthFunction<Scalar<L1 + L2, WL1 + WL2, SA1 + SA2, P1 + P2, LP1 + LP2>>
      result;
  for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
    result[i] = lhs[i] * rhs[i];
  }
  return result;
}

template<int L1, int WL1, int SA1, int P1, int LP1,
    int L2, int WL2, int SA2, int P2, int LP2>
WavelengthFunction<Scalar<L1 - L2, WL1 - WL2, SA1 - SA2, P1 - P2, LP1 - LP2>>
operator/(
    WavelengthFunction<Scalar<L1, WL1, SA1, P1, LP1>> lhs,
    WavelengthFunction<Scalar<L2, WL2, SA2, P2, LP2>> rhs) {
  WavelengthFunction<Scalar<L1 - L2, WL1 - WL2, SA1 - SA2, P1 - P2, LP1 - LP2>>
      result;
  for (unsigned int i = 0; i < spectrum::NUM_WAVELENGTH; ++i) {
    result[i] = lhs[i] / rhs[i];
  }
  return result;
}

template<typename T>
WavelengthFunction<T> min(WavelengthFunction<T> a, WavelengthFunction<T> b) {
  WavelengthFunction<T> result;
  for (unsigned int i = 0; i < result.size(); ++i) {
    result[i] = a[i] < b[i] ? a[i] : b[i];
  }
  return result;
}

template<typename T>
WavelengthFunction<T> exp(WavelengthFunction<T> wavelength_function) {
  WavelengthFunction<T> result;
  for (unsigned int i = 0; i < result.size(); ++i) {
    result[i] = exp(wavelength_function[i]);
  }
  return result;
}

// Returns the integral of the given function over its range of wavelengths,
// i.e. integral of function*dlambda from MIN_WAVELENGTH to MAX_WAVELENGTH.
template<int L, int WL, int SA, int P, int LP>
Scalar<L, WL + 1, SA, P, LP> Integral(
    WavelengthFunction<Scalar<L, WL, SA, P, LP>> wavelength_function) {
  Scalar<L, WL, SA, P, LP> sum = 0.0 * Scalar<L, WL, SA, P, LP>::Unit();
  for (unsigned int i = 0; i < wavelength_function.size(); ++i) {
    sum = sum + wavelength_function[i];
  }
  Wavelength delta_lambda =
      (spectrum::MAX_WAVELENGTH - spectrum::MIN_WAVELENGTH) /
          spectrum::NUM_WAVELENGTH;
  return sum * delta_lambda;
}

// Returns the integral of the given function over a subset of its range,
// i.e. integral of function*dlambda from min_wavelength to max_wavelength,
// using the specified number of samples.
template<int L, int WL, int SA, int P, int LP>
Scalar<L, WL + 1, SA, P, LP> Integral(
    WavelengthFunction<Scalar<L, WL, SA, P, LP>> wavelength_function,
     Wavelength min_wavelength, Wavelength max_wavelength,
     int number_of_wavelengths) {
  const Wavelength delta_lambda =
      (max_wavelength - min_wavelength) / number_of_wavelengths;
  Wavelength lambda = min_wavelength;
  Scalar<L, WL, SA, P, LP> sum = 0.0 * Scalar<L, WL, SA, P, LP>::Unit();
  for (unsigned int i = 0; i < number_of_wavelengths; ++i) {
    sum = sum + wavelength_function(lambda);
    lambda = lambda + delta_lambda;
  }
  return sum * delta_lambda;
}

// Returns the integral of the given function over a subset of its range,
// i.e. integral of function*dlambda from min_wavelength to max_wavelength.
template<int L, int WL, int SA, int P, int LP>
Scalar<L, WL + 1, SA, P, LP> Integral(
    WavelengthFunction<Scalar<L, WL, SA, P, LP>> wavelength_function,
    Wavelength min_wavelength, Wavelength max_wavelength) {
  assert(min_wavelength >= spectrum::MIN_WAVELENGTH);
  assert(max_wavelength <= spectrum::MAX_WAVELENGTH);
  Scalar<L, WL + 1, SA, P, LP> sum = 0.0 * Scalar<L, WL + 1, SA, P, LP>::Unit();
  for (unsigned int i = 0; i < wavelength_function.size(); ++i) {
    Wavelength dlambda =
        std::min(wavelength_function.GetWavelength(i + 1), max_wavelength) -
        std::max(wavelength_function.GetWavelength(i), min_wavelength);
    sum = sum + wavelength_function[i] * std::max(dlambda, 0.0 * nm);
  }
  return sum;
}

#endif  // PHYSICS_SPECTRUM_H_
