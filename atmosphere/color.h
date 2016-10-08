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
#ifndef ATMOSPHERE_COLOR_H_
#define ATMOSPHERE_COLOR_H_

#include "physics/cie.h"
#include "physics/spectrum.h"

// Converts the given spectrum to a RGB color using the "naive", and incorrect,
// method from the Nishita93, Nishita96, O'Neal and Bruneton models, i.e. by
// simply using the spectrum value at 3 wavelengths for red, green and blue (and
// assuming a wavelength independent solar spectrum). This bypasses the CIE
// color matching functions, which are needed to correctly convert spectral
// quantities to perceptual quantities.
Color GetSrgbColorNaive(const RadianceSpectrum& spectrum);

// Extracts 3 values from the given spectrum, and convert them "correctly" to an
// sRGB color, as if the whole spectrum was converted with GetSrgbColor. This
// method is based on an approximate reconstruction of the given spectrum around
// each of the 3 extracted values, using our knowledge that this spectrum is
// proportional to the solar spectrum and approximately proportional to
// lambda^(-alpha), with alpha < 4 due to aerosols. In the end, the conversion
// from the 3 samples to an sRGB color is done with a simple constant,
// precomputed 3 vector.
Color GetSrgbColorFrom3SpectrumSamples(const RadianceSpectrum& spectrum);

// Returns an approximate version of the given spectrum, reconstructed from 3 of
// its samples with the CIE S0, S1 and S2 components. The 3 samples are first
// converted to sRGB with GetSrgbColorFrom3SpectrumSamples(), then to XYZ and
// xyY, and finally to a full spectrum using the CIE S0, S1 and S2 components.
RadianceSpectrum GetApproximateSpectrumFrom3SpectrumSamples(
    const RadianceSpectrum& spectrum);

// Performs a "naive" white balance, by dividing the given color by the
// (normalized) sun color.
Color WhiteBalanceNaive(const Color& c);

#endif  // ATMOSPHERE_COLOR_H_
