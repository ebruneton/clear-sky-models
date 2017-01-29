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
#ifndef ATMOSPHERE_MODEL_BRUNETON_CORE_H_
#define ATMOSPHERE_MODEL_BRUNETON_CORE_H_

#include <algorithm>

#include "atmosphere/atmosphere.h"
#include "math/binary_function.h"
#include "math/scalar.h"
#include "math/ternary_function.h"
#include "math/vector.h"
#include "physics/units.h"

Length Limit(Length r, Number mu);

// TRANSMITTANCE --------------------------------------------------------------

// Texture size for the transmittance texture T.
const unsigned int TRANSMITTANCE_W = 256;
const unsigned int TRANSMITTANCE_H = 64;

typedef dimensional::BinaryFunction<
    TRANSMITTANCE_W, TRANSMITTANCE_H, DimensionlessSpectrum>
        TransmittanceTexture;

DimensionlessSpectrum GetTransmittance(
    const TransmittanceTexture& transmittance_sampler, Length r, Number mu);

DimensionlessSpectrum GetTransmittance(
    const TransmittanceTexture& transmittance_sampler, Length r, Number mu,
    Length d);

// SKY IRRADIANCE --------------------------------------------------------------

// Texture size for the sky irradiance textures E and DeltaE.
const unsigned int IRRADIANCE_W = 64;
const unsigned int IRRADIANCE_H = 16;

typedef dimensional::BinaryFunction<
    IRRADIANCE_W, IRRADIANCE_H, IrradianceSpectrum>
        SkyIrradianceTexture;

IrradianceSpectrum GetSkyIrradiance(
    const SkyIrradianceTexture& sky_irradiance_sampler, Length r, Number muS);

// SKY RADIANCE ----------------------------------------------------------------

// Texture size for the various 4D textures (S, DeltaS, DeltaJ, etc).
const unsigned int RES_R = 32;
const unsigned int RES_MU = 128;
const unsigned int RES_MU_S = 32;
const unsigned int RES_NU = 8;

typedef dimensional::Scalar<-3, -1, -1, 1, 0> SpectralRadianceDensity;

typedef WavelengthFunction<-3, -1, -1, 1, 0> RadianceDensitySpectrum;

typedef dimensional::TernaryFunction<
    RES_MU_S * RES_NU, RES_MU, RES_R, IrradianceSpectrum>
        IrradianceTexture;

typedef dimensional::TernaryFunction<
    RES_MU_S * RES_NU, RES_MU, RES_R, RadianceSpectrum>
        RadianceTexture;

typedef dimensional::TernaryFunction<
    RES_MU_S * RES_NU, RES_MU, RES_R, RadianceDensitySpectrum>
        RadianceDensityTexture;

// PRECOMPUTATIONS -------------------------------------------------------------

DimensionlessSpectrum ComputeTransmittance(dimensional::vec2 xy);

void ComputeTransmittance(TransmittanceTexture* transmittance_sampler);

void ComputeSkyIrradiance1(const TransmittanceTexture& transmittance_sampler,
    SkyIrradianceTexture* sky_irradiance);

void ComputeInscatter1(const TransmittanceTexture& transmittance_sampler,
    IrradianceTexture* rayleigh_single_scatter_sampler,
    IrradianceTexture* mie_single_scatter_sampler);

void ComputeInscatterS(const TransmittanceTexture& transmittance_sampler,
    const SkyIrradianceTexture& sky_irradiance_texture,
    const IrradianceTexture& rayleigh_sampler,
    const IrradianceTexture& mie_sampler, const RadianceTexture& raymie_sampler,
    bool first_iteration, RadianceDensityTexture* raymie);

void ComputeSkyIrradianceN(const IrradianceTexture& rayleigh_sampler,
    const IrradianceTexture& mie_sampler, const RadianceTexture& raymie_sampler,
    bool first_iteration, SkyIrradianceTexture* sky_irradiance);

void ComputeInscatterN(const TransmittanceTexture& transmittance_sampler,
    const RadianceDensityTexture& radiance_density_sampler,
    RadianceTexture* raymie);

// TEXTURE FETCH ---------------------------------------------------------------

template<unsigned int NX, unsigned int NY, typename T>
T texture2d(const dimensional::BinaryFunction<NX, NY, T>& table,
    const dimensional::vec2& uv) {
  return table(uv.x(), uv.y());
}

template<unsigned int NX, unsigned int NY, unsigned int NZ, typename T>
T texture4d(const dimensional::TernaryFunction<NX, NY, NZ, T>& table,
    Length r, Number mu, Number muS, Number nu) {
  const Length H =
      sqrt(AtmosphereRadius * AtmosphereRadius - EarthRadius * EarthRadius);
  Length rho = sqrt(r * r - EarthRadius * EarthRadius);
  Length rmu = r * mu;
  Area delta = rmu * rmu - r * r + EarthRadius * EarthRadius;
  Number cstx;
  Area csty;
  Length cstz;
  Number cstw;
  if (rmu < 0.0 * m && delta > 0.0 * m2) {
    cstx = 1.0;
    csty = 0.0 * m2;
    cstz = 0.0 * m;
    cstw = 0.5 - 0.5 / RES_MU;
  } else {
    cstx = -1.0;
    csty = H * H;
    cstz = H;
    cstw = 0.5 + 0.5 / RES_MU;
  }
  Number uR =
      0.5 / RES_R + rho / H * (1.0 - 1.0 / RES_R);
  Number uMu = cstw + (rmu * cstx + sqrt(delta + csty)) /
      (rho + cstz) * (0.5 - 1.0 / RES_MU);
  // paper formula
  // Number uMuS = 0.5 / RES_MU_S +
  //    max((1.0 - exp(-3.0 * muS - 0.6)) / (1.0 - exp(-3.6)), 0.0) *
  //        (1.0 - 1.0 / RES_MU_S);
  // better formula
  constexpr Angle a = 1.1 * rad;
  Number uMuS = 0.5 / RES_MU_S +
      (atan(max(muS, -0.1975) * tan(1.26 * a)) / a + (1.0 - 0.26)) *
          0.5 * (1.0 - 1.0 / RES_MU_S);
  Number lerp = (nu + 1.0) / 2.0 * (RES_NU - 1);
  Number uNu = floor(lerp);
  lerp = lerp - uNu;
  return table((uNu + uMuS)() / RES_NU, uMu(), uR()) * (1.0 - lerp) +
         table((uNu + uMuS + 1.0)() / RES_NU, uMu(), uR()) * lerp;
}

#endif  // ATMOSPHERE_MODEL_BRUNETON_CORE_H_
