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
#include "atmosphere/model/bruneton/core.h"

#include <thread>

#include "atmosphere/atmosphere.h"

using dimensional::vec2;
using dimensional::vec3;
using std::max;
using std::min;

#define IN(x) const x&
#define OUT(x) x&

constexpr Length Rg = EarthRadius;
constexpr Length Rt = AtmosphereRadius;
constexpr Length RL = Rt + 1.0 * km;

constexpr Length HR = RayleighScaleHeight;
constexpr Length HM = MieScaleHeight;

void SetLayer(int layer, Length* r, Length* dmin, Length* dmax, Length* dminp,
    Length* dmaxp) {
  double u = (layer * layer) / ((RES_R - 1.0) * (RES_R - 1.0));
  *r = sqrt(Rg * Rg + u * (Rt * Rt - Rg * Rg)) +
      (layer == 0 ? 0.01 : (layer == RES_R - 1 ? -0.001 : 0.0)) * km;
  *dmin = Rt - *r;
  *dmax = sqrt((*r) * (*r) - Rg * Rg) + sqrt(Rt * Rt - Rg * Rg);
  *dminp = *r - Rg;
  *dmaxp = sqrt((*r) * (*r) - Rg * Rg);
}

void ComputeTransmittance(TransmittanceTexture* transmittance_sampler) {
  std::cout << "Computing transmittance..." << std::endl;
  for (unsigned int j = 0; j < TRANSMITTANCE_H; ++j) {
    for (unsigned int i = 0; i < TRANSMITTANCE_W; ++i) {
      transmittance_sampler->Set(
          i, j, ComputeTransmittance(vec2(i + 0.5, j + 0.5)));
    }
  }
  std::cout << "Done." << std::endl;
}

void ComputeSkyIrradiance1(const TransmittanceTexture& transmittance_sampler,
    SkyIrradianceTexture* sky_irradiance) {
  std::cout << "Computing sky irradiance1..." << std::endl;
  for (unsigned int j = 0; j < IRRADIANCE_H; ++j) {
    for (unsigned int i = 0; i < IRRADIANCE_W; ++i) {
      sky_irradiance->Set(i, j, ComputeSkyIrradiance1(
          transmittance_sampler, vec2(i + 0.5, j + 0.5)));
    }
  }
  std::cout << "Done." << std::endl;
}

void ComputeInscatter1(const TransmittanceTexture& transmittance_sampler,
    IrradianceTexture* rayleigh_single_scatter_sampler,
    IrradianceTexture* mie_single_scatter_sampler) {
  std::cout << "Computing inscatter1..." << std::endl;
  for (unsigned int k = 0; k < RES_R; ++k) {
    std::cout << k << " / " << RES_R << std::endl;
    Length r, dmin, dmax, dminp, dmaxp;
    SetLayer(k, &r, &dmin, &dmax, &dminp, &dmaxp);
    for (unsigned int j = 0; j < RES_MU; ++j) {
      for (unsigned int i = 0; i < RES_MU_S * RES_NU; ++i) {
        IrradianceSpectrum ray;
        IrradianceSpectrum mie;
        ComputeInscatter1(transmittance_sampler, r, dmin, dmax, dminp, dmaxp,
            vec2(i + 0.5, j + 0.5), ray, mie);
        rayleigh_single_scatter_sampler->Set(i, j, k, ray);
        mie_single_scatter_sampler->Set(i, j, k, mie);
      }
    }
  }
  std::cout << "Done." << std::endl;
}

namespace {

constexpr unsigned int NUM_THREADS = 8;

struct ComputeInscatterSJob {
  const TransmittanceTexture* transmittance_sampler;
  const SkyIrradianceTexture* sky_irradiance_texture;
  const IrradianceTexture* rayleigh_sampler;
  const IrradianceTexture* mie_sampler;
  const RadianceTexture* raymie_sampler;
  bool first_iteration;
  RadianceDensityTexture* raymie;

  void run(unsigned int thread_id) const {
    for (unsigned int k = 0; k < RES_R; ++k) {
      if (k % NUM_THREADS != thread_id) {
        continue;
      }
      Length r, dmin, dmax, dminp, dmaxp;
      SetLayer(k, &r, &dmin, &dmax, &dminp, &dmaxp);
      for (unsigned int j = 0; j < RES_MU; ++j) {
        std::cout << thread_id << ": " << k << " / " << RES_R << " "
            << j << " / " << RES_MU << std::endl;
        for (unsigned int i = 0; i < RES_MU_S * RES_NU; ++i) {
          RadianceDensitySpectrum inscatterS;
          ComputeInscatterS(*transmittance_sampler, *sky_irradiance_texture,
              *rayleigh_sampler, *mie_sampler, *raymie_sampler,
              r, dmin, dmax, dminp, dmaxp, vec2(i + 0.5, j + 0.5),
              first_iteration, inscatterS);
          raymie->Set(i, j, k, inscatterS);
        }
      }
    }
  }
};

}  // anonymous namespace

void ComputeInscatterS(const TransmittanceTexture& transmittance_sampler,
    const SkyIrradianceTexture& sky_irradiance_texture,
    const IrradianceTexture& rayleigh_sampler,
    const IrradianceTexture& mie_sampler, const RadianceTexture& raymie_sampler,
    bool first_iteration, RadianceDensityTexture* raymie) {
  std::cout << "Computing inscatterS..." << std::endl;
  ComputeInscatterSJob job;
  job.transmittance_sampler = &transmittance_sampler;
  job.sky_irradiance_texture = &sky_irradiance_texture;
  job.rayleigh_sampler = &rayleigh_sampler;
  job.mie_sampler = &mie_sampler;
  job.raymie_sampler = &raymie_sampler;
  job.first_iteration = first_iteration;
  job.raymie = raymie;
  std::thread threads[NUM_THREADS];
  for (unsigned int thread_id = 0; thread_id < NUM_THREADS; ++thread_id) {
    threads[thread_id] =
        std::thread([job, thread_id]() { job.run(thread_id); });
  }
  for (unsigned int thread_id = 0; thread_id < NUM_THREADS; ++thread_id) {
    threads[thread_id].join();
  }
  std::cout << "Done." << std::endl;
}

void ComputeSkyIrradianceN(const IrradianceTexture& rayleigh_sampler,
    const IrradianceTexture& mie_sampler, const RadianceTexture& raymie_sampler,
    bool first_iteration, SkyIrradianceTexture* sky_irradiance) {
  std::cout << "Computing sky irradianceN..." << std::endl;
  for (unsigned int j = 0; j < IRRADIANCE_H; ++j) {
    for (unsigned int i = 0; i < IRRADIANCE_W; ++i) {
      sky_irradiance->Set(i, j, ComputeSkyIrradianceN(rayleigh_sampler,
          mie_sampler, raymie_sampler, vec2(i + 0.5, j + 0.5),
          first_iteration));
    }
  }
  std::cout << "Done." << std::endl;
}

void ComputeInscatterN(const TransmittanceTexture& transmittance_sampler,
    const RadianceDensityTexture& radiance_density_sampler,
    RadianceTexture* raymie) {
  std::cout << "Computing inscatterN..." << std::endl;
  for (unsigned int k = 0; k < RES_R; ++k) {
    std::cout << k << " / " << RES_R << std::endl;
    Length r, dmin, dmax, dminp, dmaxp;
    SetLayer(k, &r, &dmin, &dmax, &dminp, &dmaxp);
    for (unsigned int j = 0; j < RES_MU; ++j) {
      for (unsigned int i = 0; i < RES_MU_S * RES_NU; ++i) {
        RadianceSpectrum inscatterN;
        ComputeInscatterN(transmittance_sampler, radiance_density_sampler,
            r, dmin, dmax, dminp, dmaxp, vec2(i + 0.5, j + 0.5), inscatterN);
        raymie->Set(i, j, k, inscatterN);
      }
    }
  }
  std::cout << "Done." << std::endl;
}

#include "atmosphere/model/bruneton/core.inl"
