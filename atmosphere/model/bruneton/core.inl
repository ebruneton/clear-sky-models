/**
 * Precomputed Atmospheric Scattering
 * Copyright (c) 2008 INRIA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

// Nearest intersection of ray r,mu with ground or top atmosphere boundary
// mu=cos(ray zenith angle at ray origin).
Length Limit(Length r, Number mu) {
  Length d_out = -r * mu + sqrt(r * r * (mu * mu - 1.0) + RL * RL);
  Area delta_sq = r * r * (mu * mu - 1.0) + Rg * Rg;
  if (delta_sq >= 0.0 * m2) {
    Length d_in = -r * mu - sqrt(delta_sq);
    if (d_in >= 0.0 * m) {
      d_out = min(d_out, d_in);
    }
  }
  return d_out;
}

// TRANSMITTANCE ---------------------------------------------------------------

void GetTransmittanceRMu(vec2 xy, Length* r, Number* muS) {
  Number x = xy.x / TRANSMITTANCE_W;
  Number y = xy.y / TRANSMITTANCE_H;
  *r = Rg + (y * y) * (Rt - Rg);
  *muS = -0.15 + tan(1.5 * rad * x) / tan(1.5 * rad) * (1.0 + 0.15);
}

vec2 GetTransmittanceUV(Length r, Number mu) {
	Number uR = sqrt((r - Rg) / (Rt - Rg));
	Number uMu = atan((mu + 0.15) / (1.0 + 0.15) * tan(1.5 * rad)) / (1.5 * rad);
  return vec2(uMu, uR);
}

DimensionlessSpectrum GetTransmittance(
    const TransmittanceTexture& transmittance_sampler, Length r, Number mu) {
	vec2 uv = GetTransmittanceUV(r, mu);
  return texture2d(transmittance_sampler, uv);
}

DimensionlessSpectrum GetTransmittance(
    const TransmittanceTexture& transmittance_sampler, Length r, Number mu,
    Length d) {
  static const DimensionlessSpectrum kFullTransmittance(1.0);
  Length r1 = sqrt(r * r + d * d + 2.0 * r * mu * d);
  Number mu1 = (r * mu + d) / r1;
  if (mu > 0.0) {
    return min(GetTransmittance(transmittance_sampler, r, mu) /
        GetTransmittance(transmittance_sampler, r1, mu1), kFullTransmittance);
  } else {
    return min(GetTransmittance(transmittance_sampler, r1, -mu1) /
        GetTransmittance(transmittance_sampler, r, -mu), kFullTransmittance);
  }
}

// SKY IRRADIANCE --------------------------------------------------------------

void GetSkyIrradianceRMuS(vec2 xy, Length* r, Number* muS) {
  Number x = (xy.x - 0.5) / (IRRADIANCE_W - 1);
  Number y = (xy.y - 0.5) / (IRRADIANCE_H - 1);
  *r = Rg + y * (Rt - Rg);
  *muS = -0.2 + x * (1.0 + 0.2);
}

vec2 GetSkyIrradianceUV(Length r, Number muS) {
  Number uR = (r - Rg) / (Rt - Rg);
  Number uMuS = (muS + 0.2) / (1.0 + 0.2);
  return vec2(uMuS, uR);
}

IrradianceSpectrum GetSkyIrradiance(
    const SkyIrradianceTexture& sky_irradiance_sampler, Length r, Number muS) {
  vec2 uv = GetSkyIrradianceUV(r, muS);
  return texture2d(sky_irradiance_sampler, uv);
}

// SKY RADIANCE ----------------------------------------------------------------

void GetMuMuSNu(vec2 xy, Length r,
    Length dmin, Length dmax, Length dminp, Length dmaxp,
    Number* mu, Number* muS, Number* nu) {
  Number x = xy.x - 0.5;
  Number y = xy.y - 0.5;
  if (y < RES_MU / 2.0) {
    Length d = (1.0 - y / (RES_MU / 2.0 - 1.0)) * dmaxp;
    d = min(max(dminp, d), dmaxp * 0.999);
    *mu = (Rg * Rg - r * r - d * d) / (2.0 * r * d);
    *mu = min(*mu, -sqrt(1.0 - (Rg / r) * (Rg / r)) - 0.001);
  } else {
    Length d = ((y - RES_MU / 2.0) / (RES_MU / 2.0 - 1.0)) * dmax;
    d = min(max(dmin, d), dmax * 0.999);
    *mu = (Rt * Rt - r * r - d * d) / (2.0 * r * d);
  }
  *muS = mod(x, RES_MU_S) / (RES_MU_S - 1.0);
  // paper formula
  // muS = -(0.6 + log(1.0 - muS * (1.0 -  exp(-3.6)))) / 3.0;
  // better formula
  *muS = tan((2.0 * *muS - 1.0 + 0.26) * 1.1 * rad) / tan(1.26 * 1.1 * rad);
  *nu = -1.0 + floor(x / RES_MU_S) / (RES_NU - 1) * 2.0;
}

// PRECOMPUTATIONS -------------------------------------------------------------

const int TRANSMITTANCE_INTEGRAL_SAMPLES = 500;
const int IRRADIANCE_INTEGRAL_SAMPLES = 32;
const int INSCATTER_INTEGRAL_SAMPLES = 50;
const int INSCATTER_SPHERICAL_INTEGRAL_SAMPLES = 16;

Length ComputeOpticalLength(Length scale_height, Length r, Number mu) {
  Length result = 0.0 * m;
  Length dx = Limit(r, mu) / TRANSMITTANCE_INTEGRAL_SAMPLES;
  Number yi = exp(-(r - Rg) / scale_height);
  for (int i = 1; i <= TRANSMITTANCE_INTEGRAL_SAMPLES; ++i) {
    Length xj = i * dx;
    Length rj = sqrt(r * r + xj * xj + 2.0 * xj * r * mu);
    Number yj = exp(-(rj - Rg) / scale_height);
    result += (yi + yj) / 2.0 * dx;
    yi = yj;
  }
  Number mu_horizon = -sqrt(1.0 - (Rg / r) * (Rg / r));
  return mu < mu_horizon ? 1e9 * m : result;
}

DimensionlessSpectrum ComputeTransmittance(vec2 xy) {
  Length r;
  Number muS;
  GetTransmittanceRMu(xy, &r, &muS);
  return exp(-(
      RayleighScattering() * ComputeOpticalLength(HR, r, muS) +
      MieExtinction() * ComputeOpticalLength(HM, r, muS)));
}

IrradianceSpectrum ComputeSkyIrradiance1(
    const TransmittanceTexture& transmittance_sampler, vec2 xy) {
  Length r;
  Number muS;
  GetSkyIrradianceRMuS(xy, &r, &muS);
  return GetTransmittance(transmittance_sampler, r, muS) * max(muS, 0.0) *
      SolarSpectrum();
}

void ComputeInscatter1Integrand(
    const TransmittanceTexture& transmittance_sampler,
    Length r, Number mu, Number muS, Number nu, Length t,
    DimensionlessSpectrum* rayleigh, DimensionlessSpectrum* mie) {
  *rayleigh = DimensionlessSpectrum(0.0);
  *mie = DimensionlessSpectrum(0.0);
  Length ri = sqrt(r * r + t * t + 2.0 * r * mu * t);
  Number muSi = (nu * t + muS * r) / ri;
  ri = max(Rg, ri);
  if (muSi >= -sqrt(1.0 - Rg * Rg / (ri * ri))) {
    DimensionlessSpectrum ti =
        GetTransmittance(transmittance_sampler, r, mu, t) *
        GetTransmittance(transmittance_sampler, ri, muSi);
    *rayleigh = ti * exp(-(ri - Rg) / HR);
    *mie = ti * exp(-(ri - Rg) / HM);
  }
}

void ComputeInscatter1(
    const TransmittanceTexture& transmittance_sampler, Length r, Length dmin,
    Length dmax, Length dminp, Length dmaxp, vec2 xy,
    IrradianceSpectrum* rayleigh, IrradianceSpectrum* mie) {
  Number mu, muS, nu;
  GetMuMuSNu(xy, r, dmin, dmax, dminp, dmaxp, &mu, &muS, &nu);

  DimensionlessSpectrum ray_sum = DimensionlessSpectrum(0.0);
  DimensionlessSpectrum mie_sum = DimensionlessSpectrum(0.0);
  Length dx = Limit(r, mu) / INSCATTER_INTEGRAL_SAMPLES;
  DimensionlessSpectrum rayi;
  DimensionlessSpectrum miei;
  ComputeInscatter1Integrand(
      transmittance_sampler, r, mu, muS, nu, 0.0 * m, &rayi, &miei);
  for (int i = 1; i <= INSCATTER_INTEGRAL_SAMPLES; ++i) {
    Length xj = i * dx;
    DimensionlessSpectrum rayj;
    DimensionlessSpectrum miej;
    ComputeInscatter1Integrand(
        transmittance_sampler, r, mu, muS, nu, xj, &rayj, &miej);
    ray_sum += (rayi + rayj);
    mie_sum += (miei + miej);
    rayi = rayj;
    miei = miej;
  }
  *rayleigh = ray_sum * (dx / 2) * RayleighScattering() * SolarSpectrum();
  *mie = mie_sum * (dx / 2) * MieScattering() * SolarSpectrum();
}

void ComputeInscatterS(
    const TransmittanceTexture& transmittance_sampler,
    const SkyIrradianceTexture& sky_irradiance_sampler,
    const IrradianceTexture& rayleigh_sampler,
    const IrradianceTexture& mie_sampler,
    const RadianceTexture& raymie_sampler,
    Length r, Length dmin, Length dmax, Length dminp, Length dmaxp, vec2 xy,
    bool first_iteration, RadianceDensitySpectrum* raymie) {
  Number mu, muS, nu;
  GetMuMuSNu(xy, r, dmin, dmax, dminp, dmaxp, &mu, &muS, &nu);

  const Angle dphi = pi / INSCATTER_SPHERICAL_INTEGRAL_SAMPLES;
  const Angle dtheta = pi / INSCATTER_SPHERICAL_INTEGRAL_SAMPLES;
  r = max(Rg, min(r, Rt));
  mu = max(-1.0, min(mu, 1.0));
  muS = max(-1.0, min(muS, 1.0));
  Number var = sqrt(1.0 - mu * mu) * sqrt(1.0 - muS * muS);
  nu = max(muS * mu - var, min(nu, muS * mu + var));

  Number cthetamin = -sqrt(1.0 - (Rg / r) * (Rg / r));

  vec3 v = vec3(sqrt(1.0 - mu * mu), 0.0, mu);
  Number sx = v.x == 0.0 ? 0.0 : (nu - muS * mu) / v.x;
  vec3 s = vec3(sx, sqrt(max(0.0, 1.0 - sx * sx - muS * muS)), muS);

  *raymie = RadianceDensitySpectrum(0.0 * SpectralRadianceDensity::Unit());

  // integral over 4.PI around x with two nested loops over w directions (theta,phi) -- Eq (7)
  for (int itheta = 0; itheta < INSCATTER_SPHERICAL_INTEGRAL_SAMPLES;
      ++itheta) {
    Angle theta = (itheta + 0.5) * dtheta;
    Number ctheta = cos(theta);

    PhaseFunctionSpectrum greflectance = PhaseFunctionSpectrum(0.0 / sr);
    Length dground = 0.0 * m;
    DimensionlessSpectrum gtransp = DimensionlessSpectrum(0.0);
    if (ctheta < cthetamin) { // if ground visible in direction w
      // compute transparency gtransp between x and ground
      greflectance = GroundAlbedo() * (1.0 / (PI * sr));
      dground = -r * ctheta - sqrt(r * r * (ctheta * ctheta - 1.0) + Rg * Rg);
      gtransp = GetTransmittance(
          transmittance_sampler, Rg, -(r * ctheta + dground) / Rg, dground);
    }

    for (int iphi = 0; iphi < 2 * INSCATTER_SPHERICAL_INTEGRAL_SAMPLES;
        ++iphi) {
      Angle phi = (iphi + 0.5) * dphi;
      SolidAngle dw = dtheta.to(rad) * dphi.to(rad) * sin(theta) * sr;
      vec3 w = vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), ctheta);

      Number nu1 = dot(s, w);
      Number nu2 = dot(v, w);
      InverseSolidAngle pr2 = RayleighPhaseFunction(nu2);
      InverseSolidAngle pm2 = MiePhaseFunction(nu2);

      // compute irradiance received at ground in direction w (if ground visible) =deltaE
      vec3 gnormal = vec3(0.0, 0.0, r / Rg) + w * (dground / Rg);
      IrradianceSpectrum girradiance =
          GetSkyIrradiance(sky_irradiance_sampler, Rg, dot(gnormal, s));

      RadianceSpectrum raymie1; // light arriving at x from direction w

      // first term = light reflected from the ground and attenuated before reaching x, =T.alpha/PI.deltaE
      raymie1 = girradiance * greflectance * gtransp;

      // second term = inscattered light, =deltaS
      if (first_iteration) {
        // first iteration is special because Rayleigh and Mie were stored separately,
        // without the phase functions factors; they must be reintroduced here
        InverseSolidAngle pr1 = RayleighPhaseFunction(nu1);
        InverseSolidAngle pm1 = MiePhaseFunction(nu1);
        IrradianceSpectrum ray1 = texture4d(rayleigh_sampler, r, w.z, muS, nu1);
        IrradianceSpectrum mie1 = texture4d(mie_sampler, r, w.z, muS, nu1);
        raymie1 += ray1 * pr1 + mie1 * pm1;
      } else {
        raymie1 += texture4d(raymie_sampler, r, w.z, muS, nu1);
      }

      // light coming from direction w and scattered in direction v
      // = light arriving at x from direction w (raymie1) *
      //     SUM(scattering coefficient * phaseFunction)
      // see Eq (7)
      *raymie += raymie1 * (
          RayleighScattering() * exp(-(r - Rg) / HR) * pr2 +
          MieScattering() * exp(-(r - Rg) / HM) * pm2) * dw;
    }
  }
  // output raymie = J[T.alpha/PI.deltaE + deltaS] (line 7 in algorithm 4.1)
}

IrradianceSpectrum ComputeSkyIrradianceN(
    const IrradianceTexture& rayleigh_sampler,
    const IrradianceTexture& mie_sampler,
    const RadianceTexture& raymie_sampler,
    vec2 xy, bool first_iteration) {
  Length r;
  Number muS;
  GetSkyIrradianceRMuS(xy, &r, &muS);

  const Angle dphi = pi / IRRADIANCE_INTEGRAL_SAMPLES;
  const Angle dtheta = pi / IRRADIANCE_INTEGRAL_SAMPLES;
  vec3 s = vec3(sqrt(1.0 - muS * muS), 0.0, muS);
  IrradianceSpectrum result = IrradianceSpectrum(0.0 * SpectralIrradiance::Unit());
  // integral over 2.PI around x with two nested loops over w directions (theta,phi) -- Eq (15)
  for (int iphi = 0; iphi < 2 * IRRADIANCE_INTEGRAL_SAMPLES; ++iphi) {
    Angle phi = (iphi + 0.5) * dphi;
    for (int itheta = 0; itheta < IRRADIANCE_INTEGRAL_SAMPLES / 2; ++itheta) {
      Angle theta = (itheta + 0.5) * dtheta;
      SolidAngle dw = dtheta.to(rad) * dphi.to(rad) * sin(theta) * sr;
      vec3 w = vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
      Number nu = dot(s, w);
      if (first_iteration) {
        // first iteration is special because Rayleigh and Mie were stored separately,
        // without the phase functions factors; they must be reintroduced here
        InverseSolidAngle pr1 = RayleighPhaseFunction(nu);
        InverseSolidAngle pm1 = MiePhaseFunction(nu);
        IrradianceSpectrum ray1 = texture4d(rayleigh_sampler, r, w.z, muS, nu);
        IrradianceSpectrum mie1 = texture4d(mie_sampler, r, w.z, muS, nu);
        result += (ray1 * pr1 + mie1 * pm1) * w.z * dw;
      } else {
        result += texture4d(raymie_sampler, r, w.z, muS, nu) * w.z * dw;
      }
    }
  }
  return result;
}

RadianceDensitySpectrum ComputeInscatterNIntegrand(
    const TransmittanceTexture& transmittance_sampler,
    const RadianceDensityTexture& radiance_density_sampler,
    Length r, Number mu, Number muS, Number nu, Length t) {
  Length ri = max(sqrt(r * r + t * t + 2.0 * r * mu * t), Rg + 10.0 * m);
  Number mui = (r * mu + t) / ri;
  Number muSi = (nu * t + muS * r) / ri;
  return texture4d(radiance_density_sampler, ri, mui, muSi, nu) *
      GetTransmittance(transmittance_sampler, r, mu, t);
}

void ComputeInscatterN(
    const TransmittanceTexture& transmittance_sampler,
    const RadianceDensityTexture& radiance_density_sampler,
    Length r, Length dmin, Length dmax, Length dminp, Length dmaxp, vec2 xy,
    RadianceSpectrum* raymie) {
  Number mu, muS, nu;
  GetMuMuSNu(xy, r, dmin, dmax, dminp, dmaxp, &mu, &muS, &nu);

  *raymie = RadianceSpectrum(0.0 * SpectralRadiance::Unit());

  Length dx = Limit(r, mu) / INSCATTER_INTEGRAL_SAMPLES;
  RadianceDensitySpectrum raymiei = ComputeInscatterNIntegrand(
      transmittance_sampler, radiance_density_sampler, r, mu, muS, nu, 0.0 * m);
  for (int i = 1; i <= INSCATTER_INTEGRAL_SAMPLES; ++i) {
    Length xj = i * dx;
    RadianceDensitySpectrum raymiej = ComputeInscatterNIntegrand(
        transmittance_sampler, radiance_density_sampler, r, mu, muS, nu, xj);
    *raymie += (raymiei + raymiej) * (dx / 2.0);
    raymiei = raymiej;
  }
}
