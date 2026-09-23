#pragma once

#include "parameters.hpp"
#include "types.hpp"

#include <array>
#include <complex>

struct IntegrationCoefficients {
  SpectralField e1, e2, q1, q2, q3, q4, q5, f1, f2, f3;

  [[nodiscard]] auto fields() const {
    return std::array{&e1, &e2, &q1, &q2, &q3, &q4, &q5, &f1, &f2, &f3};
  }
};

struct IntegrationFinishValues {
  Complex initial;
  Complex stageA;
  Complex nonlinear1;
  Complex nonlinear2;
  Complex nonlinear3;
  Complex nonlinear4;
};

inline Complex integrationStageA(Integrator method, double h, std::size_t i,
                                 const IntegrationCoefficients &c, Complex w,
                                 Complex n1) {
  if (method == Integrator::integratingFactorRk2)
    return c.e1[i] * (w + h * n1);
  if (method == Integrator::etd2)
    return c.e1[i] * w + c.q1[i] * n1;
  return c.e2[i] * w + c.q1[i] * n1;
}

inline Complex integrationStageB(std::size_t i,
                                 const IntegrationCoefficients &c, Complex w,
                                 Complex n1, Complex n2) {
  return c.e2[i] * w + c.q2[i] * n1 + c.q3[i] * n2;
}

inline Complex integrationStageC(std::size_t i,
                                 const IntegrationCoefficients &c, Complex w,
                                 Complex n1, Complex n3) {
  return c.e1[i] * w + c.q4[i] * n1 + c.q5[i] * n3;
}

inline Complex integrationFinish(Integrator method, double h, std::size_t i,
                                 const IntegrationCoefficients &c,
                                 const IntegrationFinishValues &v) {
  if (method == Integrator::integratingFactorRk2)
    return c.e1[i] * v.initial +
           0.5 * h * (c.e1[i] * v.nonlinear1 + v.nonlinear2);
  if (method == Integrator::etd2)
    return v.stageA + c.f1[i] * (v.nonlinear2 - v.nonlinear1);
  return c.e1[i] * v.initial + c.f1[i] * v.nonlinear1 +
         2.0 * c.f2[i] * (v.nonlinear2 + v.nonlinear3) + c.f3[i] * v.nonlinear4;
}
