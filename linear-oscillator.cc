#include "linear-oscillator.h"
#include <math.h>

LinearOscillator::LinearOscillator(double beta, double gamma, double omega)
    : beta_(beta),
      gamma_(gamma),
      omega_(omega)
{}

LinearOscillator::~LinearOscillator()
{}

int LinearOscillator::rhs(double t, const double *x, double *fx) const {
  fx[0] = x[1];
  fx[1] = -x[0] - 2 * beta_ * x[1] + gamma_ * cos(omega_*t);
  return 0;
}
