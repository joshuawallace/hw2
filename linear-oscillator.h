#ifndef LINEAR_OSCILLATOR_H_
#define LINEAR_OSCILLATOR_H_

#include "model.h"

// Forced LinearOscillator oscillator:
//   \dot x = y
//   \dot y = - x  - 2 * beta * y + gamma * cos(omega t)
class LinearOscillator : public Model {
 public:
  LinearOscillator(double beta, double gamma, double omega);
  ~LinearOscillator();
  int rhs(double t, const double *x, double *fx) const;//used to calculate fx
  int dimen() const { return kDimen; }   //number of equations
 private:
  const double beta_;                         // damping constant (I think)
  const double gamma_;                        // amplitude of forcing
  const double omega_;                        // frequency of forcing
  static const int kDimen = 2;
};

#endif  // LINEAR_OSCILLATOR_H_
