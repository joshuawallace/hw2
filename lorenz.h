#ifndef LORENZ_H_
#define LORENZ_H_

#include "model.h"

// Lorenz oscillator:
//   \dot x = sigma * (y - x)
//   \dot y = rho * x - y - x * z
//   \dot z = -beta * z + x * y
class Lorenz : public Model {
 public:
  Lorenz(double sigma, double rho, double beta);
  ~Lorenz();
  int rhs(double t, const double *x, double *fx) const; //used to calculate fx
  int dimen() const { return kDimen; }    //number of equations
 private:
  const double sigma_;                        // I honestly don't know
  const double rho_;                          // what each of these
  const double beta_;                         // constants mean
  static const int kDimen = 3;
};

#endif  // LORENZ_H_
