#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_

#include "integrator.h"
class Model;

class Runge_Kutta : public Integrator{
 public:
  Runge_Kutta(double dt, const Model &model);
  ~Runge_Kutta();
  int Step(double t, double *x); //steps the integration forward
 private:
  const int dimen_;   //dimension of state x, which stores x and x'
  const double dt_;   //the timestep
  const Model &model_;//functor that evaluates f(x,t)
  double *fx_;        //holds f(x,t)
};

#endif   // RUNGEKUTTA_H_
