#ifndef ADAMSBASHFORTH_H_
#define ADAMSBASHFORTH_H_

#include "integrator.h"
class Model;

class Adams_Bashforth : public Integrator{
 public:
  Adams_Bashforth(double dt, const Model &model);
  ~Adams_Bashforth();
  int Step(double t, double *x);
 private:
  const int dimen_;   //dimension of state x, which stores x and x'
  const double dt_;   //the timestep
  const Model &model_;//functor that evaluates f(x,t)
  double *fx_;        //holds f(x,t)
  double *fx_previous_;//hold f(x,t) of previous timestep, unique to ab2 method
  
  int isfirstiteration; //lets us know if it is the first iteration.  Initial value is 1, meaning, yes, it is the first iteration.  This is set to 0 when the first iteration happens.
};


#endif   // ADAMSBASHFORTH_H_
