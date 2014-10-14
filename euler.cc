#include "euler.h"
#include "model.h"
#include <stdio.h>
#include <stdlib.h>

Euler::Euler(double dt, const Model &model)
  : dimen_(model.dimen()), //number of equations
    dt_(dt),             // time step size
    model_(model) {
  fx_ = new double[dimen_]; //used to hold f(x,t)
    }

Euler::~Euler() {
  delete [] fx_;
}

int Euler::Step(double t, double *x) {
  int check = model_.rhs(t, x, fx_);  //calculate f(x,t)
  if(check != 0){ //if f(x,t) failed to be calculated
    fprintf(stderr, "Looks like the code failed to calculate f(x,t).  Sorry!  Now exiting...");
    exit(1);
  }
  
  for (int i = 0; i < dimen_; ++i) {
    x[i] += dt_ * fx_[i];  //euler algorithm stepping method
  }
  
  return 0;
}
