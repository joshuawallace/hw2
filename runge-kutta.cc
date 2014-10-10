#include "runge-kutta.h"
#include "model.h"

Runge_Kutta::Runge_Kutta(double dt, const Model &model)
  : dimen_(model.dimen()),
    dt_(dt),
    model_(model) {
  fx_ = new double[dimen_];
    }

Runge_Kutta::~Runge_Kutta(){
  delete [] fx_;
}

int Runge_Kutta::Step(double t, double *x){
  model_.rhs(t,x,fx_);
  for(int i = 0; i < dimen_; i++) {
    double k[4]; //the corresponding 4 k values used in Runge_Kutta method
    double xtofeedin; //to store the x values we will need to use to find the k's
    //first, calculate k1
    k[0]=fx[i];

    //find the x value we need to plug in to find k2
    xtofeedin=x[i] + dt_/2. * k[0];

    //now find k2, first evaluate rhs then get k2
    model_.rhs(t,x,fx
  }
  return 0;
}
