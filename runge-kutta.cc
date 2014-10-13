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

  double k[dimen_][4]; //to store the k values as we go along
  model_.rhs(t,x,fx_);
  
  //calculates k1, based on fx_ above
  for(int i=0;i<dimen_;i++)
    {
      k[i][0]=fx_[i];
    }
  
  //allows us to set x to feed in to the rhs function as x + h/2 * k_1
  //since k_1 is going to be different for each of the n elements of x
  double xtofeedin[dimen_];
  
  for(int i=0;i<dimen_;i++)
    {
      xtofeedin[i]=x[i]+dt_/2 * k[i][0];
    }
  
  model_.rhs(t+dt_/2.,xtofeedin,fx_);
  
  //calculates k2, based on fx above, which was based on xtofeedin above
  for(int i=0;i<dimen_;i++)
    {
      k[i][1]=fx_[i];
    }
  
  
  //Now for k_3
  //First, xtofeedin
  for(int i=0;i<dimen_;i++)
    {
      xtofeedin[i]=x[i]+dt_/2. * k[i][1];
    }
  
  model_.rhs(t+dt_/2.,xtofeedin,fx_);
  
  //calculates k3, based on fx above, which was based on xtofeedin above
  for(int i=0;i<dimen_;i++)
    {
      k[i][2]=fx_[i];
    }
  
  
  //Now for k_4
  //First, xtofeedin
  for(int i=0;i<dimen_;i++)
    {
      xtofeedin[i]=x[i]+dt_ * k[i][2];
    }
  
  model_.rhs(t+dt_,xtofeedin,fx_);
  
  //calculates k3, based on fx above, which was based on xtofeedin above
  for(int i=0;i<dimen_;i++)
    {
      k[i][3]=fx_[i];
    }
  
  /////////////////Now to actually update the x values
  
  for(int i=0;i<dimen_;i++)
    {
      x[i] = x[i] + dt_/6.*k[i][0] + dt_/3.*(k[i][1] + k[i][2]) + dt_/6.*k[i][3];
    }
  return 0;
}
