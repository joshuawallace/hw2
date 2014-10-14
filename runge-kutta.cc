#include "runge-kutta.h"
#include "model.h"
#include <stdio.h>
#include <stdlib.h>

Runge_Kutta::Runge_Kutta(double dt, const Model &model)
  : dimen_(model.dimen()), //number of equations
    dt_(dt),               //time step size
    model_(model) {
  fx_ = new double[dimen_]; //used to hold f(x,t) of current iteration
    }

Runge_Kutta::~Runge_Kutta(){
  delete [] fx_;
}

int Runge_Kutta::Step(double t, double *x){
  
  int check = -1; //used to check the output of the f(x,t) calculation
  double k[dimen_][4]; //to store the k values as we go along
  
  check = model_.rhs(t,x,fx_); //calculate f(x,t)
  if(check != 0){ //if f(x,t) failed to be calculated
    fprintf(stderr, "Looks like the code failed to calculate f(x,t).  Sorry!  Now exiting...");
    exit(1);
  }
  
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
  
  check = model_.rhs(t+dt_/2.,xtofeedin,fx_); //calculate f(x,t)
  if(check != 0){ //if f(x,t) failed to be calculated
    fprintf(stderr, "Looks like the code failed to calculate f(x,t).  Sorry!  Now exiting...");
    exit(1);
  }
  
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
  
  check = model_.rhs(t+dt_/2.,xtofeedin,fx_); //calculate f(x,t)
  if(check != 0){ //if f(x,t) failed to be calculated
    fprintf(stderr, "Looks like the code failed to calculate f(x,t).  Sorry!  Now exiting...");
    exit(1);
  }
  
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
  
  check = model_.rhs(t+dt_,xtofeedin,fx_); //calculate f(x,t)
  if(check != 0){ //if f(x,t) failed to be calculated
    fprintf(stderr, "Looks like the code failed to calculate f(x,t).  Sorry!  Now exiting...");
    exit(1);
  }
  
  //calculates k3, based on fx above, which was based on xtofeedin above
  for(int i=0;i<dimen_;i++)
    {
      k[i][3]=fx_[i];
    }
  
  /////////////////Now to actually update the x values
  
  for(int i=0;i<dimen_;i++)
    {
      x[i] = x[i] + dt_/6.*k[i][0] + dt_/3.*(k[i][1] + k[i][2]) + dt_/6.*k[i][3]; //Uses the Runge-Kutta algorithm to step x
    }

  return 0;
}
