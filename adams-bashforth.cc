#include "adams-bashforth.h"
#include "model.h"

Adams_Bashforth::Adams_Bashforth(double dt, const Model &model)
  : dimen_(model.dimen()), //number of equations
    dt_(dt),               //time step size
    model_(model) {
  fx_ = new double[dimen_];  //used to hold f(x,t) of current iteration
  fx_previous_ = new double[dimen_]; //used to hold the f(x,t) values of the previous iteration, which is necessary for the Adams-Bashforth algorithm
  isfirstiteration = 1; //tracks whether it is the first iteration or not, defined here initially to be 1, which means "yes"

    }

Adams_Bashforth::~Adams_Bashforth(){
  delete [] fx_;
  delete [] fx_previous_;
}

int Adams_Bashforth::Step(double t, double *x){
  
  if(isfirstiteration == 1){ //if it is the first iteration
    //implement Runge-Kutta for first iteration
    
    double k[dimen_][4]; //to store the k values as we go along
    model_.rhs(t,x,fx_); //calculate f(x,t)
    
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
    
    model_.rhs(t+dt_/2.,xtofeedin,fx_); //calculate f(x',t')
    
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
    
    model_.rhs(t+dt_/2.,xtofeedin,fx_); //calculate f(x',t')
    
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
    
    model_.rhs(t+dt_,xtofeedin,fx_); //calculate f(x',t')
    
    //calculates k3, based on fx above, which was based on xtofeedin above
    for(int i=0;i<dimen_;i++)
      {
        k[i][3]=fx_[i];
      }
    
    /////////////////Now to actually update the x values
    
    for(int i=0;i<dimen_;i++)
      {
        x[i] = x[i] + dt_/6.*k[i][0] + dt_/3.*(k[i][1] + k[i][2]) + dt_/6.*k[i][3]; //First step uses Runge-Kutta method to step x
        fx_previous_[i]=fx_[i]; //update fx_previous_ for the next iteration
      }

    isfirstiteration = 0;//so that next time the computer will no it's no longer the first iteration
  }

  else{ //All but the first iteration
    
    model_.rhs(t,x,fx_); //find f(x,t)

    for(int i=0;i<dimen_;i++){
      x[i] = x[i] + 1.5 * dt_ * fx_[i] - 0.5 * dt_ * fx_previous_[i]; //Uses the Adams-Bashforth algorithm to step x
      fx_previous_[i]=fx_[i]; //update fx_previous_ for the next iteration
    }
  }
  return 0;
}
