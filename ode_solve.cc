#include <stdio.h>
#include <stdlib.h>
#include <cstring>

//diff eq header files
#include "duffing.h"
#include "lorenz.h"
#include "linear-oscillator.h"

//integration scheme header files
#include "euler.h"
#include "runge-kutta.h"
#include "adams-bashforth.h"

#include "model.h"

void intromessage(char *name)
{
  //triggers if there aren't the right number of command line arguments
  printf("Usage: %s <equation> <integrator> <timestep> <num_steps>\n",name);
  exit(0);
}

void PrintState(double n, double t, const double *x) {
  // Print a line
  //    time x[0] x[1] ...
  // to standard out
  printf("%15.8f", t);
  for (int i = 0; i < n; ++i) {
    printf("%15.8f", x[i]);
  }
  printf("\n");
}

int main(int argc, char **argv){
  
  if(argc != 5) {
    intromessage(argv[0]); //triggers if the wrong number of command line arguments is entered
  }
  
  
  const double dt = atof(argv[3]); //pull timestep in from command line
  const int nsteps = atoi(argv[4]); //pull number of steps in from command line
  Model *model; //make room for the info for the ODE to solve
  double *x; //make room for x, x', ...etc.

  if(strncmp(argv[1],  "duffing",10)==0 ){ //if duffing is chosen
    const double delta = 0.2; //parameter for equation
    const double gamma = 0.3; //parameter for equation
    const double omega = 1.0; //parameter for equation
    model = new Duffing(delta, gamma, omega); //Set up stuff to run duffing
    x = new double[model->dimen()]; //make room in memory for x
    x[0]=0.0; //initial condition
    x[1]=0.0; //initial condition
  }
  else if(strncmp(argv[1],  "lorenz",10)==0 ){ //if lorenz is chosen
    const double sigma = 10.; //parameter for equation
    const double rho = 28.; //parameter for equation
    const double beta = 8./3.; //parameter for equation
    model = new Lorenz(sigma, rho, beta); //Set up stuff to run lorenz
    x = new double[model->dimen()]; //make room in memory for x
    x[0]=0.0; //initial conditio
    x[1]=0.01; //initial conditio
    x[2]=0.0; //initial conditio
  }
  else if(strncmp(argv[1], "linear",10)==0 ){ //if linear is chosen
    const double beta = 0.1; //parameter for equation
    const double gamma = 1.; //parameter for equation
    const double omega = 0.9; //parameter for equation
    model = new LinearOscillator(beta, gamma, omega);//Set up stuff to run linear oscillator
    x = new double[model->dimen()]; //make room in memory for x
    x[0]=0.0; //initial conditio
    x[1]=0.0; //initial conditio
  }

  else{ //if the entered equation name wasn't recognized, let user know and quit
    printf("Didn't recognize the EQUATION name you entered, which was %s\nNow exiting...\n",argv[1]);
    return -1;
  }


  Integrator *integrator; //make room for the integrator

  if(strncmp(argv[2], "euler",10)==0 ){ //if euler is selected
    integrator = new Euler(dt, *model); //use the Euler algorithm
  }

  else if(strncmp(argv[2], "rk4",10)==0 ){ //if rk4 is selected
    integrator = new Runge_Kutta(dt, *model); //use the Runge-Kutta algorithm
  }

  else if(strncmp(argv[2], "ab2",10)==0 ){ //if ab2 is selected
    integrator = new Adams_Bashforth(dt, *model); //use the Adams-Bashforth algorithm
  }

  else{ //if entered integrator wasn't recognized, let user know and quit
    printf("Didn't recognize the INTEGRATOR name you entered, which was %s\nNow exiting...\n",argv[2]);
    delete model;
    delete [] x;
    return -1;
  }

  double t=0.; //The time
  PrintState(model->dimen(), t, x); //Print initial conditions
  for(int i = 0; i < nsteps; ++i){
    integrator->Step(t, x); //step the integration
    t = (i+1) * dt; //update t
    PrintState(model->dimen(), t, x); //print out the results of the step
  }

  delete integrator;
  delete model;
  delete [] x;
  return 0;
}
