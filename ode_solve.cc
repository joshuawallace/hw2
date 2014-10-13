#include <stdio.h>
#include <stdlib.h>

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
  printf("Usage: %s <equation> <integrator> <timestep> <num_steps>\n",name);
  exit(0);
}

// Print a line
//    time x[0] x[1] ...
// to standard out
void PrintState(double n, double t, const double *x) {
  printf("%15.8f", t);
  for (int i = 0; i < n; ++i) {
    printf("%15.8f", x[i]);
  }
  printf("\n");
}

int main(int argc, char **argv){

  //  const char *equation_name[120];
  //  const char *integrator_name[120];


  if(argc != 5) {
    intromessage(argv[0]);
  }

  const double dt = atof(argv[3]);
  const int nsteps = atoi(argv[4]);

  if(argv[1] == "duffing"){
    const double delta = 0.2;
    const double gamma = 0.3;
    const double omega = 1.0;
    Model *model = new Duffing(delta, gamma, omega);
  }

  else if(argv[1] == "lorenz"){
    const double sigma = 10.;
    const double rho = 28.;
    const double beta = 8./3.;
    Model *model = new Lorenz(sigma, rho, beta);
  }

  else if(argv[1] == "linear-oscillator"){
    const double beta = 0.1;
    const double gamma = 1.;
    const double omega = 0.9;
    Model *model = new LinearOscillator(beta, gamma, omega);
  }

  else{
    printf("Didn't recognize the EQUATION name you entered, which was %s\nNow exiting...\n",argv[1]);
    return -1;
  }

  if(argv[2] == "euler"){
    Integrator *integrator = new Euler(dt, *model);
  }

  else if(argv[2] == "runge-kutta"){
    Integrator *integrator = new Runge_Kutta(dt, *model);
  }

  else if(argv[2] == "euler"){
    Integrator *integrator = new Adams_Bashforth(dt, *model);
  }

  else{
    printf("Didn't recognize the INTEGRATOR name you entered, which was %s\nNow exiting...\n",argv[2]);
    delete model;
    return -1;
  }

  double t=0.;
  PrintState(dimen, t, x);
  for(int i = 0; i < nsteps; ++i){
    integrator->Step(t, x);
    t = (i+1) * dt;
    PrintState(dimen, t, x);
  }

  delete integrator;
  delete model;
  return 0;
}
