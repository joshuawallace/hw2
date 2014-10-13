#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <cmath>

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
  printf("Usage: %s <integrator> <timestep> <num_steps>\n",name);
  exit(0);
}

double x_of_linear_osc(const double t, const double beta, const double omega, const double gamma, const double c1, const double c2, const double A, const double B){
  //const double A = ( (1 - omega*omega) * gamma) / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);
  //const double B = 2 * beta * omega * gamma / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);
  const double omega_d = sqrt(1 - beta*beta);

  return exp(-beta*t) * (c1 * cos(omega_d * t) + c2 * sin(omega_d * t)) + A * cos(omega*t) + B * sin(omega*t);
}

double x_prime_of_linear_osc(const double t, const double beta, const double omega, const double gamma, const double c1, const double c2, const double A, const double B){
  // const double A = ( (1 - omega*omega) * gamma) / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);
  //const double B = 2 * beta * omega * gamma / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);
  const double omega_d = sqrt(1 - beta*beta);

  return exp(-beta*t) * ( -beta * (c1 * cos(omega_d * t) + c2 * sin(omega_d * t)) + omega_d * ( -c1 * sin(omega_d * t) + c2 * cos(omega_d*t)) ) - omega * ( A * sin(omega*t) - B * cos(omega*t));
}

int main(int argc, char **argv){
  
  if(argc != 4) {
    intromessage(argv[0]);
  }

  const double dt = atof(argv[2]);
  const int nsteps = atoi(argv[3]);

    const double beta = 0.1;
    const double gamma = 1.;
    const double omega = 0.9;
    Model *model = new LinearOscillator(beta, gamma, omega);

    double x[model->dimen()];
    x[0]=0.0;
    x[1]=0.0;

  Integrator *integrator;

  if(strncmp(argv[1], "euler",10)==0 ){
    integrator = new Euler(dt, *model);
  }

  else if(strncmp(argv[1], "rk4",10)==0 ){
    integrator = new Runge_Kutta(dt, *model);
  }

  else if(strncmp(argv[1], "ab2",10)==0 ){
    integrator = new Adams_Bashforth(dt, *model);
  }

  else{
    printf("Didn't recognize the INTEGRATOR name you entered, which was %s\nNow exiting...\n",argv[1]);
    delete model;
    return -1;
  }


  const double A = ( (1 - omega*omega) * gamma) / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);
  const double B = 2 * beta * omega * gamma / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);

  const double c1 = x[0] - A;
  const double c2 = (x[1] + beta * (x[0] - A) - omega * B)/sqrt(1 - beta * beta);

  double error_norm = 0.0;

  double t=0.;

  error_norm += dt * pow( (x[0] - x_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2) * pow( (x[1] - x_prime_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2);  //For the zeroeth step of the iteration, starting to calcualte the error_norm.
  //this likely will add 0 to the error_norm, since the itegration hasn't even started yet and so
  //there has been no divergence between itegration and actual.

  for(int i = 0; i < nsteps; ++i){
    integrator->Step(t, x);
    t = (i+1) * dt;
    error_norm += dt * pow( (x[0] - x_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2) * pow( (x[1] - x_prime_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2);
  }
  
  printf("%lf\n",sqrt(error_norm));
  
  delete integrator;
  delete model;
  return 0;
}
