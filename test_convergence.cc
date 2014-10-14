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

void intromessage(char *name, int argc)
{
  //triggers if there aren't the right number of command line arguments
  printf("Usage: %s <integrator> <timestep> <num_steps>, you have %d arguments\n",name,argc);
  exit(0);
}

double x_of_linear_osc(const double t, const double beta, const double omega, const double gamma, const double c1, const double c2, const double A, const double B){
  //calculates the exact x value from the analytic solution

  const double omega_d = sqrt(1 - beta*beta); //one of the convenient parameters used

  return exp(-beta*t) * (c1 * cos(omega_d * t) + c2 * sin(omega_d * t)) + A * cos(omega*t) + B * sin(omega*t); //return the calculated x value
}

double x_prime_of_linear_osc(const double t, const double beta, const double omega, const double gamma, const double c1, const double c2, const double A, const double B){
  //calculates the exact x' value from the analytic solution

  const double omega_d = sqrt(1 - beta*beta); //one of the convenient parameters used

  return exp(-beta*t) * ( -beta * (c1 * cos(omega_d * t) + c2 * sin(omega_d * t)) + omega_d * ( -c1 * sin(omega_d * t) + c2 * cos(omega_d*t)) ) - omega * ( A * sin(omega*t) - B * cos(omega*t)); //return the calculated x' value
}

int main(int argc, char **argv){
  
  if(argc != 4) {
    intromessage(argv[0],argc); //triggers if the wrong number of command line arguments is entered
  }

  const double dt = atof(argv[2]);  //pull timestep in from command line
  const int nsteps = atoi(argv[3]); //pull number of steps in from command line

  const double beta = 0.1;  //parameter for the ODE
  const double gamma = 1.;  //parameter for the ODE
  const double omega = 0.9; //parameter for the ODE
  Model *model = new LinearOscillator(beta, gamma, omega); //Set up stuff to run the linear oscillator
  
  double x[model->dimen()]; //make array for x and x'
  x[0]=0.0; //initial conditions
  x[1]=0.0; //initial conditions
  
  Integrator *integrator; //allocate some memory for the integrator
    
  if(strncmp(argv[1], "euler",10)==0 ){ //if euler is selected
    integrator = new Euler(dt, *model); //use the Euler method
  }
  
  else if(strncmp(argv[1], "rk4",10)==0 ){ //if rk4 is selected
    integrator = new Runge_Kutta(dt, *model); //use the Runge-Kutta algorithm
  }
  
  else if(strncmp(argv[1], "ab2",10)==0 ){ //if ab2 is selected
    integrator = new Adams_Bashforth(dt, *model); //use the Adams-Bashforth algorithm
  }
  
  else{
    //if none of the above integrators were selected, let user know and quit
    printf("Didn't recognize the INTEGRATOR name you entered, which was %s\nNow exiting...\n",argv[1]);
    delete model;
    return -1;
  }
  
    
  const double A = ( (1 - omega*omega) * gamma) / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega); //convenient parameter for analytic solution
  const double B = 2 * beta * omega * gamma / ( pow(1-omega*omega,2) + 4 * beta*beta * omega*omega);//convenient parameter for analytic solution
  
  const double c1 = x[0] - A; //parameter, based on x, for analytic solution
  const double c2 = (x[1] + beta * (x[0] - A) - omega * B)/sqrt(1 - beta * beta); //parameter, based on x and x', for analytic solution
  
  double error_norm = 0.0; //This variable keeps track of the error norm as the integration progresses
  
  double t=0.; //The time
  
  error_norm += dt * pow( (x[0] - x_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2) * pow( (x[1] - x_prime_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2);  //For the zeroeth step of the iteration, starting to calcualte the error_norm.
  //this likely will add 0 to the error_norm, since the itegration hasn't even started yet and so
  //there has been no divergence between itegration and actual.
  
  for(int i = 0; i < nsteps; ++i){
    integrator->Step(t, x); //step the integration forward
    t = (i+1) * dt; //update the time
    error_norm += dt * pow( (x[0] - x_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2) * pow( (x[1] - x_prime_of_linear_osc(t,beta,gamma,omega,c1,c2,A,B)), 2); //update the error norm
  }
  
  printf("%lf\n",sqrt(error_norm)); //prints out the error norm
  
  delete integrator;
  delete model;
  return 0;
}
