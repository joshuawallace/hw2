integrators = euler.o
equations = duffing.o
objects = duffing_solve.o $(integrators) $(equations)

integrators2 = euler.o adams-bashforth.o runge-kutta.o
equations2 = duffing.o lorenz.o linear-oscillator.o
objects2 = ode_solve.o $(integrators2) $(equations2)

CXXFLAGS = -g -Wall

all: duffing_solve ode_solve test_convergence

duffing_solve : $(objects)
	$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend

ode_solve: $(objects2)
	$(CXX) -o $@ $^

test_convergence: test_convergence.o linear-oscillator.o $(integrators2)
	$(CXX) -o $@ $^
