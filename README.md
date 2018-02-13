# python_ODEsolver
Made to try to simplify ODE-solving by using dictionaries. Should translate and export equations fairly easily to LaTeX. 

Basics: 
The goal is to be able to quickly set up and solve a compounded ODE (i.e. multiple dimensions) 
when you only have the general formulas for each step. 

What's included: 
- A general ODE-solver that uses RungeKutta (2 or 4) , or ForwardEuler 
- a function-generator that creates an f(u,t) for the ODE-solver class 
  from a set of descriptions in a dictionary. 
 
This should work for most easily set up problems, as long as you can describe 
the differential equations they operate under. 
i.e. 
z(t) = v_z(t) + b 
v_z(t) = F_z/m
becomes: 
'z' : [z_initial_value, 'v_z(t) + b(t)'] 
'v_z' : [v_initial_value, 'F(t)/%s'%m ]
and a general formula for the force: 
'F(t)' : '%s' %(m * g)

You can pass multiple variables as well as these will simply become an 
anonymous function (lambda (x,y,z,...): return_value) 

From here on I'll provide 3 examples you can test out or be inspired by: 
- Particle on an elastic band/spring. 
- The Lorenz attractor (with animation suggestion) 
- SIZR (https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model
        with Zombies) 

#-------------------------- Example: ---------------------------------- # 
# Particle on a Spring. 
g = 9.81 # gravitatinoal force. 
m = 0.1  # mass of particle
k = 10.  # spring-constant. 
b = 0.08 # friction. 


z_0 = 0.1 # initial position 
v_0 = 0   # initial velocity

T = 50    # time (seconds) 
dt = 0.01 # time-step 

Variable_Equations = { 
            'z' : [z_0, 'v'], 
            'v' : [v_0, 'g(t) + r(v) + k(z)]
            } 
# This translates to mean z updates every time-step by whatever v is . 
# and v updates as the sum of forces acting on it. 
 
you_want_F = false #See below
Acceleration_Equations = {
            'g(t)': '-%s' %g , 
            'r(v)': '%s * (- v)' % (b/m) , 
            'k(z)': '%s * (- z)' % (k/m) 
              }
# These will not be stored and returned as they are only used to calculate 
# the update of our variables (the ones we want to track i.e. z and v)

# This adds another Force that only acts in some time-interval. 
# should also be possible to use a PiecewiseConstant, but as noted in the code, 
# I don't think I tested that. 
if you_want_F: 
  Variable_Equations['v'][1] += ' + F(t)'
  Acceleration_Equations['F(t)'] ='(%s) * cos(%s * t) if (t > %s and t < %s)  else 0' %((F/m),omega_f,F_start,F_stop)


f = FromDict([Variable_Equations,  Acceleration_Equations], T) 

# Creates a f(u,t), which is the general formula for the complete ODE. 
# that is; it contains all the updates to each value, and computes the 
# change for each time-step up to T. (The T-argument may be outdated from 
# an earlier version, I'm not actually sure why I pass it anymore. been using 
# this on and off for about the last 3 years)

#Next we set up a solver. 
solver = RungeKutta4(f) # generally works well, but there are other solvers too. 
solver.set_initial_conditions(f.ic) # could be more verbosely named, but I'm lazy. 
time_ = numpy.linspace(0,T,T/dt)     # should give you enough steps

u,t = solver.solve(time_) #solves the equations and returns to u. 
#Note: At this point the order of the variables is all jumbled up, 
#      so if we want to use them later on, we will need to extract the 
#      variable-names by running: 

module = sys.modules[__name__]
for val in range(len(f.key_list)):
	setattr(module, f.key_list[val], u[:,val])
# This will now allocate z and v as variables which we can use to plot, or analyse further. 
#i.e. 
from Your_favourite_plotter import plot  
plot(t, z, legend='z(t)', xlabel='t', ylabel='z(t) / v(t)')
plot(z, v, legend='phase-plot',xlabel='z', ylabel='v' )
#------------------------ end example ------------------------------__# 

Some notes: you can easily add or remove coefficients, or alter the 
way it acts by appending the string describing the behaviour to match what you 
want to model. However that (of course) means you'll have to generate the f(u,t) all over 
again. 
#------------------------  Example: ----------------------------------- #
#Simulate chaos: 
#Lorenz Attractor:
Lorenz_eqs = {
            'sigma(t)': 10.0,
            'beta(t)': '8.0/3',
            'rho(t)': 28.0
            }
Lorenz_vars  ={
            'x': [1, 'sigma(t)*(y - x)'],
            'y': [1, 'x * ( rho(t) - z ) - y'],            
            'z': [1, 'x*y - beta(t)*z']
            } 
Lorenz = [Lorenz_eqs, Lorenz_vars]
T = 20
dt = 0.2
f = (Lorenz, T)
solver = RungeKutta4(f)
solver.set_initial_conditions(f.ic)  
time_ = numpy.linspace(0,T,T/dt)     

u,t = solver.solve(time_)
module = sys.modules[__name__]
for val in range(len(f.key_list)):
	setattr(module, f.key_list[val], u[:,val])

#To see an animation of the return run something like: 
from mpl_toolkits.mplot3d import Axes3D
pylab.ion()
limits = (- np.max(u),np.max(u))
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False,
                          xlim=limits,
                          ylim=limits, 
                          zlim=limits,
                          projection = '3d')
pylab.show()
chaosline = ax.plot(x[0:1], y[0:1], z[0:1])
print 'animation ready, return to play'

raw_input('RETURN to run anumation')
for i in xrange(0,len(t),100):
  for line in (chaosline):    
    line.set_data(x[0:i],y[0:i])
    line.set_3d_properties(z[0:i])
    pylab.draw()
  pylab.draw()

raw_input('end')
# -------------------------- End Example ---------------------------- #


# ------------------------- SIZR- model. -----------------------------# 
##---------------SIZR - model, Disease and Zombies-------------------##
# Here's just a quick-setup of the dictionary for running a SIR-model 
# with Zombies included. 
# S - Susceptibles. 
# I - Infected. 
# Z - Zombies 
# R - Recovered (aka. Dead zombies) 

SIZR_Parameters = {'alpha(t)':'0.0016',
                'rho(t)': '1.0',
                'beta(t)': '0.003',
                'dS(t)': 0.001,
                'dI(t)': 0.014,
                'p(t)': '0.1 if t <= 15 and t >= 6 else 0',
                'nu(t)': 0.1,
                'summa(t)': 2,}               

SIZR_Equations = {'S': [10,  'summa(t) - beta(t)*S*Z - dS(t)*S'], 
                  'I': [0,   'beta(t)*S*Z - rho(t)*I - dI(t)*I'],
                  'Z': [100, 'rho(t)*I - alpha(t)*S*Z'],
                  'R': [0,   'dS(t)*S + dI(t)*I + alpha(t)*S*Z']}
SIZR = [SIZR_Parameters, SIZR_Equations]

solver = RungeKutta4(f) 
solver.set_initial_conditions(f.ic) 
time_ = numpy.linspace(0,T,T/dt)    
u,t = solver.solve(time_) 
module = sys.modules[__name__]
for val in range(len(f.key_list)):
	setattr(module, f.key_list[val], u[:,val])
  
# Plot using plot(t,S); plot(t,I); plot(t,Z); plot(t,R) #or similar 
# to see the effect of a zombie outbreak and the lagging effect of 
# the contagion. of course set your appropriate T and dt here. 
# something like T = 60, dt = 0.1 should do. 
