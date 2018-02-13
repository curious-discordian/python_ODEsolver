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

From I'll provide 3 examples you can test out or be inspired by: 
- Particle on an elastic band/spring. 
- The Lorenz attractor (with animation suggestion) 
- SIZR (https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model
        with Zombies) 


