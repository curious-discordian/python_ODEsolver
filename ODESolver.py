
"""
Modified ODESolver superclass, 
last modified 06.05.2017 ?
Addon 06.05.17:
Modify to make universal ODESolver access

Addons (13.02.2017) :
Added FromDict option, for creating f(u,t) from dictionaries. 
Based on the earlier ProblemGEN, but without odespy. 
Rewritten in to use the dictionaries for creating 
general f(u,t) type problems. Should still work for 
PDE- though solver class in ODESolver needs rewriting 
for that specific purpose. 

Also, figured out the trick;
    #retrieve variables from within the Problem-class(f) and make global: 

module = sys.modules[__name__]
for val in range(len(f.key_list)):
    setattr(module, f.key_list[val], u[:,val])


useful for later manipulation of data after ODE is solved, e.g. plotting.

Addons (28.12.2014) :
Timer, which shows quadrants of completion. 
augmented with 'clear' command, for clean display. 

"""

###Note to self: 
"""
Use snippet below to export variables by name. 

module = sys.modules[__name__]
for val in range(len(f.key_list)):
    setattr(module, f.key_list[val], u[:,val])

"""

import numpy as np

class ODESolver:
    """
    Superclass for numerical methods solving scalar and vector ODEs

      du/dt = f(u, t)

    Attributes:
    t: array of time values
    u: array of solution values (at time points t)
    k: step number of the most recently computed solution
    f: callable object implementing f(u, t)
    timer: if 'on' shows percentage points 25, 50, 75, 99 complete
    """
    def __init__(self, f):      #, timer='off'):
        if not callable(f):
            raise TypeError('f is %s, not a function' % type(f))
        # For ODE systems, f will often return a list, but
        # arithmetic operations with f in numerical methods
        # require that f is an array. Let self.f be a function
        # that first calls f(u,t) and then ensures that the
        # result is an array of floats.
        self.f = lambda u, t: np.asarray(f(u, t), float)
        #self.timer = timer 
        
    def advance(self):
        """Advance solution one time step."""
        raise NotImplementedError

    def set_initial_condition(self, U0):
        if isinstance(U0, (float,int)):  # scalar ODE
            self.neq = 1
            U0 = float(U0)
        else:                            # system of ODEs
            U0 = np.asarray(U0)          # (assume U0 is sequence)
            self.neq = U0.size
        self.U0 = U0

        # Check that f returns correct length:
        try:
            f0 = self.f(self.U0, 0)
        except IndexError:
            raise IndexError('Index of u out of bounds in f(u,t) func. Legal indices are %s' % (str(range(self.neq))))
        if f0.size != self.neq:
            raise ValueError('f(u,t) returns %d components, while u has %d components' % (f0.size, self.neq))

    def solve(self, time_points, terminate=None):
        """
        Compute solution u for t values in the list/array
        time_points, as long as terminate(u,t,step_no) is False.
        terminate(u,t,step_no) is a user-given function
        returning True or False. By default, a terminate
        function which always returns False is used.
        """
        if terminate is None:
            terminate = lambda u, t, step_no: False

        if isinstance(time_points, (float,int)):
            raise TypeError('solve: time_points is not a sequence')
        if time_points.size <= 1:
            raise ValueError('ODESolver.solve requires time_points array with at least 2 time points')

        self.t = np.asarray(time_points)
        n = self.t.size
        if self.neq == 1:  # scalar ODEs
            self.u = np.zeros(n)
        else:              # systems of ODEs
            self.u = np.zeros((n,self.neq))

        # Assume that self.t[0] corresponds to self.U0
        self.u[0] = self.U0
        #set up for clear screen
        #Loc.1
        
        # Time loop
        a,b,c,stop = False, False,False,False
        for k in range(n-1):
            self.k = k
            self.u[k+1] = self.advance()

            #Removed Timer. Loc.2 
                
            if terminate(self.u, self.t, self.k+1):
                break  # terminate loop over k
        return self.u[:k+2], self.t[:k+2]


class ForwardEuler(ODESolver):
    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        u_new = u[k] + dt*f(u[k], t[k])
        return u_new

class RungeKutta2(ODESolver):
    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        dt2 = dt/2.0
        K1 = dt*f(u[k], t[k])
        K2 = dt*f(u[k] + 0.5*K1, t[k] + dt2)
        u_new = u[k] + K2
        return u_new

class RungeKutta4(ODESolver):
    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        dt2 = dt/2.0
        K1 = dt*f(u[k], t[k])
        K2 = dt*f(u[k] + 0.5*K1, t[k] + dt2)
        K3 = dt*f(u[k] + 0.5*K2, t[k] + dt2)
        K4 = dt*f(u[k] + K3, t[k] + dt)
        u_new = u[k] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
        return u_new


from numpy import *
class FromDict: 
    """
    f(u,t) Generator for solving ODE. 
    
    It will generate the f(u,t) from a set of dictionaries as 
    described in __init__. 
    """

    def __init__(self, Formula_list_of_dict , T, Show_Identified=True,\
                 Starttime=0):
        """ 
        
        Variables format:
            'Variable_Name' : ['Initial Condition', 'Differential Formula']
        Parameters format:
            'Parameter(variables)': 'value or formula'


        Keywords: 
          Show_Identified, bool - shows the interpreted input. 
          Starttime, int/float - for starting at different value t.
                                useful for chaining ODE using terminate().
        Usage example: 
            Variables = {
                        'x' : [0, 'vx'],
                        'vx' : [1, '10 - P(x,vx)']
                        }
            Parameters = {
                    'P(x,vx)' : '1/x * (1 + vx**2)' 
                        }
            
            Formula = [Parameters, Variables]
            problem = FromDict(Formula, T=30, Show_Identified=True, Starttime=0)
            
        returns: u,t (list of indices for u[i] from 'Show_Identified' parameter)

        Use snippet below to export variables by name. 

        module = sys.modules[__name__]
        for val in range(len(f.key_list)):
            setattr(module, f.key_list[val], u[:,val])

        """
        Parameters_Dict = Formula_list_of_dict[0]
        Variables_Dict = Formula_list_of_dict[1]
        parameter_list = []
        if Show_Identified:
            print 'Identified Parameters:'
        for factor in Parameters_Dict:
            if Show_Identified:
                print '%s : %s'\
                   %(factor, Parameters_Dict[factor])
            param_start = factor.find('(')
            parameter_list.append(str(factor[:param_start+1]))
            # PiecewiseConstant not tested. 
            if isinstance(Parameters_Dict[factor], str):
                if 'PiecewiseConstant' in Parameters_Dict[factor]: 
                    start = factor.find('(')
                    stop = factor.find(')')
                    variables = factor[start+1:stop]
                    variables.replace(',', ' ')
                    setattr(self, factor[:start],
                        eval(Parameters_Dict[factor]))

            #----Parameter settings-----#
                else:
                    start = factor.find('(')
                    stop = factor.find(')')
                    variables.replace(',', ' ')
                    variables = factor[start+1:stop]
                    setattr(self, factor[:start],
                            eval('lambda %s: %s'
                            % (variables, Parameters_Dict[factor]))) 
                  
            else: 
                start = factor.find('(')
                stop = factor.find(')')
                variables.replace(',', ' ')
                variables = factor[start+1:stop]
                setattr(self, factor[:start],
                         eval('lambda %s: %s'
                         % (variables, Parameters_Dict[factor]))) 

        key_list = []
        initial_cond_list = []
        formula_list = []
        for key in Variables_Dict:
            key_list.append(str(key))
            formula_list.append(Variables_Dict[key][1])
            initial_cond_list.append(str(key)+'0')
            setattr(self, str(key)+'0', float(Variables_Dict[key][0]))
   
        self.key_list = key_list 
        self.initial_conditions = initial_cond_list 
        self.T = T
        self.Starttime= Starttime #may be obsolete
        
        
        # appending parameters with self.prefix 
        # if they're in the parameter list
        # and creating the list of formulas.
 
        formula_list_mod = []
        for formula in formula_list:
            formula_elements = formula.split()
            inner_formula_elements = []

            for element in formula_elements:  

                if any([element.startswith(x) for x in parameter_list]):
                    inner_formula_elements.append('self.%s'%(element))           
                else:
                    inner_formula_elements.append(element)

            inner_formula = ['']
            for i in inner_formula_elements:
                inner_formula[0] += i #puts all into one.
                inner_formula[0] += ' ' #added whitespace, reads easier. 
            formula_list_mod += inner_formula #extracts to the list.

        self.formula_list = formula_list_mod #formulas in same order
        self.ic = list()
        #Snazz: Just to show what's been interpreted 
        #       and how the formulas look 
        if Show_Identified:
            print 'Identified Variables and formulas:'
            for i in range(len(key_list)):
                print key_list[i], ':', formula_list_mod[i]
            for n in self.initial_conditions: 
                print n, ':', eval('self.%s'%n)
                self.ic.append(eval('self.%s'%n))

        ## ------ Try to set up __call__ underneath here.-----##
        """
        for key in range(len(key_list)):
            name = key_list[key]
            setattr(self.__call__, locals()[str(name)], u[key])
            #setattr(self.__call__, 
            #        eval('vars()[key_list[key]]'+'_new'), formula_list[key])
            #setattr(
        """
    def __call__(self,u,t):
        """
        Call works as a "standard" f(u,t)
        NOTE: Try to set this up otherwise, creating the Call functino
              in the __init__, Would make faster ? 
        """

        # extracts the names of variables for the formulas
        for key in range(len(self.key_list)):
            name = self.key_list[key]
            value = u[key]
            vars()[name] = value      
          
        # generates next instance. 
        listed_returns = []
        for i_ in range(len(self.key_list)):                     
            listed_returns.append(eval(self.formula_list[i_]))
        
        return listed_returns  





class Derivative:
    def __init__(self, f, h=1E-9):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h      # make short forms
        return (f(x+h) - f(x-h))/(2*h)


registered_solver_classes = [
    ForwardEuler, RungeKutta2, RungeKutta4]

def test_exact_numerical_solution():
    a = 0.2; b = 3

    def f(u, t):
        return a + (u - u_exact(t))**5

    def u_exact(t):
        """Exact u(t) corresponding to f above."""
        return a*t + b

    U0 = u_exact(0)
    T = 8
    n = 10
    tol = 1E-15
    t_points = np.linspace(0, T, n)
    for solver_class in registered_solver_classes:
        solver = solver_class(f)
        solver.set_initial_condition(U0)
        u, t = solver.solve(t_points)
        u_e = u_exact(t)
        max_error = (u_e - u).max()
        msg = '%s failed with max_error=%g' % \
              (solver.__class__.__name__, max_error)
        print '%s: %s' %(solver_class, max_error < tol)
        assert max_error < tol, msg

if __name__ == '__main__':
    test_exact_numerical_solution()


"""
Skratches: 
    #From loc.1
    if self.timer == 'on':
            from os import system
        else: 
            pass
    #From loc.2
    #Timer; shows percentage complete:
            if self.timer == 'on': #use checkpoints, to minimize strain
                if (float(k)/n - 0.25) > 1e-10 and a == False :
                    system('clear')
                    print '25%'
                    a = True
                elif (float(k)/n - 0.5) > 1e-10 and b == False:
                    system('clear')
                    print '50%'
                    b = True
                elif (float(k)/n - 0.75) > 1e-10 and c == False:
                    system('clear')
                    print '75%'
                    c = True
                elif (float(k)/n - 0.99) > 1e-10 and stop == False:
                    system('clear')
                    print '99%'
                    stop = True
"""