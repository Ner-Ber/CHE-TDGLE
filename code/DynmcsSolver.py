import numpy as np
from  Grid import Grid
from tqdm import tqdm


class DynmcsSolver:
    def __init__(self, initial_condition, dt, Nsteps, polynom, D, equation_name):
        # set values to instance
        self.initial_condition = initial_condition
        self.Nsteps = Nsteps
        self.eq_name = equation_name
        self.polynom = polynom
        self.dt = dt
        self.D = D
        self.prep4sol()

    def prep4sol(self):
        # create grounds for solving
        # self.results_array = np.empty(([self.Nsteps] + list(self.initial_condition.shape) ))
        self.results_array = np.empty( [self.Nsteps] + list(self.initial_condition.shape) )
        self.results_array[0] = self.initial_condition
        self.start_grid = Grid(self.initial_condition)
        self.t_vec = range(1, self.Nsteps)
    
    def solve(self):
        """ a simple and straigh forward numeric solution (probably unstable)  """
        current_grid = self.start_grid
        for t in tqdm(self.t_vec):
            if 'GL' in self.eq_name:
                current_grid = current_grid.return_TDGLE_step(self.polynom, self.D, self.dt)
            elif self.eq_name=='CH':
                current_grid = current_grid.return_CH_step(self.polynom, self.D, self.dt)
            self.results_array[t] = current_grid.base_grid

    def solver_rungekutta(self):
        """ a fourthe order Runge-Kutta numerical solver.
        see: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods"""
        
        assert (('GL' in self.eq_name) | (self.eq_name=='CH')) , 'inadequate equatio name'

        current_solution = self.start_grid
        for t in tqdm(self.t_vec):
            #-- select current rhs function
            if 'GL' in self.eq_name:
                rhs_function = lambda X: current_solution.return_TDGLE_rhs(X, self.polynom, self.D)
            elif self.eq_name=='CH':
                rhs_function = lambda X: current_solution.return_CH_rhs(X, self.polynom, self.D)

            #-- create 4th order RK step
            current_step = current_solution.base_grid
            #--- crop to boundaries of +/-1
            # current_step[current_step>1] = 1
            # current_step[current_step<-1] = -1
            #--- normalize betwen 1 and -1
            # current_step -= current_step.min()
            # current_step /= current_step.max()
            # current_step *= 2
            # current_step -= 1
            
            k1 = (rhs_function(current_step));
            k2 = (rhs_function(current_step+self.dt*0.5*k1));
            k3 = (rhs_function(current_step+self.dt*0.5*k2));
            k4 = (rhs_function(current_step+self.dt*k3));
            next_step = current_step + (self.dt/6)*(k1 + 2*k2 + 2*k3 + k4)
            self.results_array[t] = next_step
            current_solution = Grid(next_step)

        pass





if __name__ == '__main__':

    axis_vec = np.arange(-100,100)
    X,Y = np.meshgrid(axis_vec, axis_vec)
    R = 15
    init_cond = (X**2 + Y**2 < R**2).astype(np.float)
    # init_cond = np.random.rand(20,20)


    dt = 1e-3
    Nsteps = int(1e3)
    # polynom = lambda x: x*(1-x**2)
    polynom = lambda u: u**3-u
    D = 1
    equation_name = 'GL'
    solver = DynmcsSolver(init_cond,  dt,  Nsteps,  polynom,  D,  equation_name)

    solver.solve()



    pass