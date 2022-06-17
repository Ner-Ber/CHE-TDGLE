import numpy as np
from scipy import signal

class Grid:
    def __init__(self, base_grid:np.ndarray):
        self.base_grid = base_grid
        pass

    def apply_dbl_laplacian(self, target_grid):
        # grid_dbl_laplaced = signal.convolve2d(target_grid, DBL_LAPLACIAN, mode='same', boundary='symm')
        grid_dbl_laplaced = self.apply_laplacian(self.apply_laplacian(target_grid))
        return grid_dbl_laplaced

    def apply_laplacian(self, target_grid):
        x, h = np.linspace(0, 1, target_grid.shape[0], endpoint=False, retstep=True)
        grid_laplaced = signal.convolve2d(target_grid, LAPLACIAN, mode='same', boundary='symm')/h**2
        return grid_laplaced
    
    def apply_polinomial_2grid(self, target_grid, polynom):
        return polynom(target_grid)
        pass


    def return_CH_rhs(self, base_grid, polynom, D):
        """will return numerical the right-hand-side of CH equation"""
        # return self.apply_laplacian(self.apply_polinomial_2grid(base_grid, polynom)) - D*self.apply_dbl_laplacian(base_grid)
        return self.apply_laplacian(self.apply_polinomial_2grid(base_grid, polynom)) - D*self.apply_laplacian(base_grid)
    
    def return_TDGLE_rhs(self, base_grid, polynom, D):
        """will return numerical the right-hand-side of TDGLE"""
        return D*self.apply_laplacian(base_grid) + self.apply_polinomial_2grid(base_grid, polynom)

    def return_CH_step(self, polynom, D, dt):
        # rhs = self.apply_laplacian(self.apply_polinomial_2grid(self.base_grid, polynom)) - D*self.apply_dbl_laplacian(self.base_grid)
        next_step_grid = dt*self.return_CH_rhs(self.base_grid, polynom, D) + self.base_grid
        # next_step_grid[next_step_grid>1] = 1
        # next_step_grid[next_step_grid<-1] = -1
        return Grid(next_step_grid)
    
    def return_TDGLE_step(self, polynom, D, dt):
        # rhs = self.apply_laplacian(self.apply_polinomial_2grid(self.base_grid, polynom)) - D*self.apply_dbl_laplacian(self.base_grid)
        # rhs = D*self.apply_laplacian(self.base_grid) + self.apply_polinomial_2grid(self.base_grid, polynom)
        next_step_grid = dt*self.return_TDGLE_rhs(self.base_grid, polynom, D) + self.base_grid
        # next_step_grid[next_step_grid>1] = 1
        # next_step_grid[next_step_grid<-1] = -1
        return Grid(next_step_grid)


    

LAPLACIAN = np.array([
    [0,1,0],
    [1,-4,1],
    [0,1,0],
])

DBL_LAPLACIAN = np.array([
    [1,-4,1],
    [-4,-2,-4],
    [1,-4,1],
])