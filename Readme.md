# Spectral Element Method - 2D Poisson Solver

A Julia implementation of the Spectral Element Method (SEM) for solving the 2D Poisson equation:

-∇²u = f in Ω

with Dirichlet boundary conditions.

## Overview

This solver uses high-order spectral elements with Gauss-Legendre quadrature for accurate numerical solutions. The implementation features:

- Higher-order polynomial basis functions (configurable order r)
- Gauss-Legendre quadrature for numerical integration
- Structured rectangular mesh generation
- Support for Dirichlet and periodic boundary conditions
- Error computation against known analytical solutions

## Files

### Core Components
- `Poisson.jl` - Main solver implementation
- `Mesh2Drect.jl` - Rectangular mesh generator
- `intNodes.jl` - Internal node generation for higher-order elements
- `ID.jl` - Global degree of freedom assignment

### Basis Functions
- `basis.jl` - Lagrange polynomial basis functions
- `basisx.jl` - Derivatives of basis functions
- `GLpw.jl` - Gauss-Lobatto points and weights computation

### Matrix Assembly
- `KMtrx.jl` - Stiffness matrix assembly
- `MMtrx.jl` - Mass matrix assembly
- `intWeights.jl` - Integration weights computation

### Utilities
- `ChkNodeEx.jl` - Node existence checker
- `bubble_sort.jl` - Sorting utility
- `f_ex.jl` - Source term function
- `ex_sol.jl` - Exact solution for error computation
- `femerror.jl` - Error calculation functions

## Usage

1. Set problem parameters in `Poisson.jl`:
   ```julia
   # Problem size
   Nx = 4  # Number of elements in x direction
   Ny = 4  # Number of elements in y direction
   r = 4   # Polynomial order
   
   # Domain size
   Lx = 1.0  # Length in x direction
   Ly = 1.0  # Length in y direction
   ```

2. Run the solver:
   ```julia
   include("Poisson.jl")
   ```

3. The solution and error metrics will be computed and displayed.

## Example Results



## License

This code is released under the MIT License. See LICENSE file for details.
