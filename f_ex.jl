"""
Source term function for the Poisson equation.
Represents the right-hand side f in -∇²u = f.

Parameters:
- n: Oscillation frequency parameter
- x: x-coordinate
- y: y-coordinate

Returns:
- f(x,y) = 2n²π² sin(nπx)sin(nπy)
"""
function f_ex(n, x, y)
    return 2*sin.(n*pi*x).*sin.(n*pi*y)*n^2*pi^2
end
	