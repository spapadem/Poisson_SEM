"""
Compute the value of a Lagrange polynomial basis function.

Parameters:
- a: Index of shape function (1 to r+1)
- j: Index of evaluation point in the set of Gauss-Legendre points [-1,1]
- r: Order of polynomial basis

Returns:
- Na: Value of the a-th Lagrange polynomial at point x_j

Note: Uses Gauss-Legendre points as interpolation nodes
"""
function basis(a, j, r)
    # Get Gauss-Legendre points
    glp, = GLpw(r)
    
    # Initialize basis function value
    Na = 1
    
    # Compute Lagrange polynomial using product formula
    for p = 1:r+1
        if p == a
            continue
        end
        # Na = ∏(x - x_p)/(x_a - x_p) for p ≠ a
        Na = Na * (glp[j] - glp[p])/(glp[a] - glp[p])
    end
    
    return Na
end
