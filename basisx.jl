"""
Compute the derivative of a Lagrange polynomial basis function.

Parameters:
- a: Index of shape function (1 to r+1)
- j: Index of evaluation point in the set of Gauss-Legendre points [-1,1]
- r: Order of polynomial basis

Returns:
- Na_x: Value of the derivative of a-th Lagrange polynomial at point x_j

Note: Uses product rule to compute derivative of Lagrange polynomial
"""
function basisx(a, j, r)
    # Get Gauss-Legendre points
    glp, = GLpw(r)
    
    # Initialize derivative
    Na_x = 0
    
    # Apply product rule for derivatives
    for p = 1:r+1
        if p != a
            k = 1/(glp[a] - glp[p])
            # Multiply by other terms in product
            for m = 1:r+1
                if m != a && m != p
                    k = k * (glp[j] - glp[m])/(glp[a] - glp[m])
                end
            end
            Na_x = Na_x + k
        end
    end
    
    return Na_x
end
