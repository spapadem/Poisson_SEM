"""
Compute error metrics for finite element solution.

Parameters:
- nodes: Array of node coordinates and labels
- el: Element connectivity array
- uh: Computed finite element solution
- hx: Element size in x direction
- hy: Element size in y direction
- r: Order of polynomial basis
- osc_f: Oscillation frequency parameter for exact solution

Returns:
- errt: L2 norm of error between computed and exact solution
- nu: L2 norm of exact solution
- nuh: L2 norm of computed solution
"""
function femerror(nodes, el, uh, hx, hy, r, osc_f)
    # Get Gauss-Legendre points and weights
    glp, glw = GLpw(r)
    
    # Get quadrature weights
    gj = intWeights(r, 2, 2)
    
    # Initialize error metrics for each element
    err = zeros(size(el,1), 1)  # Error norm
    nu = zeros(size(el,1), 1)   # Exact solution norm
    nuh = zeros(size(el,1), 1)  # Computed solution norm
    
    # Loop over elements
    for l = 1:size(el,1)
        # Get element nodes
        pt = nodes[el[l,:], 1:2]
        
        # Integrate over element using quadrature
        for j = 1:(r+1)^2
            xp = pt[j,1]
            yp = pt[j,2]
            indx = el[l,j]
            
            # Compute error at quadrature point
            err[l] += 0.25*hx*hy*gj[j] * (ex_sol(osc_f,xp,yp) - uh[indx])^2
            
            # Compute norms at quadrature point
            nu[l] += 0.25*hx*hy*gj[j] * (ex_sol(osc_f,xp,yp))^2
            nuh[l] += 0.25*hx*hy*gj[j] * (uh[indx])^2
        end
    end
    
    # Compute global error metrics
    errt = sqrt(sum(err))  # Total error norm
    nu = sqrt(sum(nu))     # Total exact solution norm
    nuh = sqrt(sum(nuh))   # Total computed solution norm
    
    return errt, nu, nuh
end

