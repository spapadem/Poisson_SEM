"""
Assemble element mass matrix for spectral element method.

Parameters:
- r: Order of polynomial basis
- hx: Element size in x direction
- hy: Element size in y direction
- gj: Gaussian quadrature weights

Returns:
- M: Element mass matrix of size ((r+1)², (r+1)²)
"""
function MMtrx(r, hx, hy, gj)
    # Initialize mass matrix
    M = zeros((r+1)^2, (r+1)^2)
    
    # Get Gauss-Legendre points and weights
    glp, glw = GLpw(r)
    
    # Reference element vertices
    pref = [-1 -1 0;
             1 -1 0;
             1  1 0;
            -1  1 0]
    elref = [1 2 3 4]
    
    # Get internal nodes for reference element
    pref, elref = intNodes(elref, pref, r, 2, 2, [-1 1], [-1 1])
    elref = elref[1,:]
    
    # Loop over test functions
    for i = 1:size(M,1)
        # Find x,y indices for test function
        itomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[i]),1] .- glp)))
        itomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[i]),2] .- glp)))
        
        # Loop over trial functions
        for j = 1:size(M,2)
            jtomnx = findall(x -> x .< 1e-10, vec(abs.(pref[Int(elref[j]),1] .- glp)))
            jtomny = findall(y -> y .< 1e-10, vec(abs.(pref[Int(elref[j]),2] .- glp)))
            
            # Numerical integration using Gaussian quadrature
            for p = 1:(r+1)^2
                xp = pref[p,1]
                yp = pref[p,2]
                
                p2x = findall(x -> x.<1e-10, vec(abs.(xp .- glp)))
                q2y = findall(x -> x.<1e-10, vec(abs.(yp .- glp)))
                
                # Add contribution from basis function product
                M[i,j] = M[i,j] + 0.25*hy*hx*gj[p] * 
                    (itomnx[1] == p2x[1]) * (jtomnx[1] == p2x[1]) *
                    (itomny[1] == q2y[1]) * (jtomny[1] == q2y[1])
            end
        end
    end
    
    return M
end
