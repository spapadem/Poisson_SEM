"""
Generate internal nodes for higher-order elements.

Parameters:
- el: Element connectivity array
- p: Initial node coordinates and labels array
- r: Order of polynomial basis
- hx: Element size in x direction
- hy: Element size in y direction
- x: Vector of x-coordinates defining mesh divisions
- y: Vector of y-coordinates defining mesh divisions

Returns:
- pr: Array of all nodes (including internal nodes) with coordinates and labels
- elr: Updated element connectivity array including internal nodes
"""
function intNodes(el, p, r, hx, hy, x, y)
    Nel = size(el,1)  # Number of elements
    N = size(p,1)     # Number of initial nodes
    
    # Get Gauss-Legendre points for internal node placement
    glp, = GLpw(r)
    
    # Initialize arrays for new nodes
    pr = zeros(Nel*(r+1)^2, 3)  # All nodes including internal
    elr = zeros(size(el,1), (r+1)^2)  # Updated connectivity
    
    ind = 1  # Counter for new nodes
    
    # Loop over elements
    for l = 1:Nel
        # First handle vertices of element
        for indel = 1:4
            p1 = p[Int(el[l,indel]),:]
            # Check if node already exists
            AlEx, ind1 = ChkNodeEx(p1, pr, ind, hx, hy)
            if AlEx == 0
                # Add new node
                pr[ind,:] = p1
                elr[l,indel] = ind
                ind += 1
            else
                # Use existing node
                elr[l,indel] = ind1
            end
        end
        
        # Now add internal nodes
        indel = 5
        # Loop over Gauss points in y direction
        for k = 1:r+1
            # Loop over Gauss points in x direction
            for j = 1:r+1
                # Get coordinates of bottom-left vertex
                pt = p[Int(el[l,1]),1:2]
                # Map reference element coordinates to physical space
                px = 0.5*hx*(glp[j]+1) + pt[1]
                py = pt[2] + 0.5*hy*(glp[k]+1)
                
                # Check if node already exists
                AlEx, ind1 = ChkNodeEx([px py], pr, ind, hx, hy)
                if AlEx == 0
                    # Add new node
                    pr[ind,1:2] = [px py]
                    elr[l,indel] = ind
                    indel += 1
                    
                    # Assign boundary label
                    if px == x[1]
                        plab = 2      # Left boundary
                    elseif px == x[end]
                        plab = 4      # Right boundary
                    elseif py == y[1]
                        plab = 1      # Bottom boundary
                    elseif py == y[end]
                        plab = 3      # Top boundary
                    else
                        plab = 5      # Interior node
                    end
                    pr[ind,3] = plab
                    ind += 1
                else
                    # Use existing node if not already in element
                    if !any(ind1 .== elr[l,:])
                        elr[l,indel] = ind1
                        indel += 1
                    end
                end
            end
            
            # Add nodes along element edges
            for i = 1:r+1
                # Right edge
                pt = p[Int(el[l,2]),1:2]
                px = pt[1] - 0.5*hx*(glp[k]+1)
                py = 0.5*hy*(glp[i]+1) + pt[2]
                
                AlEx, ind1 = ChkNodeEx([px py], pr, ind, hx, hy)
                if AlEx == 0
                    pr[ind,1:2] = [px py]
                    elr[l,indel] = ind
                    indel += 1
                    
                    # Assign boundary label
                    if px == x[1]
                        plab = 2
                    elseif px == x[end]
                        plab = 4
                    elseif py == y[1]
                        plab = 1
                    elseif py == y[end]
                        plab = 3
                    else
                        plab = 5
                    end
                    pr[ind,3] = plab
                    ind += 1
                else
                    if !any(ind1 .== elr[l,:])
                        elr[l,indel] = ind1
                        indel += 1
                    end
                end
            end
            
            # Similar process for top and left edges...
            # (Code continues with similar pattern for remaining edges)
        end
    end
    
    # Trim unused entries from pr array
    pr = pr[1:ind-1,:]
    
    return pr, elr
end
