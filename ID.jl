"""
Assign global degrees of freedom (DOF) to mesh nodes.

Parameters:
- p: Node coordinates and boundary labels array
     p[:,1:2] contains x,y coordinates
     p[:,3] contains boundary labels:
     1: Bottom boundary
     2: Left boundary  
     3: Top boundary
     4: Right boundary
     5: Interior node

Returns:
- DOFS: Array containing global DOF number for each node
        0 indicates a constrained DOF (e.g., Dirichlet boundary)
"""
function ID(p)
    # Boundary condition types
    D = [1 2 3 4]  # Dirichlet boundaries
    P = [0 0]      # Periodic boundary pairs (none by default)
    
    # Initialize DOF array and counter
    DOFS = zeros(size(p,1), 1)
    ind = 1
    
    # Loop over all nodes
    for i = 1:size(p,1)
        # Check if node is on Dirichlet boundary
        if !isempty(findall(x->x!=0, vec(p[i,3] .== D)))
            DOFS[i] = 0  # Constrained DOF
        else
            # Check for periodic boundary conditions
            if P == [1 3]  # If bottom-top periodicity is enabled
                if p[i,3] == 3  # Top boundary node
                    # Find matching node on bottom boundary
                    inds1 = findall(p[:,3] .== 1)
                    xind = findall(abs.(p[i,1] .- p[inds1,1]) .< 1e-10)
                    if !isempty(xind)
                        # Use same DOF as bottom node
                        DOFS[i] = DOFS[inds1[xind[1]]]
                        continue
                    end
                end
            end
            # Assign new DOF number
            DOFS[i] = ind
            ind += 1
        end
    end
    
    return DOFS
end

