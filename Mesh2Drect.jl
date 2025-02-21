"""
Generate a structured rectangular mesh.

Parameters:
- x: Vector of x-coordinates defining mesh divisions
- y: Vector of y-coordinates defining mesh divisions

Returns:
- p: Array of node coordinates and boundary labels
     p[:,1:2] = x,y coordinates
     p[:,3] = boundary label
- e: Array of edge connectivity and labels
     e[:,1:2] = node indices forming edge
     e[:,3] = edge label
- el: Array of element connectivity (4 nodes per element)

Note: Boundary/Edge labels:
1: Bottom boundary
2: Left boundary
3: Top boundary
4: Right boundary
5: Interior
"""
function Mesh2Drect(x, y)
    # Get number of elements in each direction
    Nx = length(x) - 1
    Ny = length(y) - 1
    
    # Create nodes
    p = zeros((Nx+1)*(Ny+1), 3)
    ind = 1
    for j = 1:Ny+1
        for i = 1:Nx+1
            # Set node coordinates
            p[ind,1:2] = [x[i] y[j]]
            
            # Assign boundary labels
            if x[i] == x[1]
                plab = 2      # Left boundary
            elseif x[i] == x[end]
                plab = 4      # Right boundary
            elseif y[j] == y[1]
                plab = 1      # Bottom boundary
            elseif y[j] == y[end]
                plab = 3      # Top boundary
            else
                plab = 5      # Interior node
            end
            p[ind,3] = plab
            ind += 1
        end
    end
    
    # Create edges
    e = zeros(4*Nx*Ny, 3)
    ind = 1
    
    # Horizontal edges
    for j = 1:Ny+1
        for i = 1:Nx
            e[ind,1:2] = [(j-1)*(Nx+1) + i, (j-1)*(Nx+1) + (i + 1)]
            ind += 1
        end
    end
    
    # Vertical edges
    for i = 1:Nx+1
        for j = 1:Ny
            e[ind,1:2] = [(j-1)*(Nx+1)+i, j*(Nx+1)+i]
            ind += 1
        end
    end
    
    # Assign edge labels
    for i = 1:ind-1
        ed = e[i,1:2]
        p1 = p[Int(ed[1]),1:2]
        p2 = p[Int(ed[2]),1:2]
        
        # Label edges based on location
        if p1[1] == x[1] && p2[1] == x[1]
            elab = 2      # Left boundary
        elseif p1[1] == x[end] && p2[1] == x[end]
            elab = 4      # Right boundary
        elseif p1[2] == y[1] && p2[2] == y[1]
            elab = 1      # Bottom boundary
        elseif p1[2] == y[end] && p2[2] == y[end]
            elab = 3      # Top boundary
        else
            elab = 5      # Interior edge
        end
        e[i,3] = elab
    end
    e = e[1:ind-1,:]
    
    # Create elements (counter-clockwise node ordering)
    Nel = Nx*Ny
    el = zeros(Nel, 4)
    for i = 1:Nel
        rem = floor((i-1)/Nx)
        el[i,:] = [i+rem, i+1+rem, i+Nx+2+rem, i+Nx+1+rem]
    end
    
    return p, e, el
end

    
    
    
    
