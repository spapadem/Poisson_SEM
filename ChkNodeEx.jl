"""
Check if a node already exists in the mesh within a given tolerance.

Parameters:
- p: Coordinates of node to check [x, y, label]
- pr: Array of existing node coordinates and labels
- ind: Current number of nodes in pr
- hx: Mesh spacing in x direction (used for tolerance)
- hy: Mesh spacing in y direction (used for tolerance)

Returns:
- AlEx: Boolean indicating if node exists (1 if exists, 0 if not)
- indx: Index of matching node if found, 0 otherwise
"""
function ChkNodeEx(p, pr, ind, hx, hy)
    # Initialize return values
    indx = 0
    AlEx = 0
    
    # Extract coordinates of node to check
    px = p[1]
    py = p[2]
    
    # Set tolerance based on mesh spacing
    xtol = 1e-6 * hx
    ytol = 1e-6 * hy
    
    # Find nodes matching x coordinate within tolerance
    ind1 = findall(x -> x .< xtol, vec(abs.(px .- pr[1:ind-1,1])))
    # Find nodes matching y coordinate within tolerance
    ind2 = findall(y -> y .< ytol, vec(abs.(py .- pr[1:ind-1,2])))
    
    # Check for nodes matching in both x and y
    flag = 0
    for m = 1:length(ind1)
        for n = 1:length(ind2)
            if ind1[m] == ind2[n]
                AlEx = 1
                indx = ind1[m]
                flag = 1
                break
            end
        end
        if flag == 1
            break
        end
    end
    
    return [AlEx, indx]
end
