# Main file implementing Spectral Element Method (SEM) for 2D Poisson equation
# ∇²u = f in Ω
# with Dirichlet boundary conditions

using SparseArrays
using LinearAlgebra
using Printf
using Plots

"""
Solves the Poisson equation using Spectral Element Method (SEM) with higher-order elements.
Parameters:
- r: Order of the finite elements (r=1 for linear elements, r=2 for quadratic, etc.)
- hx: Grid spacing in x direction
- hy: Grid spacing in y direction
"""
# Import required FEM utility functions
include("basis.jl")        # Lagrange basis functions
include("basisx.jl")       # Derivatives of basis functions
include("bubble_sort.jl")  # Sorting utility
include("ChkNodeEx.jl")    # Node existence checker
include("ex_sol.jl")       # Exact solution for error computation
include("femerror.jl")     # Error calculation functions
include("f_ex.jl")         # Source term function
include("GLpw.jl")         # Gauss-Legendre quadrature points and weights
include("ID.jl")           # Global DOF assignment
include("intNodes.jl")     # Interior nodes generation
include("intWeights.jl")   # Integration weights
include("KMtrx.jl")        # Stiffness matrix assembly
include("Mesh2Drect.jl")   # Rectangular mesh generator
include("MMtrx.jl")        # Mass matrix assembly

# Domain parameters
b  = 1   # Width of the domain
xF = 1   # Length of the domain

# Discretization parameters
r = 3;
hx = 0.025;
hy = 0.025;

# Generate mesh points
x = Array(0 : hx : xF)    
y = Array(0 : hy : b)

osc_f = 1  # Oscillation frequency parameter for the source term

# Generate initial coarse mesh
println("Creating initial mesh.")
p_in,_,el_in = Mesh2Drect(x,y)  # Points and elements of coarse mesh

# Refine mesh based on element order r
println("Creating additional DOFs.")
p,el = intNodes(el_in,p_in,r,hx,hy,x,y)  # p: node coordinates, el: element connectivity
Nel = size(el,1)  # Number of elements

# Set boundary conditions
ind1 = findall(x -> x .== 1, p[:,1])  # Find nodes on x=1 boundary
p[ind1[end-r+1],3] = 3  # Set Dirichlet BC type at specific node
p[ind1[1],3] = 1        # Set another BC type at first node

el = Int.(el);
DOF = ID(p)  # Assign global degrees of freedom

# Remove duplicate DOFs
findind = DOF
for i in eachindex(findind)
    tempind = [findind[1:i-1]; NaN; findind[i+1:end]]
    ind3 = findind[i] .== tempind
    if any(ind3)
        ind4 = findall(ind3)
        findind[ind4[end]] = 0
    end
end

Ndof = maximum(DOF)  # Total number of degrees of freedom
gj = intWeights(r,2,2)  # Gaussian quadrature weights

# Compute element matrices
Kel = KMtrx(r,hx,hy,gj)  # Element stiffness matrix
Mel = MMtrx(r,hx,hy,gj)  # Element mass matrix

# Assemble global matrices
println("Assembling matrices.")
K = spzeros(Int(Ndof),Int(Ndof))  # Global stiffness matrix
M = spzeros(Int(Ndof),Int(Ndof))  # Global mass matrix

# Assembly loop over elements
for l = 1:Nel
    for i = 1:(r+1)^2
        m = Int(DOF[Int(el[l,i])])
        if m != 0  # Skip constrained DOFs
            for j = 1:(r+1)^2
                n = Int(DOF[Int(el[l,j])])
                if n != 0
                    K[m,n] = K[m,n] + Kel[i,j]
                    M[m,n] = M[m,n] + Mel[i,j]
                end
            end
        end
    end
end 

println("Matrices Assembled")

# Initialize solution vectors
uel = zeros((r+1)^2, 1)    # Element solution vector
fel = zeros((r+1)^2, 1)    # Element force vector
ufull = zeros(size(p,1), 1)  # Global solution vector
F = zeros(Int(Ndof), 1)      # Global force vector

# Assemble right-hand side force vector
for l = 1:Nel
    # Compute element contributions
    for m = 1:(r+1)^2
        local ind1 = Int(DOF[Int(el[l,m])])
        if ind1 != 0
            uel[m] = 0  # Initialize element solution
            fel[m] = f_ex(osc_f, p[el[l,m],1], p[el[l,m],2])  # Source term
        else
            fel[m] = f_ex(osc_f, p[el[l,m],1], p[el[l,m],2])
        end
    end
    
    # Add element contribution to global force vector
    Fel = Mel * fel  
    for k = 1:(r+1)^2
        local ind1 = Int(DOF[el[l,k]])
        if ind1 != 0
            F[ind1] = F[ind1] + Fel[k]
        end
    end
end

# Solve linear system Ku = F
println("Right hand side assembled")
u = K\F  # Solve using backslash operator
println("Linear system solved")

# Map solution back to full domain including boundary nodes
ufull[findall(x -> x != 0, vec(DOF))] = u

# Compute error metrics
errt, nu, nuh = femerror(p, el, ufull, hx, hy, r, osc_f)
rel_err = errt/nuh
@printf("| Norm Sol.: %3.5g | Norm. App.: %3.5g | Rel. err. %3.5e\n", nu, nuh, rel_err)

# Save mesh points to file
f = open("points.dat", "w")
for i = 1:size(p,1)
    for j = 1:3
        print(f, p[i,j], " ")
    end
    println(f, "")
end
close(f)

# Save solution to file
f = open("sol.dat", "w")
for i = 1:size(ufull,1)
    println(f, ufull[i])
end
close(f)

# Visualize solution using 3D scatter plot
scatter3d(p[:,1], p[:,2], vec(ufull),
    marker_z=vec(ufull),
    markersize=2,
    markerstrokewidth=0,
    xlabel="x",
    ylabel="y", 
    zlabel="u(x,y)",
    title="Solution",
    colorbar=true
)
