using Gridap
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using LineSearches: BackTracking, Static, MoreThuente
using FillArrays
using Gridap.Fields: meas

"""
LidDrivenCavityFlow
"""

#Parameters
function stretching_y_function(x)
  gamma1 = 2.5
  -tanh.(gamma1 .* (x)) ./ tanh.(gamma1)
end


function stretching(x::Point)
  m = zeros(length(x))
  m[1] = stretching_y_function(x[1])

  
  m[2] = stretching_y_function(x[2])
  Point(m)
end


Re = 100
L = 0.5
u0 = 1

ν = u0*L/Re #m2/s 

order = 1 #Order of pressure and velocity

N = 100; #cells per dimensions
hf = VectorValue(0,0)


#ODE settings
t0 = 0.0
dt = 5e-4
tF = 5*dt

Ntimestep = (tF-t0)/dt
θ = 0.5

initial_condition = false #print model of initial condition



#MESH DEFINITION
domain = (-L, L, -L, L)
partition = (N, N)
model = CartesianDiscreteModel(domain, partition,map=stretching)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri1",[5,])
add_tag_from_tags!(labels,"diri0",[1,2,3,5,6,7,8])
add_tag_from_tags!(labels,"p",[4,])



writevtk(model,"model")

#ANALITICAL SOLUTION, used also for initial condition

u_wall(x, t) = VectorValue(0,0)
u_wall(t::Real) = x -> u_wall(x, t)

u_top(x, t) = VectorValue(u0, 0)
u_top(t::Real) = x -> u_top(x, t)

p0(x, t) = 0
p0(t::Real) = x -> p0(x, t)




reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["diri0","diri1"])
reffeₚ = ReferenceFE(lagrangian,Float64, order)
#reffeₚ = ReferenceFE(lagrangian,Float64,order-1; space=:P)
#reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)
#Q = TestFESpace(model,reffeₚ, conformity=:L2, constraint=:zeromean)
#Q = TestFESpace(model,reffeₚ, conformity=:L2, dirichlet_tags="interior")
Q = TestFESpace(model,reffeₚ, conformity=:H1, dirichlet_tags=["p"])

#Since we impose Dirichlet boundary conditions on the entire boundary ∂Ω, the mean value of the pressure is constrained to zero in order have a well posed problem
#Q = TestFESpace(model, reffeₚ)


#Transient is just for the fact that the boundary conditions change with time
U = TransientTrialFESpace(V, [u_wall, u_top])
#P = TrialFESpace(Q) #?transient
P = TransientTrialFESpace(Q, p0) #?transient



Y = MultiFieldFESpace([V, Q]) #?transient
X = TransientMultiFieldFESpace([U, P])

degree = 4*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)


h = lazy_map(h->h^(1/2),get_cell_measure(Ω))



# Momentum residual, without the viscous term
Rm(t,(u,p)) = ∂t(u) + u⋅∇(u) + ∇(p) - hf

# Continuity residual
Rc(u) = ∇⋅u


function τ(u,h)
    r = 2
    c₁= 1
    c₂=16
    β=u0
    τ₂ = h^2/(4*ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
    u = val(norm(u))
    
    if iszero(u)
        return τ₂
        
    end
    τ₃ =  dt/2 #h/(2*u) #0  dt/2 #

    τ₁ = h/(2*β) #h/(2*u) #
    return (1/τ₁^r + c₁^2/τ₂^r + c₂/τ₃^r)^(-1/r)
    
end


#τb(u,h) = (u⋅u)*τ(u,h)
τb(u,h) = (u⋅u)*τ(u,h)

var_equations(t,(u,p),(v,q)) = ∫(
    ν*∇(v)⊙∇(u) # Viscous term
    + v⊙Rm(t,(u,p)) # Other momentum terms
    + q*Rc(u)
 )dΩ # Continuity


stab_equations(t,(u,p),(v,q)) = ∫(  (τ∘(u,h)*(u⋅∇(v) + ∇(q)))⊙Rm(t,(u,p)) # First term: SUPG, second term: PSPG
    +τb∘(u,h)*(∇⋅v)⊙Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ


res(t,(u,p),(v,q)) = var_equations(t,(u,p),(v,q)) + stab_equations(t,(u,p),(v,q))
op = TransientFEOperator(res,X,Y)
nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)

solver = FESolver(nls)




U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

uh0 = interpolate_everywhere(VectorValue(0,0), U0)
ph0 = interpolate_everywhere(0, P0)

xh0 = interpolate_everywhere([uh0, ph0], X0)


tau =τ∘(uh0, h)

writevtk(Ω, "Start_", cellfields=["tau" => tau, "uh0"=>uh0])


ode_solver = ThetaMethod(nls, dt, θ)

sol_t = solve(ode_solver, op, xh0, t0, tF)


_t_nn = t0
iteration = 0
createpvd("TV_2d") do pvd
  for (xh_tn, tn) in sol_t
    global _t_nn
    _t_nn += dt
    global iteration
    iteration += 1
    println("it_num = $iteration\n")
    uh_tn = xh_tn[1]
    ph_tn = xh_tn[2]
    ωh_tn = ∇ × uh_tn
  
      pvd[tn] = createvtk(Ω, "Results/TV_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn])
    
  end

end

