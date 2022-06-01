using Gridap
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using LineSearches: BackTracking, Static, MoreThuente
using FillArrays
using Gridap.Fields: meas
using JLD2
"""
LidDrivenCavityFlow
VMS implementation, based on Bazilevs et al. DOI: 10.1016/j.cma.2007.07.016
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


Re = 100000
L = 0.5
u0 = 1

ν = 2*u0*L/Re #m2/s 

order = 1 #Order of pressure and velocity

N = 100; #cells per dimensions
hf = VectorValue(0,0)


#ODE settings
t0 = 0.0
dt = 0.1
tF = 100*dt

Ntimestep = (tF-t0)/dt
ρ∞=0.8

initial_condition = false #print model of initial condition



#MESH DEFINITION
domain = (-L, L, -L, L)
partition = (N, N)
model = CartesianDiscreteModel(domain, partition,map=stretching)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri1",[5,])
add_tag_from_tags!(labels,"diri0",[1,2,4,3,6,7,8])
add_tag_from_tags!(labels,"p",[4,])



writevtk(model,"model")

#BOUNDARY CONDITIONS FUNCTIONS

u_wall(x, t) = VectorValue(0,0)
u_wall(t::Real) = x -> u_wall(x, t)

u_top(x, t) = VectorValue(u0, 0)
u_top(t::Real) = x -> u_top(x, t)

p0(x, t) = 0
p0(t::Real) = x -> p0(x, t)




reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["diri0","diri1"])
reffeₚ = ReferenceFE(lagrangian,Float64, order)
Q = TestFESpace(model,reffeₚ, conformity=:H1, dirichlet_tags=["p"])


U = TransientTrialFESpace(V, [u_wall, u_top])
P = TransientTrialFESpace(Q, p0) 



Y = MultiFieldFESpace([V, Q]) 
X = TransientMultiFieldFESpace([U, P])

degree = 4*order
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)


#Integration Ω
Ωint(f) = ∫(f)dΩ


# Momentum residual, without viscous term
Rm(t,(u,p)) = ∂t(u) + u⋅∇(u) + ∇(p) #- ν*Δ(u) -hf

# Continuity residual
Rc(u) = ∇⋅u


Bᴳ(t,(u,p),(v,q)) = Ωint(∂t(u)⋅v) + Ωint((u⋅∇(u))⋅v) - Ωint((∇⋅v)*p) + Ωint(q*(∇⋅u)) + ν*Ωint(∇(v)⊙∇(u)) #Variational equations

B_SUPG(t,(u,p),(v,q)) =  Ωint((u⋅∇(v) + ∇(q))⋅(τm∘(u,G,GG).*Rm(t,(u,p)))) + Ωint((∇⋅v)⋅(τc∘(u,gg,G,GG) *Rc(u))) #SUPG terms 

B_VMS1(t,(u,p),(v,q)) = Ωint((u⋅∇(v)')⊙(τm∘(u,G,GG)*Rm(t,(u,p)))) #first VMS term

TRm(t,(u,p)) = τm∘(u,G,GG)*Rm(t,(u,p)) #Function used in the second VMS term

B_VMS2(t,(u,p),(v,q)) = -1*Ωint(∇(v)⊙(outer(TRm(t,(u,p)),TRm(t,(u,p)))))# second VMS term (To be added, it does not work at the moment)



Bᴹ(t,(u,p),(v,q)) = Bᴳ(t,(u,p),(v,q)) + B_SUPG(t,(u,p),(v,q))# + B_VMS1(t,(u,p),(v,q)) + B_VMS2(t,(u,p),(v,q))


h = lazy_map(h->h^(1/2),get_cell_measure(Ω))

ξₖ = get_cell_map(Ω)

Jt     = lazy_map(Broadcasting(∇),ξₖ)
inv_Jt = lazy_map(Operation(inv),Jt)

d = evaluate(inv_Jt, [Point(0.0,0.0)]) #x2

#G = inv_Jt .⋅ inv_Jt
#evaluate(GG, [Point(0.0,0.0)])
#gg = inv_Jt .⋅ inv_Jtt

d ⋅ d

G = (d .⋅ d')[1,:]

GG = G .⊙ G


gg = zeros(num_cells(Ω))
for i  = 1:1:num_cells(Ω)
g = (d[i][1] +  d[i][2])^2+(d[i][3] +  d[i][4])^2
gg[i] = (d[i][1] +  d[i][2])^2 + (d[i][3] +  d[i][4])^2
end

gg
"""
τ₁ = (2/dt).^2 
τ₃ = (ν^2 *GG)
1 ./τ₃
tm = (τ₁ .+ τ₃).^(-1/2)
1 ./(tm.*gg)
"""



function τm(uu,G,GG)
    Cᵢ = 1
    τ₁ = Cᵢ * (2/dt)^2 
    τ₃ = (ν^2 *GG)
    val(x) = x
    function val(x::Gridap.Fields.ForwardDiff.Dual)
        x.value
    end
    uu1 = val(uu[1])
    uu2 = val(uu[2])
    uu_new = VectorValue(uu1,uu2)
    

    if iszero(norm(uu_new))
        return (τ₁ .+ τ₃).^(-1/2)      
    end

    τ₂ = uu_new⋅G⋅uu_new
    #print("$τ₂\n")
    return (τ₁ .+  τ₂ .+ τ₃).^(-1/2)     
end



function τc(uu,gg,G,GG)
   return 1/(τm(uu,G,GG)⋅ gg)
end




res(t,(u,p),(v,q)) = Bᴹ(t,(u,p),(v,q))
op = TransientFEOperator(res,X,Y)
nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)
solver = FESolver(nls)



U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

uh0 = interpolate_everywhere(VectorValue(0,0), U0)
ph0 = interpolate_everywhere(0, P0)
xh0 = interpolate_everywhere([uh0, ph0], X0)


vuh0 = interpolate_everywhere(VectorValue(0,0), U0)
vph0 = interpolate_everywhere(0, P0)
vxh0 = interpolate_everywhere([vuh0, vph0], X0)

#testing the stabilization parameters
taum = τm∘(uh0,G,GG)
tauc = τc∘(uh0,gg,G,GG)
uh1 = interpolate_everywhere(VectorValue(u0,u0), U0)
taum1 = τm∘(uh1,G,GG)
writevtk(Ω, "Start_", cellfields=["taum" => taum, "taum1" => taum1,   "tauc" => tauc, "uh0"=>uh0])



ode_solver = GeneralizedAlpha(nls,dt,ρ∞)

sol_t = solve(ode_solver,op,(xh0,vxh0),t0,tF)


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



