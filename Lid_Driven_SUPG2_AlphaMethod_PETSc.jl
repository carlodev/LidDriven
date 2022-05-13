using Gridap
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using LineSearches: BackTracking, Static, MoreThuente
using FillArrays
using Gridap.Fields: meas
using Statistics
using JLD2
using GridapPETSc


"""
LidDrivenCavityFlow
The domain is a square refined close to the walls
Non slip condition on the 3 walls, constant velocity on the top wall
The pressure is set to zero in the bottom left corner

"""


  #Parameters
  function stretching_y_function(x)
    gamma1 = 2.5
    S = 0.5815356159649889 #for rescaling the function over the domain -0.5 -> 0.5
    -tanh.(gamma1 .* (x)) ./ tanh.(gamma1) .* S
  end


  function stretching(x::Point)
    m = zeros(length(x))
    m[1] = stretching_y_function(x[1])
    m[2] = stretching_y_function(x[2])
    Point(m)
  end


  Re = 1000
  L = 0.5 #1/2 of the domain dimensions
  u0 = 1 #Top velocity

  ν = u0 * 2 * L / Re #m2/s 

  order = 1 #Order of pressure and velocity

  N = 100 #cells per dimensions
  hf = VectorValue(0, 0)


  #ODE settings
  t0 = 0.0
  dt = 0.5
  tF = 5*dt

  Ntimestep = (tF - t0) / dt
  initial_condition = false #print model of initial condition



  #MESH DEFINITION
  domain = (-L, L, -L, L)
  partition = (N, N)
  model = CartesianDiscreteModel(domain, partition, map=stretching)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri1", [5,])
  add_tag_from_tags!(labels, "diri0", [1, 2, 3, 4, 6, 7, 8])
  add_tag_from_tags!(labels, "p", [4,])



  writevtk(model, "model")

  #Function of velocity and pressure on the boundaries
  u_wall(x, t) = VectorValue(0, 0)
  u_wall(t::Real) = x -> u_wall(x, t)

  u_top(x, t) = VectorValue(u0, 0)
  u_top(t::Real) = x -> u_top(x, t)

  p0(x, t) = 0
  p0(t::Real) = x -> p0(x, t)




  reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
  V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["diri0", "diri1"]) #diri0 union of the walls, diri1 top wall
  reffeₚ = ReferenceFE(lagrangian, Float64, order)
  Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["p"])

  U = TransientTrialFESpace(V, [u_wall, u_top])
  P = TransientTrialFESpace(Q, p0)



  Y = MultiFieldFESpace([V, Q])
  X = TransientMultiFieldFESpace([U, P])

  degree = 4 * order
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)


  h = lazy_map(h -> h^(1 / 2), get_cell_measure(Ω)) #get the area of each cell (because in 2D), then ^1/2



  # Momentum residual, without the viscous term
  Rm(t, (u, p)) = ∂t(u) + u ⋅ ∇(u) + ∇(p) - hf

  # Continuity residual
  Rc(u) = ∇ ⋅ u


  function τ(u, h)

    β = u0
    τ₂ = h^2 / (4 * ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
    u = val(norm(u))

    if iszero(u)
      return τ₂

    end
    τ₃ = dt / 2 #h/(2*u) #0  dt/2 #

    τ₁ = h / (2 * u) #h/(2*u) #
    return 1 / (1 / τ₁ + 1 / τ₂ + 1 / τ₃)

  end


  #τb(u,h) = (u⋅u)*τ(u,h)
  τb(u, h) = (u ⋅ u) * τ(u, h)

  var_equations(t, (u, p), (v, q)) = ∫(
    ν * ∇(v) ⊙ ∇(u) # Viscous term
    + v ⊙ Rm(t, (u, p)) # Other momentum terms
    + q * Rc(u)
  )dΩ # Continuity


  stab_equations(t, (u, p), (v, q)) = ∫((τ ∘ (u, h) * (u ⋅ ∇(v) + ∇(q))) ⊙ Rm(t, (u, p)) # First term: SUPG, second term: PSPG
                                        +
                                        τb ∘ (u, h) * (∇ ⋅ v) ⊙ Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
  )dΩ


  res(t, (u, p), (v, q)) = var_equations(t, (u, p), (v, q)) + stab_equations(t, (u, p), (v, q))
  op = TransientFEOperator(res, X, Y)



  nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)

  solver = FESolver(nls)




  U0 = U(0.0)
  P0 = P(0.0)
  X0 = X(0.0)

  #initial condition
  uh0 = interpolate_everywhere(VectorValue(0, 0), U0)
  ph0 = interpolate_everywhere(0, P0)
  xh0 = interpolate_everywhere([uh0, ph0], X0)

  #initial condition - derivative; It is not ideal, the first iteration are not perfect
  vuh0 = interpolate_everywhere(VectorValue(0, 0), U0)
  vph0 = interpolate_everywhere(0, P0)
  vxh0 = interpolate_everywhere([vuh0, vph0], X0)

  #For extracting results, for fft
  N_samples = 100
  U_vector = zeros(length(uh0.free_values), N_samples)

  ρ∞ = 0.8 #ρ∞=1 no dissipation, ρ∞=0 max dissipation, ρ∞=0.5 quite good 
  ode_solver = GeneralizedAlpha(nls, dt, ρ∞)
  sol_t = solve(ode_solver, op, (xh0, vxh0), t0, tF)


  _t_nn = t0
  iteration = 0
  s_iteration = 1

  createpvd("TV_2d") do pvd
    for (xh_tn, tn) in sol_t
      global _t_nn
      _t_nn += dt
      global iteration
      global s_iteration
      global U_vector
      iteration += 1
      println("it_num = $iteration\n")
      uh_tn = xh_tn[1]
      ph_tn = xh_tn[2]
      ωh_tn = ∇ × uh_tn
      if tn > 100
        U_vector[:, s_iteration] = uh_tn.free_values
        s_iteration += 1
      end

      if mod(iteration, 1) < 1 #useful for printing not every iteration
        pvd[tn] = createvtk(Ω, "Results/TV_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn])
      end

    end

  end

  @save "LidDriven.jld2" U_vector
