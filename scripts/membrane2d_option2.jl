using Gridap, LineSearches
import Unicode: graphemes

include("utilities.jl")


λ, μ = [0.01], 1.0

Te((x, y)) = x ≈ 48.0 ? VectorValue(λ[1], 0.0) : VectorValue(0.0, 0.0)

XA = Point(48.0, 60.0)

filename = joinpath(DATA_DIR, "membrane2d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "f,UxA,UyA,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")


function run_membrane2d_option2(reffe_names, ncell)
  chars = collect(graphemes(reffe_names))
  ords = [parse(Int, c) for c in chars if all(isdigit, c)]
  typs = [c for c in chars if !any(isdigit, c)]

  simpl = 'P' ∈ typs[1]
  change_type(typ) = ('P' ∈ typ || 'Q' ∈ typ) ? (simpl ? "Q" : typ) : typ
  reffe_U = symbol_to_space[change_type(typs[1])](VectorValue{2,Float64}, ords[1])
  rff = symbol_to_space[change_type(typs[2])]
  reffe_K = ('P' ∈ typs[2] || 'Q' ∈ typs[2]) ? rff(VectorValue{2,Float64}, ords[2]) : rff(ords[2])
  reffe_P = symbol_to_space[change_type(typs[3])](ords[3])
  reffe_p = symbol_to_space[change_type(typs[4])](Float64, ords[4])

  model = _init_membrane_model1x1()
  factors = factorise_integer(ncell)
  for ftr in factors
    model = Gridap.Adaptivity.refine(model, ftr)
  end
  simpl && (model = model |> simplexify)

  confs = [symbol_to_conformity[typ] for typ in typs]
  L2U = FESpace(model, reffe_U, conformity=confs[1])
  L2K = FESpace(model, reffe_K, conformity=confs[2])
  DvP = TestFESpace(model, reffe_P, conformity=confs[3], dirichlet_tags="neumann")
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([L2U, L2K, L2K, DvP, DvP, L2p])

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)

  function a((U, K1, K2, P1, P2, p), (V, κ1, κ2, ψ1, ψ2, q))
    K, P = vectors_to_tensor(K1, K2), vectors_to_tensor(P1, P2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)
    divP = scalars_to_vector(∇ ⋅ P1, ∇ ⋅ P2)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2)

    F = K + I2
    J = det(F)
    pQ = p * inv(F)'

    ∫(V ⋅ divP -
      κ ⊙ (μ * F + pQ - P) +
      ψ ⊙ K + divψ ⋅ U -
      q * (log ∘ (J))
    )dΩ
  end

  function j((U, K1, K2, P1, P2, p), (dU, dK1, dK2, dP1, dP2, dp), (V, κ1, κ2, ψ1, ψ2, q))
    K = vectors_to_tensor(K1, K2)
    dK, dP = vectors_to_tensor(dK1, dK2), vectors_to_tensor(dP1, dP2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)
    divdP = scalars_to_vector(∇ ⋅ dP1, ∇ ⋅ dP2)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2)

    F = K + I2
    iF = inv(F)
    iFdK = iF ⋅ dK

    ∫(V ⊙ divdP -
      κ ⊙ (μ * dK - p * (iFdK ⋅ iF)' + dp * iF' - dP) +
      ψ ⊙ dK + divψ ⋅ dU -
      q * tr(iFdK)
    )dΩ
  end

  r(x, y) = a(x, y)

  Udof, Kdof = num_free_dofs(L2U), 2num_free_dofs(L2K)
  Pdof, pdof = 2num_free_dofs(DvP), num_free_dofs(L2p)
  xdof = num_free_dofs(Y)

  meshtype, meshidx = "cart", ncell
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = sqrt(max_cell_area)
  simpl && (h = sqrt(4 * max_cell_area / sqrt(3)))
  ncells = num_cells(model)

  function correct_displacement(Kh)
    reffe_U1 = ReferenceFE(lagrangian, VectorValue{2,Float64}, ords[2] + 1)
    H10 = FESpace(model, reffe_U1, dirichlet_tags="dirichlet")
    H1U = TrialFESpace(H10, VectorValue(0., 0.))
    au(UU, UV) = ∫(∇(UV) ⊙ ∇(UU))dΩ
    lu(UV) = ∫(∇(UV) ⊙ Kh)dΩ
    solve(AffineFEOperator(au, lu, H1U, H10))
  end

  function write_record(λi, xh)
    Uh, Kh1, Kh2, Ph1, Ph2, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2), vectors_to_tensor(Ph1, Ph2)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2)
    UxA, UyA = Uh(XA)
    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ⁺), compute_l2_norm(Kh, dΩ⁺)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ⁺), compute_l2_norm(ph, dΩ⁺)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ⁺)^2)
    write_line(filename,
      "$reffe_names,$meshtype,$meshidx,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λi,$UxA,$UyA,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    Uh1 = correct_displacement(Kh)
    UxA1, UyA1 = Uh1(XA)
    Ul2norm1 = compute_l2_norm(Uh1, dΩ⁺)
    Udof1 = length(get_free_dof_values(Uh1))
    xdof1 = Udof1 + Kdof + Pdof + pdof
    write_line(filename,
      "$reffe_names(corr),$meshtype,$meshidx,$ncells,$h,$Udof1,$Kdof,$Pdof,$pdof,$xdof1,",
      "$λi,$UxA1,$UyA1,$Ul2norm1,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = [
      "displacement" => Uh, "corrected displacement" => Uh1,
      "strain" => Kh, "stress" => Ph, "div stress" => ∇Ph,
      "pressure" => ph, "Jacobian (f = $λi)" => det(Kh + I2)
    ]
    writevtk(Ω, joinpath(DATA_DIR, "membrane2d_$(reffe_names)_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang(), iterations=500)
  solver = FESolver(nls)
  cache = nothing
  xh = zero(Y)

  println("\n[$reffe_names] $ncells cells, $xdof DoFs optimisation begins...")
  for λi in [0:0.02:0.3..., 0.3125:0.0175:0.4...]
    λ[1] = λi
    X = MultiFieldFESpace([L2U, L2K, L2K, DvP, TrialFESpace(DvP, Te), L2p])
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println("λ = $λi optimisation starting using ln(J)...\n")
    xh, cache = solve!(xh, solver, op, cache)
    λi in 0.2:0.1:0.4 && write_record(λi, xh)
    println_now("λ = $λi optimisation finished\n..................................\n")
  end
end

pairs = ["P̄0d1d1P1", "P̄1d2d2P2"]
try
  run_membrane2d_option2(pairs[parse(Int, ARGS[1])], parse(Int, ARGS[2]))
catch e
  println("[run_membrane2d_option2] error: $e")
end
