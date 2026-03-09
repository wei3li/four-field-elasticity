using Gridap, LineSearches
import Unicode: graphemes

include("utilities.jl")


λ, μ = [0.01], 1.0

T̄((x, y)) = x ≈ 48.0 ? VectorValue(0.0, λ[1]) : VectorValue(0.0, 0.0)

XA = Point(48.0, 60.0)

filename = joinpath(DATA_DIR, "membrane2d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "f,UxA,UyA,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")

function run_membrane2d_option0(reffe_names, ncell)
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

  model = init_membrane_model1x1()
  factors = factorise_integer(ncell)
  for ftr in factors
    model = Gridap.Adaptivity.refine(model, ftr)
  end
  simpl && (model = model |> simplexify)

  confs = [symbol_to_conformity[typ] for typ in typs]
  H1U = FESpace(model, reffe_U, conformity=confs[1], dirichlet_tags="dirichlet")
  CuK = FESpace(model, reffe_K, conformity=confs[2])
  DvP = TestFESpace(model, reffe_P, conformity=confs[3],)
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([H1U, CuK, CuK, DvP, DvP, L2p])
  X = Y

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)
  Γ = BoundaryTriangulation(Ω, tags="neumann")
  dΓ = Measure(Γ, qdeg)

  function l((V, κ1, κ2, ψ1, ψ2, q))
    ∫(T̄ ⋅ V)dΓ
  end

  function a((U, K1, K2, P1, P2, p), (V, κ1, κ2, ψ1, ψ2, q))
    K, P = vectors_to_tensor(K1, K2), vectors_to_tensor(P1, P2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)

    F = K + I2
    J = det(F)
    pQ = p * inv(F)'

    ∫(∇(V) ⊙ P +
      κ ⊙ (μ * F + pQ - P) +
      ψ ⊙ (∇(U) - K) +
      q * (log ∘ (J)))dΩ
  end

  function j((U, K1, K2, P1, P2, p), (dU, dK1, dK2, dP1, dP2, dp), (V, κ1, κ2, ψ1, ψ2, q))
    K = vectors_to_tensor(K1, K2)
    dK, dP = vectors_to_tensor(dK1, dK2), vectors_to_tensor(dP1, dP2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)

    F = K + I2
    iF = inv(F)
    iFdK = iF ⋅ dK

    ∫(∇(V) ⊙ dP +
      κ ⊙ (μ * dK - p * (iFdK ⋅ iF)' + dp * iF' - dP) +
      ψ ⊙ (∇(dU) - dK) +
      q * tr(iFdK))dΩ
  end

  r(x, y) = a(x, y) - l(y)

  Udof, Kdof = num_free_dofs(H1U), 2num_free_dofs(CuK)
  Pdof, pdof = 2num_free_dofs(DvP), num_free_dofs(L2p)
  xdof = num_free_dofs(Y)
  meshtype, meshidx = "cart", ncell
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = sqrt(max_cell_area)
  simpl && (h = sqrt(4 * max_cell_area / sqrt(3)))
  ncells = num_cells(model)

  function write_record(λi, xh, dΩ)
    Uh, Kh1, Kh2, Ph1, Ph2, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2), vectors_to_tensor(Ph1, Ph2)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2)
    UxA, UyA = Uh(XA)
    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ), compute_l2_norm(Kh, dΩ)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ), compute_l2_norm(ph, dΩ)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ)^2)

    write_line(filename,
      "$reffe_names,$meshtype,$meshidx,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λi,$UxA,$UyA,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = [
      "displacement" => Uh, "strain" => Kh, "stress" => Ph, "div stress" => ∇Ph,
      "pressure" => ph, "Jacobian (f = $λi)" => det(Kh + I2)
    ]
    writevtk(Ω, joinpath(DATA_DIR, "membrane2d_$(reffe_names)_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang(), iterations=500)
  solver = FESolver(nls)
  cache = nothing
  xh = zero(X)

  println("\n[$reffe_names] $ncells cells, $xdof DoFs optimisation begins...")
  for λi in [0:0.02:0.3..., 0.3125:0.0175:0.4...]
    λ[1] = λi
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println_now("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)
    λi in 0.2:0.1:0.4 && write_record(λi, xh, dΩ⁺)
    println("λ = $λi optimisation finished\n..................................\n")
  end

end

pairs = ["P1c1d̄0P̄0", "P2c2d̄1P̄1"]
try
  run_membrane2d_option0(pairs[parse(Int, ARGS[1])], parse(Int, ARGS[2]))
catch e
  println("[run_membrane2d_option0] error: $e")
end
