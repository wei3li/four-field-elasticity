using Gridap, LineSearches
import Unicode: graphemes

include("utilities.jl")


λ, μ = [1.1], 1.0
Rin, Rout = 0.5, 1.0

r(R) = sqrt(R^2 + (λ[1]^2 - 1) * Rout^2)

function Ue(X)
  R = sqrt(X ⋅ X)
  (r(R) / R - 1) * X
end

Ke(X) = ∇(Ue)(X)

function pe(X)
  R = sqrt(X ⋅ X)
  -μ * R^2 / (r(R))^2 +
  μ * (λ[1]^2 - 1) * Rout^2 / 2 * (1 / r(Rin)^2 - 1 / r(R)^2) +
  μ * log(r(Rin) * R / (Rin * r(R)))
end

_F(U) = ∇(U) + I2
_F(U, X) = ∇(U)(X) + I2

function Pe(X)
  F = _F(Ue, X)
  J = det(F)
  μ * F + pe(X) * J * inv(F)'
end

Ue_sym(X) = VectorValue(0.0, 0.0)

function transform(x::Point)
  r, θ = x[1], x[2]
  l = r + Rin
  Point(l * cospi(θ), l * sinpi(θ))
end

get_col(S, i) = VectorValue(S[1, i], S[2, i])
col1(S) = get_col(S, 1)
col2(S) = get_col(S, 2)

filename = joinpath(DATA_DIR, "inflation2d_record.csv")
create_file_with_header(filename,
  "comb,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "lambda,fem_Ul2err,fem_Kl2err,fem_Pl2err,fem_Phdiverr,fem_pl2err")


function run_inflation2d_option2(reffe_names, ncell)
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

  d = Rout - Rin
  model = CartesianDiscreteModel((0, d, 0, 0.5), (ncell, 2ncell), map=transform)
  simpl && (model = model |> simplexify)

  confs = [symbol_to_conformity[typ] for typ in typs]
  L2U = FESpace(model, reffe_U, conformity=confs[1])
  DvK = FESpace(model, reffe_K, conformity=confs[2])
  DP1 = TestFESpace(model, reffe_P, conformity=confs[3], dirichlet_tags=[1, 3, 5, 7])
  DP2 = TestFESpace(model, reffe_P, conformity=confs[3], dirichlet_tags=[1, 3, 6, 7])
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([L2U, DvK, DvK, DP1, DP2, L2p])
  X = Y

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)
  Γ = BoundaryTriangulation(Ω, tags=[2, 4, 8])
  dΓ = Measure(Γ, qdeg)
  N = get_normal_vector(Γ)

  function a((U, K1, K2, P1, P2, p), (V, κ1, κ2, ψ1, ψ2, q))
    K, P = vectors_to_tensor(K1, K2), vectors_to_tensor(P1, P2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)
    divP = scalars_to_vector(∇ ⋅ P1, ∇ ⋅ P2)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2)

    F = K + I2
    J = det(F)
    pQ = p * J * inv(F)'

    ∫(V ⋅ divP -
      κ ⊙ (μ * F + pQ - P) +
      ψ ⊙ K + divψ ⋅ U -
      q * (J - 1))dΩ
  end

  function l((V, κ1, κ2, ψ1, ψ2, q))
    ψ = vectors_to_tensor(ψ1, ψ2)
    ∫(N ⋅ ψ ⋅ Ue)dΓ
  end

  function j((U, K1, K2, P1, P2, p), (dU, dK1, dK2, dP1, dP2, dp), (V, κ1, κ2, ψ1, ψ2, q))
    K = vectors_to_tensor(K1, K2)
    dK, dP = vectors_to_tensor(dK1, dK2), vectors_to_tensor(dP1, dP2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)
    divdP = scalars_to_vector(∇ ⋅ dP1, ∇ ⋅ dP2)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2)

    F = K + I2
    iF = inv(F)
    J = det(F)
    iFdK = iF ⋅ dK
    JiFdK = J * tr(iFdK)

    ∫(V ⊙ divdP -
      κ ⊙ (μ * dK + p * (JiFdK * iF' - J * (iFdK ⋅ iF)') + dp * J * iF' - dP) +
      ψ ⊙ dK + divψ ⋅ dU -
      q * JiFdK)dΩ
  end

  r(x, y) = a(x, y) - l(y)

  Udof, Kdof = num_free_dofs(L2U), 2num_free_dofs(DvK)
  Pdof, pdof = num_free_dofs(DP1) + num_free_dofs(DP2), num_free_dofs(L2p)
  xdof = num_free_dofs(Y)
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = sqrt(max_cell_area)
  simpl && (h = sqrt(4 * max_cell_area / sqrt(3)))

  function correct_displacement(Kh)
    reffe_U1 = ReferenceFE(lagrangian, VectorValue{2,Float64}, ords[2] + 1)
    t2m = Dict(
      2 => (true, true), 4 => (true, true), 8 => (true, true),
      1 => (false, true), 5 => (false, true),  # y displacement is zero
      3 => (true, false), 6 => (true, false))  # x displacement is zero

    H10 = TestFESpace(model, reffe_U1,
      dirichlet_tags=[2, 4, 8, 1, 3, 5, 6],
      dirichlet_masks=[t2m[2], t2m[4], t2m[8], t2m[1], t2m[3], t2m[5], t2m[6]])
    H1U = TrialFESpace(H10, [Ue, Ue, Ue, Ue_sym, Ue_sym, Ue_sym, Ue_sym])
    au(UU, UV) = ∫(∇(UV) ⊙ ∇(UU))dΩ
    lu(UV) = ∫(∇(UV) ⊙ Kh)dΩ
    solve(AffineFEOperator(au, lu, H1U, H10))
  end

  function record_results(λi, xh, dΩ)
    Uh, Kh1, Kh2, Ph1, Ph2, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2), vectors_to_tensor(Ph1, Ph2)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2)

    fem_Ul2err, fem_Kl2err = compute_l2_norm(Uh - Ue, dΩ), compute_l2_norm(Kh - Ke, dΩ)
    fem_Pl2err, fem_pl2err = compute_l2_norm(Ph - Pe, dΩ), compute_l2_norm(ph - pe, dΩ)
    fem_Phdiverr = sqrt(fem_Pl2err^2 + compute_l2_norm(∇Ph - ∇ ⋅ Pe, dΩ)^2)

    write_line(filename,
      "$reffe_names,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λi,$fem_Ul2err,$fem_Kl2err,$fem_Pl2err,$fem_Phdiverr,$fem_pl2err")

    Uh1 = correct_displacement(Kh)
    fem_Ul2err = compute_l2_norm(Uh1 - Ue, dΩ)
    Udof1 = length(Uh1.free_values)
    xdof1 = Udof1 + Kdof + Pdof + pdof
    write_line(filename,
      "$(reffe_names)(corr),$ncells,$h,$Udof1,$Kdof,$Pdof,$pdof,$xdof1,",
      "$λi,$fem_Ul2err,$fem_Kl2err,$fem_Pl2err,$fem_Phdiverr,$fem_pl2err")

    fields = [
      "displacement" => Uh, "true displacement" => Ue, "displacement error" => Uh - Ue,
      "corrected displacement" => Uh1, "corrected displacement error" => Uh1 - Ue,
      "strain" => Kh, "true strain" => Ke, "strain error" => Kh - Ke,
      "stress" => Ph, "true stress" => Pe, "stress error" => Ph - Pe,
      "pressure" => ph, "true pressure" => pe, "pressure error" => abs(ph - pe),
      "Jacobian (u = $λi)" => det(Kh + I2)
    ]
    writevtk(Ω, joinpath(DATA_DIR, "inflation2d_$(reffe_names)_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang())
  solver = FESolver(nls)
  cache = nothing
  xh = zero(X)

  ncells = num_cells(model)
  println("\n[$reffe_names] $ncells cells optimisation begins...")
  for λi in 1.0:0.5:3.0
    λ[1] = λi
    TDP1, TDP2 = TrialFESpace(DP1, col1 ∘ Pe), TrialFESpace(DP2, col2 ∘ Pe)
    X = MultiFieldFESpace([L2U, DvK, DvK, TDP1, TDP2, L2p])
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)
    λi in (2.0, 3.0) && record_results(λi, xh, dΩ⁺)
    println("λ = $λi optimisation finished\n..................................\n")
  end
end

pairs = ["P̄0d1d1P1", "P̄1d2d2P2"]
ncells = [6, 9, 14, 20, 30, 45]
for ncell in ncells
  run_inflation2d_option2(pairs[parse(Int, ARGS[1])], ncell)
end
