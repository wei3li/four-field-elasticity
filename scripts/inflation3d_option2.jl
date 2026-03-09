using Gridap, GridapGmsh, LineSearches
import Unicode: graphemes
import BSON

include("utilities.jl")


λ, μ = [1.1], 1.0
Rin, Rout = 0.5, 1.0

r(R) = cbrt(R^3 + (λ[1]^3 - 1) * Rin^3)
g(R) = R * (3 * (r(R)^3) + (λ[1]^3 - 1) * Rin^3) / (r(R)^4)

function Ue(X)
  R = sqrt(X ⋅ X)
  (r(R) / R - 1) * X
end

Ke(X) = ∇(Ue)(X)

function pe(X)
  R = sqrt(X ⋅ X)
  -μ * (Rout / r(Rout))^4 + μ / 2 * (g(R) - g(Rout))
end

_F(U) = ∇(U) + I3
_F(U, X) = ∇(U)(X) + I3

function Pe(X)
  F = _F(Ue, X)
  J = det(F)
  μ * F + pe(X) * J * inv(F)'
end

Ue_sym(X) = VectorValue(0.0, 0.0, 0.0)

function gen_diri_mask_and_values(tag)
  mask = [false for _ in 1:3]
  for (i, s) in enumerate(["x", "y", "z"])
    contains(tag, s) && (mask[i] = true)
  end
  vals = sum(mask) == 3 ? Ue : Ue_sym
  Tuple(mask), vals
end

get_col(S, i) = VectorValue(S[1, i], S[2, i], S[3, i])
col1(S) = get_col(S, 1)
col2(S) = get_col(S, 2)
col3(S) = get_col(S, 3)


filename = joinpath(DATA_DIR, "inflation3d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "lambda,fem_Ul2err,fem_Kl2err,fem_Pl2err,fem_Phdiverr,fem_pl2err")

function run_inflation3d_option2(reffe_names, meshname)
  chars = collect(graphemes(reffe_names))
  ords = [parse(Int, c) for c in chars if all(isdigit, c)]
  typs = [c for c in chars if !any(isdigit, c)]

  simpl = 'P' ∈ typs[1]
  change_type(typ) = ('P' ∈ typ || 'Q' ∈ typ) ? (simpl ? "Q" : typ) : typ
  reffe_U = symbol_to_space[change_type(typs[1])](VectorValue{3,Float64}, ords[1])
  rff = symbol_to_space[change_type(typs[2])]
  reffe_K = ('P' ∈ typs[2] || 'Q' ∈ typs[2]) ? rff(VectorValue{3,Float64}, ords[2]) : rff(ords[2])
  reffe_P = symbol_to_space[change_type(typs[3])](ords[3])
  reffe_p = symbol_to_space[change_type(typs[4])](Float64, ords[4])

  model = GmshDiscreteModel(joinpath(MESH_DIR, "inflation3d/inflation3d_$(meshname).msh"))

  confs = [symbol_to_conformity[typs[i]] for i in 1:4]
  L2U = FESpace(model, reffe_U, conformity=confs[1])
  L2K = FESpace(model, reffe_K, conformity=confs[2])
  DP1 = FESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=["neum_surfaces", "diri_surfaces_y", "diri_surfaces_z"])
  DP2 = FESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=["neum_surfaces", "diri_surfaces_x", "diri_surfaces_z"])
  DP3 = FESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=["neum_surfaces", "diri_surfaces_x", "diri_surfaces_y"])
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([L2U, L2K, L2K, L2K, DP1, DP2, DP3, L2p])
  X = Y

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)
  Γ = BoundaryTriangulation(Ω, tags=["diri_points_xyz", "diri_curves_xyz", "diri_surfaces_xyz"])
  dΓ = Measure(Γ, qdeg)
  N = get_normal_vector(Γ)

  function a((U, K1, K2, K3, P1, P2, P3, p), (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    K, P = vectors_to_tensor(K1, K2, K3), vectors_to_tensor(P1, P2, P3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)
    divP = scalars_to_vector(∇ ⋅ P1, ∇ ⋅ P2, ∇ ⋅ P3)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2, ∇ ⋅ ψ3)

    F = K + I3
    J, iF = det(F), inv(F)
    pQ = p * J * iF'

    ∫(V ⋅ divP -
      κ ⊙ (μ * F + pQ - P) +
      ψ ⊙ K + divψ ⋅ U -
      q * (J - 1))dΩ
  end

  function l((V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    ψ = vectors_to_tensor(ψ1, ψ2, ψ3)
    ∫(N ⋅ ψ ⋅ Ue)dΓ
  end

  function j(
    (U, K1, K2, K3, P1, P2, P3, p),
    (dU, dK1, dK2, dK3, dP1, dP2, dP3, dp),
    (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))

    K = vectors_to_tensor(K1, K2, K3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)
    dK = vectors_to_tensor(dK1, dK2, dK3)
    dP = vectors_to_tensor(dP1, dP2, dP3)
    divdP = scalars_to_vector(∇ ⋅ dP1, ∇ ⋅ dP2, ∇ ⋅ dP3)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2, ∇ ⋅ ψ3)

    F = K + I3
    J, iF = det(F), inv(F)
    iFdK = iF ⋅ dK
    JiFdK = J * tr(iFdK)

    ∫(V ⊙ divdP -
      κ ⊙ (μ * dK + p * (JiFdK * iF' - J * (iFdK ⋅ iF)') + dp * J * iF' - dP) +
      ψ ⊙ dK + divψ ⋅ dU -
      q * JiFdK)dΩ
  end

  Udof, Kdof = num_free_dofs(L2U), 3num_free_dofs(L2K)
  Pdof = num_free_dofs(DP1) + num_free_dofs(DP2) + num_free_dofs(DP3)
  pdof, xdof = num_free_dofs(L2p), num_free_dofs(Y)
  ncells = num_cells(model)
  meshtype = meshname[1:findlast("_", meshname).start-1]
  meshidx = meshname[findlast("_", meshname).start+1:end]
  volumes = sort!(collect((∫(1)Measure(Ω, 2))[Ω]), rev=true)
  maxidx = max(length(volumes) ÷ 50, 1)
  h = cbrt(6volumes[maxidx] / π)

  diri_tags = filter(s -> startswith(s, "diri"), model.face_labeling.tag_to_name)
  diri_masks = [gen_diri_mask_and_values(tag)[1] for tag in diri_tags]
  diri_values = [gen_diri_mask_and_values(tag)[2] for tag in diri_tags]

  function correct_displacement(Kh)
    reffe_U1 = ReferenceFE(lagrangian, VectorValue{3,Float64}, ords[2] + 1)
    H10 = TestFESpace(model, reffe_U1, dirichlet_tags=diri_tags, dirichlet_masks=diri_masks)
    H1U = TrialFESpace(H10, diri_values)
    au(UU, UV) = ∫(∇(UV) ⊙ ∇(UU))dΩ
    lu(UV) = ∫(∇(UV) ⊙ Kh)dΩ
    solve(AffineFEOperator(au, lu, H1U, H10))
  end

  function write_record(λi, xh, dΩ)
    Uh, Kh1, Kh2, Kh3, Ph1, Ph2, Ph3, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2, Kh3), vectors_to_tensor(Ph1, Ph2, Ph3)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2, ∇ ⋅ Ph3)

    Ul2err, Kl2err = compute_l2_norm(Uh - Ue, dΩ), compute_l2_norm(Kh - Ke, dΩ)
    Pl2err, pl2err = compute_l2_norm(Ph - Pe, dΩ), compute_l2_norm(ph - pe, dΩ)
    Phdiverr = sqrt(Pl2err^2 + compute_l2_norm(∇Ph - ∇ ⋅ Pe, dΩ)^2)

    write_line(filename,
      "$reffe_names,$meshtype,$meshidx,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λi,$Ul2err,$Kl2err,$Pl2err,$Phdiverr,$pl2err")

    Uh1 = correct_displacement(Kh)
    Ul2err = compute_l2_norm(Uh1 - Ue, dΩ)
    Udof1 = length(Uh1.free_values)
    xdof1 = Udof1 + Kdof + Pdof + pdof
    write_line(filename,
      "$(reffe_names)(corr),$meshtype,$meshidx,$ncells,$h,$h̄,$Udof1,$Kdof,$Pdof,$pdof,$xdof1,",
      "$λi,$Ul2err,NaN,NaN,NaN,NaN")

    fields = [
      "displacement" => Uh, "true displacement" => Ue, "displacement error" => Uh - Ue,
      "corrected displacement" => Uh1, "corrected displacement error" => Uh1 - Ue,
      "strain" => Kh, "true strain" => Ke, "strain error" => Kh - Ke,
      "stress" => Ph, "true stress" => Pe, "stress error" => Ph - Pe,
      "pressure" => ph, "true pressure" => pe, "pressure error" => abs(ph - pe)
    ]
    writevtk(Ω, joinpath(DATA_DIR, "inflation3d_$(reffe_names)_results"), cellfields=fields)
  end

  r(x, y) = a(x, y) - l(y)

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang())
  solver = FESolver(nls)
  cache = nothing
  xh = zero(X)

  println("\n[$reffe_names, $ncells cells] optimisation begins...")
  for λi in 1:0.5:3.0
    λ[1] = λi
    TDP1 = TrialFESpace(DP1, col1 ∘ Pe)
    TDP2 = TrialFESpace(DP2, col2 ∘ Pe)
    TDP3 = TrialFESpace(DP3, col3 ∘ Pe)
    X = MultiFieldFESpace([L2U, L2K, L2K, L2K, TDP1, TDP2, TDP3, L2p])
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println_now("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)
    println("λ = $λi optimisation finished\n..................................\n")
    λi in (2, 3) && write_record(λi, xh, dΩ⁺)
    GC.gc()
  end
end

pairs = ["P̄0d1d1P1", "P̄1d2d2P2"]

pair = pairs[parse(Int, ARGS[1])]
mshidx = ARGS[2]
run_inflation3d_option2(pair, "delaunay_$mshidx")
