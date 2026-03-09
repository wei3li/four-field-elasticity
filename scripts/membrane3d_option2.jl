using Gridap, LineSearches
import GridapGmsh: GmshDiscreteModel
import Unicode: graphemes

include("utilities.jl")


λ, μ = [0.01], 1.0

zv = VectorValue(0.0, 0.0, 0.0)
T1(X) = X[1] ≈ 48 ? zv : zv
T2(X) = X[1] ≈ 48 ? VectorValue(λ[1], 0.0, 0.0) : zv
T3(X) = X[1] ≈ 48 ? zv : zv

XA = Point(48.0, 60.0, 5.0)

filename = joinpath(DATA_DIR, "membrane3d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "f,UxA,UyA,UzA,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")


function run_membrane3d_option2(reffe_names, meshname)
  chars = collect(graphemes(reffe_names))
  ords = [parse(Int, c) for c in chars if all(isdigit, c)]
  typs = [c for c in chars if !any(isdigit, c)]

  simpl = 'P' ∈ typs[1]
  change_type(typ) = ('P' ∈ typ || 'Q' ∈ typ) ? (simpl ? "Q" : typ) : typ

  rfu = symbol_to_space[change_type(typs[1])]
  reffe_U = ('P' ∈ typs[1] || 'Q' ∈ typs[1]) ? rfu(VectorValue{3,Float64}, ords[1]) : rfu(ords[1])
  rff = symbol_to_space[change_type(typs[2])]
  reffe_K = ('P' ∈ typs[2] || 'Q' ∈ typs[2]) ? rff(VectorValue{3,Float64}, ords[2]) : rff(ords[2])
  reffe_P = symbol_to_space[change_type(typs[3])](ords[3])
  reffe_p = symbol_to_space[change_type(typs[4])](Float64, ords[4])

  model = GmshDiscreteModel(joinpath(MESH_DIR, "membrane3d/membrane3d_$(meshname).msh"))
  neum_tags = filter(s -> contains(s, "neum"), get_face_labeling(model).tag_to_name)

  confs = [symbol_to_conformity[typs[i]] for i in 1:4]
  L2U = FESpace(model, reffe_U, conformity=confs[1])
  L2K = FESpace(model, reffe_K, conformity=confs[2])
  DP1 = TestFESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=[neum_tags..., "diri_surfaces_z", "diri_points_z", "diri_lines_z"])
  DP2 = TestFESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=[neum_tags..., "diri_surfaces_z", "diri_points_z", "diri_lines_z"])
  DP3 = TestFESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=neum_tags)
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([L2U, L2K, L2K, L2K, DP1, DP2, DP3, L2p])

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)

  function a((U, K1, K2, K3, P1, P2, P3, p), (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    K, P = vectors_to_tensor(K1, K2, K3), vectors_to_tensor(P1, P2, P3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)
    divP = scalars_to_vector(∇ ⋅ P1, ∇ ⋅ P2, ∇ ⋅ P3)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2, ∇ ⋅ ψ3)

    F = K + I3
    J = det(F)
    pQ = p * inv(F)'

    ∫(V ⋅ divP -
      κ ⊙ (μ * F + pQ - P) +
      ψ ⊙ K + divψ ⋅ U -
      q * (log ∘ (J)))dΩ
  end

  function j(
    (U, K1, K2, K3, P1, P2, P3, p),
    (dU, dK1, dK2, dK3, dP1, dP2, dP3, dp),
    (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))

    K = vectors_to_tensor(K1, K2, K3)
    dK, dP = vectors_to_tensor(dK1, dK2, dK3), vectors_to_tensor(dP1, dP2, dP3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)
    divdP = scalars_to_vector(∇ ⋅ dP1, ∇ ⋅ dP2, ∇ ⋅ dP3)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2, ∇ ⋅ ψ3)

    F = K + I3
    iF = inv(F)
    iFdK = iF ⋅ dK

    ∫(V ⊙ divdP -
      κ ⊙ (μ * dK - p * (iFdK ⋅ iF)' + dp * iF' - dP) +
      ψ ⊙ dK + divψ ⋅ dU -
      q * tr(iFdK))dΩ
  end

  r(x, y) = a(x, y)

  Udof, Kdof = num_free_dofs(L2U), 3num_free_dofs(L2K)
  Pdof = num_free_dofs(DP1) + num_free_dofs(DP2) + num_free_dofs(DP3)
  pdof, xdof = num_free_dofs(L2p), num_free_dofs(Y)
  ncells = num_cells(model)
  meshtype = meshname[1:findlast("_", meshname).start-1]
  meshidx = meshname[findlast("_", meshname).start+1:end]
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = cbrt(max_cell_area)
  simpl && (h = cbrt(6max_cell_area / π))

  diri_tags = filter(s -> startswith(s, "diri"), labels.tag_to_name)
  add_tag_from_tags!(labels, "diri_xyz", filter(t -> contains(t, "_xyz"), diri_tags))
  add_tag_from_tags!(labels, "diri_z", filter(t -> contains(t, "_z"), diri_tags))

  function correct_displacement(Kh)
    reffe_U1 = ReferenceFE(lagrangian, VectorValue{3,Float64}, ords[2] + 1)
    H10 = TestFESpace(model, reffe_U1, conformity=:H1,
      dirichlet_tags=["diri_xyz", "diri_z"],
      dirichlet_masks=[(true, true, true), (false, false, true)])
    au(UU, UV) = ∫(∇(UV) ⊙ ∇(UU))dΩ
    lu(UV) = ∫(∇(UV) ⊙ Kh)dΩ
    solve(AffineFEOperator(au, lu, H10, H10))
  end

  function record_norms(λi, xh, dΩ)
    Uh, Kh1, Kh2, Kh3, Ph1, Ph2, Ph3, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2, Kh3), vectors_to_tensor(Ph1, Ph2, Ph3)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2, ∇ ⋅ Ph3)
    UxA, UyA, UzA = Uh(XA)
    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ), compute_l2_norm(Kh, dΩ)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ), compute_l2_norm(ph, dΩ)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ)^2)

    write_line(filename,
      "$reffe_names,$meshtype,$meshidx,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λi,$UxA,$UyA,$UzA,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    Uh1 = correct_displacement(Kh)
    UxA1, UyA1, UzA1 = Uh1(XA)
    Udof1 = length(Uh1.free_values)
    xdof1 = Udof1 + Kdof + Pdof + pdof
    Ul2norm1 = compute_l2_norm(Uh1, dΩ)
    write_line(filename,
      "$(reffe_names)(corr),$meshtype,$meshidx,$ncells,$h,$Udof1,$Kdof,$Pdof,$pdof,$xdof1,",
      "$λi,$UxA1,$UyA1,$UzA1,$Ul2norm1,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = [
      "displacement" => Uh, "corrected displacement" => Uh1,
      "strain" => Kh, "stress" => Ph, "div stress" => ∇Ph,
      "pressure" => ph, "Jacobian (f = $λi)" => det(Kh + I3)]
    writevtk(Ω, joinpath(DATA_DIR, "membrane3d_$(reffe_names)_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang())
  solver = FESolver(nls)
  cache = nothing
  xh = zero(Y)

  println_now("\n[$reffe_names] $ncells cells optimisation begins...")
  for λi in [0:0.05:0.2..., 0.225:0.025:0.35..., 0.36:0.02:0.4...]
    λ[1] = λi
    DPT1, DPT2, DPT3 = TrialFESpace(DP1, T1), TrialFESpace(DP2, T2), TrialFESpace(DP3, T3)
    X = MultiFieldFESpace([L2U, L2K, L2K, L2K, DPT1, DPT2, DPT3, L2p])
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println_now("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)

    λi in (0.2, 0.4) && record_norms(λi, xh, dΩ⁺)
    println("λ = $λi optimisation finished\n..................................\n")
  end

end

pairs = ["P̄0d1d1P1", "P̄1d2d2P2"]
pair, meshidx = pairs[parse(Int, ARGS[1])], ARGS[2]
try
  run_membrane3d_option2(pair, "delaunay_$meshidx")
catch e
  println("[run_membrane3d_option2] pair $pair mesh $meshidx error: $e")
end
