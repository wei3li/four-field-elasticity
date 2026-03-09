using Gridap, LineSearches
import GridapGmsh: GmshDiscreteModel
import Unicode: graphemes

include("utilities.jl")


λ, μ = [0.1], 1.0
lc = [1]

Ū((x, y, z)) = VectorValue(0.0, 0.0, z > 0.498 ? λ[1] : 0.0)
Ũ(_) = VectorValue(0.0, 0.0, 0.0)

filename = joinpath(DATA_DIR, "stretch3d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "u,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")


function run_stretch3d_option2(reffe_names, meshname)
  chars = collect(graphemes(reffe_names))
  ords = [parse(Int, c) for c in chars if all(isdigit, c)]
  typs = [c for c in chars if !any(isdigit, c)]

  simpl = true  # 'P' ∈ typs[1]
  change_type(typ) = ('P' ∈ typ || 'Q' ∈ typ) ? (simpl ? "Q" : typ) : typ

  rfu = symbol_to_space[change_type(typs[1])]
  reffe_U = ('P' ∈ typs[1] || 'Q' ∈ typs[1]) ? rfu(VectorValue{3,Float64}, ords[1]) : rfu(ords[1])
  rff = symbol_to_space[change_type(typs[2])]
  reffe_K = ('P' ∈ typs[2] || 'Q' ∈ typs[2]) ? rff(VectorValue{3,Float64}, ords[2]) : rff(ords[2])
  reffe_P = symbol_to_space[change_type(typs[3])](ords[3])
  reffe_p = symbol_to_space[change_type(typs[4])](Float64, ords[4])

  model = GmshDiscreteModel(joinpath(MESH_DIR, "stretch3d/stretch3d_$(meshname).msh"))
  labels = get_face_labeling(model)
  tag_to_name = labels.tag_to_name
  add_tag_from_tags!(labels, "diri_x", filter(s -> endswith(s, "_x"), tag_to_name))
  add_tag_from_tags!(labels, "diri_y", filter(s -> endswith(s, "_y"), tag_to_name))
  add_tag_from_tags!(labels, "diri_z", filter(s -> endswith(s, "_z"), tag_to_name))
  add_tag_from_tags!(labels, "diri_xy", filter(s -> endswith(s, "_xy"), tag_to_name))
  add_tag_from_tags!(labels, "diri_xz", filter(s -> endswith(s, "_xz"), tag_to_name))
  add_tag_from_tags!(labels, "diri_yz", filter(s -> endswith(s, "_yz"), tag_to_name))
  add_tag_from_tags!(labels, "diri_xyz", filter(s -> endswith(s, "_xyz"), tag_to_name))

  neum_tags = filter(s -> startswith(s, "neum"), tag_to_name)
  confs = [symbol_to_conformity[typs[i]] for i in 1:4]
  L2U = FESpace(model, reffe_U, conformity=confs[1])
  L2K = FESpace(model, reffe_K, conformity=confs[2])
  DP1 = TestFESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=[neum_tags..., "diri_y", "diri_z", "diri_yz"])
  DP2 = TestFESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=[neum_tags..., "diri_x", "diri_z", "diri_xz"])
  DP3 = TestFESpace(model, reffe_P, conformity=confs[3],
    dirichlet_tags=[neum_tags..., "diri_x", "diri_y", "diri_xy"])
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([L2U, L2K, L2K, L2K, DP1, DP2, DP3, L2p])
  X = Y

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)

  Γ = BoundaryTriangulation(Ω, tags="diri_xyz")
  dΓ = Measure(Γ, qdeg)
  N = get_normal_vector(Γ)

  function a((U, K1, K2, K3, P1, P2, P3, p), (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    K, P = vectors_to_tensor(K1, K2, K3), vectors_to_tensor(P1, P2, P3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)
    divP = scalars_to_vector(∇ ⋅ P1, ∇ ⋅ P2, ∇ ⋅ P3)
    divψ = scalars_to_vector(∇ ⋅ ψ1, ∇ ⋅ ψ2, ∇ ⋅ ψ3)

    F = K + I3
    J = det(F)

    if lc[1] == 1
      pQ = p * inv(F)'

      ∫(V ⋅ divP -
        κ ⊙ (μ * F + pQ - P) +
        ψ ⊙ K + divψ ⋅ U -
        q * (log ∘ (J))
      )dΩ
    else
      pQ = p * J * inv(F)'

      ∫(V ⋅ divP -
        κ ⊙ (μ * F + pQ - P) +
        ψ ⊙ K + divψ ⋅ U -
        q * (J - 1))dΩ
    end
  end

  function l((V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    ψ = vectors_to_tensor(ψ1, ψ2, ψ3)
    ∫(N ⋅ ψ ⋅ Ū)dΓ
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
    J, iF = det(F), inv(F)
    iFdK = iF ⋅ dK

    if lc[1] == 1
      ∫(V ⊙ divdP -
        κ ⊙ (μ * dK - p * (iFdK ⋅ iF)' + dp * iF' - dP) +
        ψ ⊙ dK + divψ ⋅ dU -
        q * tr(iFdK)
      )dΩ
    else
      JiFdK = J * tr(iFdK)

      ∫(V ⊙ divdP -
        κ ⊙ (μ * dK + p * (JiFdK * iF' - J * (iFdK ⋅ iF)') + dp * J * iF' - dP) +
        ψ ⊙ dK + divψ ⋅ dU -
        q * JiFdK)dΩ
    end
  end

  r(x, y) = a(x, y) - l(y)

  Udof, Kdof = num_free_dofs(L2U), 3num_free_dofs(L2K)
  Pdof = num_free_dofs(DP1) + num_free_dofs(DP2) + num_free_dofs(DP3)
  pdof, xdof = num_free_dofs(L2p), num_free_dofs(Y)
  meshtype = meshname[1:findlast("_", meshname).start-1]
  meshidx = meshname[findlast("_", meshname).start+1:end]
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = cbrt(max_cell_area)
  simpl && (h = cbrt(6max_cell_area / π))
  ncells = num_cells(model)

  function correct_displacement(Kh)
    diri_tags = ["diri_x", "diri_y", "diri_z", "diri_xy", "diri_xz", "diri_yz", "diri_xyz"]
    diri_masks = [(true, false, false), (false, true, false), (false, false, true),
      (true, true, false), (true, false, true), (false, true, true), (true, true, true)]
    reffe_U1 = ReferenceFE(lagrangian, VectorValue{3,Float64}, ords[2] + 1)
    H10 = FESpace(model, reffe_U1, conformity=:H1;
      dirichlet_tags=diri_tags, dirichlet_masks=diri_masks)
    H1g = TrialFESpace(H10, [Ũ, Ũ, Ũ, Ũ, Ũ, Ũ, Ū])
    au(UU, UV) = ∫(∇(UV) ⊙ ∇(UU))dΩ
    lu(UV) = ∫(∇(UV) ⊙ Kh)dΩ
    solve(AffineFEOperator(au, lu, H1g, H10))
  end

  function record_norms(λ, xh, dΩ)
    Uh, Kh1, Kh2, Kh3, Ph1, Ph2, Ph3, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2, Kh3), vectors_to_tensor(Ph1, Ph2, Ph3)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2, ∇ ⋅ Ph3)
    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ), compute_l2_norm(Kh, dΩ)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ), compute_l2_norm(ph, dΩ)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ)^2)

    write_line(filename,
      "$reffe_names,$meshtype,$meshidx,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λ,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    Uh1 = correct_displacement(Kh)
    Udof1 = length(Uh1.free_values)
    xdof1 = Udof1 + Kdof + Pdof + pdof
    Ul2norm1 = compute_l2_norm(Uh1, dΩ⁺)
    write_line(filename,
      "$(reffe_names)(corr),$meshtype,$meshidx,$ncells,$h,$Udof1,$Kdof,$Pdof,$pdof,$xdof1,",
      "$λi,$Ul2norm1,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = [
      "displacement" => Uh, "corrected displacement" => Uh1,
      "strain" => Kh, "stress" => Ph, "div stress" => ∇Ph,
      "pressure" => ph, "Jacobian (u = $λ)" => det(Kh + I3)
    ]
    writevtk(Ω, joinpath(DATA_DIR, "stretch3d_$(reffe_names)_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang(), iterations=50)
  solver = FESolver(nls)
  cache = nothing
  xh = zero(X)

  mo = maximum(ords)
  λs = sort!(unique([0:0.1/mo:0.3..., 0.3:0.05/mo:0.4..., 0.4:0.025/mo:0.5...]))
  println("\n[$reffe_names] $ncells cells optimisation begins...")
  for λi in λs
    λ[1] = λi
    lc[1] = λi < 0.4 ? 1 : 2
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)

    λi in (0.2, 0.5) && record_norms(λi, xh, dΩ⁺)
    println_now("λ = $λi optimisation finished\n..................................\n")
  end
end

pairs = ["P̄0d1d1P1", "P̄1d2d2P2"]
pair, meshidx = pairs[parse(Int, ARGS[1])], ARGS[2]
try
  run_stretch3d_option2(pair, "delaunay_$meshidx")
catch e
  println("[run_stretch3d_option2] pair $pair mesh $meshidx error: $e")
end
