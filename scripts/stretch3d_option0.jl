using Gridap
import GridapGmsh: GmshDiscreteModel
import LineSearches: HagerZhang, MoreThuente
import Unicode: graphemes

include("utilities.jl")
include("MixedNedelec3d.jl")


λ, μ = [0.1], 1.0
Ū((x, y, z)) = VectorValue(0.0, 0.0, z > 0.498 ? λ[1] : 0.0)
Ũ(_) = VectorValue(0.0, 0.0, 0.0)

reffe_U = ReferenceFE(lagrangian, VectorValue{3,Float64}, 2)
reffe_K = NedelecMixedRefFE(Float64)
reffe_P = ReferenceFE(raviart_thomas, Float64, 0)
reffe_p = ReferenceFE(lagrangian, Float64, 0)

filename = joinpath(DATA_DIR, "stretch3d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "u,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")


function run_stretch3d_option0(meshname, α)
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

  diri_tags = ["diri_x", "diri_y", "diri_z", "diri_xy", "diri_xz", "diri_yz", "diri_xyz"]
  diri_masks = [(true, false, false), (false, true, false), (false, false, true),
    (true, true, false), (true, false, true), (false, true, true), (true, true, true)]

  H10 = TestFESpace(model, reffe_U, conformity=:H1; dirichlet_tags=diri_tags, dirichlet_masks=diri_masks)
  Hc = FESpace(model, reffe_K)
  Hd = FESpace(model, reffe_P)
  L2 = FESpace(model, reffe_p, conformity=:L2)
  Y = MultiFieldFESpace([H10, Hc, Hc, Hc, Hd, Hd, Hd, L2])

  qdeg = 6
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)

  function a((U, K1, K2, K3, P1, P2, P3, p), (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    K, P = vectors_to_tensor(K1, K2, K3), vectors_to_tensor(P1, P2, P3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)

    F = K + I3
    J = det(F)
    pQ = p * J * inv(F)'
    δK = ∇(U) - K

    ∫(∇(V) ⊙ (P + α * δK) +
      κ ⊙ (μ * F + pQ - P - α * δK) +
      ψ ⊙ δK +
      q * (J - 1)
    )dΩ
  end

  r(x, y) = a(x, y)

  function j(
    (U, K1, K2, K3, P1, P2, P3, p),
    (dU, dK1, dK2, dK3, dP1, dP2, dP3, dp),
    (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))

    K = vectors_to_tensor(K1, K2, K3)
    dK, dP = vectors_to_tensor(dK1, dK2, dK3), vectors_to_tensor(dP1, dP2, dP3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)

    F = K + I3
    J, iF = det(F), inv(F)
    iFdK = iF ⋅ dK
    JiFdK = J * tr(iFdK)
    δdK = ∇(dU) - dK

    ∫(∇(V) ⊙ (dP + α * δdK) +
      κ ⊙ (μ * dK + p * (JiFdK * iF' - J * (iFdK ⋅ iF)') + dp * J * iF' - dP - α * δdK) +
      ψ ⊙ δdK +
      q * JiFdK
    )dΩ
  end

  Udof, Kdof = num_free_dofs(H10), 3num_free_dofs(Hc)
  Pdof, pdof = 3num_free_dofs(Hd), num_free_dofs(L2)
  xdof = num_free_dofs(Y)
  meshtype = meshname[1:findlast("_", meshname).start-1]
  meshidx = meshname[findlast("_", meshname).start+1:end]
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = cbrt(6max_cell_area / π)
  ncells = num_cells(model)

  function record_norms(λ, xh, dΩ)
    Uh, Kh1, Kh2, Kh3, Ph1, Ph2, Ph3, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2, Kh3), vectors_to_tensor(Ph1, Ph2, Ph3)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2, ∇ ⋅ Ph3)

    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ), compute_l2_norm(Kh, dΩ)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ), compute_l2_norm(ph, dΩ)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ)^2)

    write_line(filename,
      "MixedNedelec,$meshtype,$meshidx,$ncells,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λ,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = ["displacement" => Uh, "strain" => Kh, "stress" => Ph,
      "div stress" => ∇Ph, "pressure" => ph, "Jacobian (u = $λ)" => det(Kh + I3)]
    writevtk(Ω, joinpath(DATA_DIR, "stretch3d_mixednedelec_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang(), ftol=1e-10)
  solver = FESolver(nls)
  cache = nothing
  xh = zero(Y)

  λs = sort!(unique([0:0.1:0.3..., 0.3:0.05:0.4..., 0.4:0.025:0.5...]))
  println("\n[MixedNedelec] $ncells cells optimisation begins...")
  for λi in λs
    λ[1] = λi
    X = MultiFieldFESpace([TrialFESpace(H10, [Ũ, Ũ, Ũ, Ũ, Ũ, Ũ, Ū]), Hc, Hc, Hc, Hd, Hd, Hd, L2])
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)
    λi in (0.2, 0.5) && record_norms(λi, xh, dΩ⁺)
    println_now("λ = $λi optimisation finished\n..................................\n")
  end
end

meshidx = ARGS[1]
try
  run_stretch3d_option0("delaunay_$meshidx", 1e5)
catch e
  println("[run_stretch3d_option0] mesh $meshidx error: $e")
end
