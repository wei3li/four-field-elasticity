using Gridap, LineSearches
import GridapGmsh: GmshDiscreteModel
import Unicode: graphemes

include("utilities.jl")
include("MixedNedelec3d.jl")


λ, μ = [0.01], 1.0
T̄((x, y, z)) = x ≈ 48 ? VectorValue(0, λ[1], 0) : VectorValue(0, 0.0, 0)

XA = Point(48.0, 60.0, 5.0)

reffe_U = ReferenceFE(lagrangian, VectorValue{3,Float64}, 2)
reffe_K = NedelecMixedRefFE(Float64)
reffe_P = ReferenceFE(raviart_thomas, Float64, 0)
reffe_p = ReferenceFE(lagrangian, Float64, 0)

filename = joinpath(DATA_DIR, "membrane3d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "f,UxA,UyA,UzA,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")


function run_membrane3d_option0(meshname, α)
  model = GmshDiscreteModel(joinpath(MESH_DIR, "membrane3d/membrane3d_$(meshname).msh"))
  diri_tags = filter(s -> startswith(s, "diri"), model.face_labeling.tag_to_name)
  neum_tags = filter(s -> startswith(s, "neum"), model.face_labeling.tag_to_name)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri_xyz", filter(t -> contains(t, "_xyz"), diri_tags))
  add_tag_from_tags!(labels, "diri_z", filter(t -> contains(t, "_z"), diri_tags))

  H10 = TestFESpace(model, reffe_U, conformity=:H1,
    dirichlet_tags=["diri_xyz", "diri_z"], dirichlet_masks=[(true, true, true), (false, false, true)])
  Hc = FESpace(model, reffe_K)
  Hd = FESpace(model, reffe_P)
  L2 = FESpace(model, reffe_p, conformity=:L2)
  Y = MultiFieldFESpace([H10, Hc, Hc, Hc, Hd, Hd, Hd, L2])
  X = Y

  qdeg = 6
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)
  Γ = BoundaryTriangulation(Ω, tags=neum_tags)
  dΓ = Measure(Γ, qdeg)

  function a((U, K1, K2, K3, P1, P2, P3, p), (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))
    K, P = vectors_to_tensor(K1, K2, K3), vectors_to_tensor(P1, P2, P3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)

    F = K + I3
    J = det(F)
    pQ = p * inv(F)'
    δK = ∇(U) - K

    ∫(∇(V) ⊙ (P + α * δK) +
      κ ⊙ (μ * F + pQ - P - α * δK) +
      ψ ⊙ δK +
      q * (log ∘ (J))
    )dΩ
  end

  l((V, κ1, κ2, ψ1, ψ2, q)) = ∫(T̄ ⋅ V)dΓ
  r(x, y) = a(x, y) - l(y)

  function j(
    (U, K1, K2, K3, P1, P2, P3, p),
    (dU, dK1, dK2, dK3, dP1, dP2, dP3, dp),
    (V, κ1, κ2, κ3, ψ1, ψ2, ψ3, q))

    K = vectors_to_tensor(K1, K2, K3)
    dK, dP = vectors_to_tensor(dK1, dK2, dK3), vectors_to_tensor(dP1, dP2, dP3)
    κ, ψ = vectors_to_tensor(κ1, κ2, κ3), vectors_to_tensor(ψ1, ψ2, ψ3)

    F = K + I3
    iF = inv(F)
    iFdK = iF ⋅ dK
    δdK = ∇(dU) - dK

    ∫(∇(V) ⊙ (dP + α * δdK) +
      κ ⊙ (μ * dK - p * (iFdK ⋅ iF)' + dp * iF' - dP - α * δdK) +
      ψ ⊙ δdK +
      q * tr(iFdK)
    )dΩ
  end

  Udof, Kdof = num_free_dofs(H10), 3num_free_dofs(Hc)
  Pdof, pdof = 3num_free_dofs(Hd), num_free_dofs(L2)
  xdof = num_free_dofs(Y)
  meshtype = meshname[1:findlast("_", meshname).start-1]
  meshidx = meshname[findlast("_", meshname).start+1:end]
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = cbrt(6max_cell_area / π)

  function record_norms(λ, xh, dΩ)
    Uh, Kh1, Kh2, Kh3, Ph1, Ph2, Ph3, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2, Kh3), vectors_to_tensor(Ph1, Ph2, Ph3)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2, ∇ ⋅ Ph3)
    UxA, UyA, UzA = Uh(XA)
    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ), compute_l2_norm(Kh, dΩ)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ), compute_l2_norm(ph, dΩ)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ)^2)

    write_line(filename,
      "MixedNedelec,$meshtype,$meshidx,$ncell,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λ,$UxA,$UyA,$UzA,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = ["displacement" => Uh, "strain" => Kh, "stress" => Ph,
      "div stress" => ∇Ph, "pressure" => ph]
    writevtk(Ω, joinpath(DATA_DIR, "membrane3d_mixednedelec_results"), cellfields=fields)
  end

  nls = NLSolver(show_trace=true, method=:newton, linesearch=HagerZhang())
  solver = FESolver(nls)
  cache = nothing
  xh = zero(Y)

  ncell = num_cells(model)
  println("\n[MixedNedelec] $ncell cells optimisation begins...")
  for λi in [0:0.025:0.35..., 0.36:0.02:0.4...]
    λ[1] = λi
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)
    λi in (0.2, 0.4) && record_norms(λi, xh, dΩ⁺)
    println_now("λ = $λi optimisation finished\n..................................\n")
  end
end

meshidx = ARGS[1]
try
  run_membrane3d_option0("delaunay_$meshidx", 1e5)
catch e
  println("[run_membrane3d_option0] mesh $meshidx error: $e")
end
