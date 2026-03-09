using Gridap, LineSearches
import GridapGmsh: GmshDiscreteModel
import Unicode: graphemes

include("utilities.jl")


λ, μ = [1.1], 1.0

Ū((x, y)) = VectorValue(x < 0.5 ? 0.0 : λ[1], 0.0)

filename = joinpath(DATA_DIR, "stretch2d_record.csv")
create_file_with_header(filename,
  "comb,meshtype,meshidx,ncell,h,Udof,Kdof,Pdof,pdof,xdof,",
  "u,Ul2norm,Kl2norm,Pl2norm,Phdivnorm,pl2norm")


function run_stretch2d_option0(reffe_names, meshname)
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

  model = GmshDiscreteModel(joinpath(MESH_DIR, "stretch2d/stretch2d_$(meshname).msh"))

  tag_to_name = get_face_labeling(model).tag_to_name
  diri_tags = filter(s -> startswith(s, "diri"), tag_to_name)
  confs = [symbol_to_conformity[typs[i]] for i in 1:4]
  H1U = FESpace(model, reffe_U, conformity=confs[1]; dirichlet_tags=diri_tags)
  CuK = FESpace(model, reffe_K, conformity=confs[2])
  DvP = TestFESpace(model, reffe_P, conformity=confs[3])
  L2p = FESpace(model, reffe_p, conformity=confs[4])
  Y = MultiFieldFESpace([H1U, CuK, CuK, DvP, DvP, L2p])

  qdeg = 2max(maximum(ords), 1) + 1
  Ω = Triangulation(model)
  dΩ, dΩ⁺ = Measure(Ω, qdeg), Measure(Ω, 2qdeg)

  function a((U, K1, K2, P1, P2, p), (V, κ1, κ2, ψ1, ψ2, q))
    K, P = vectors_to_tensor(K1, K2), vectors_to_tensor(P1, P2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)

    F = K + I2
    J = det(F)
    pQ = p * J * inv(F)'

    ∫(∇(V) ⊙ P +
      κ ⊙ (μ * F + pQ - P) +
      ψ ⊙ (∇(U) - K) +
      q * (J - 1))dΩ
  end

  function j((U, K1, K2, P1, P2, p), (dU, dK1, dK2, dP1, dP2, dp), (V, κ1, κ2, ψ1, ψ2, q))
    K = vectors_to_tensor(K1, K2)
    dK, dP = vectors_to_tensor(dK1, dK2), vectors_to_tensor(dP1, dP2)
    κ, ψ = vectors_to_tensor(κ1, κ2), vectors_to_tensor(ψ1, ψ2)

    F = K + I2
    J, iF = det(F), inv(F)
    iFdK = iF ⋅ dK
    JiFdK = J * tr(iFdK)

    ∫(∇(V) ⊙ dP +
      κ ⊙ (μ * dK + p * (JiFdK * iF' - J * (iFdK ⋅ iF)') + dp * J * iF' - dP) +
      ψ ⊙ (∇(dU) - dK) +
      q * JiFdK)dΩ
  end

  r(x, y) = a(x, y)

  Udof, Kdof = num_free_dofs(H1U), 2num_free_dofs(CuK)
  Pdof, pdof = 2num_free_dofs(DvP), num_free_dofs(L2p)
  xdof = num_free_dofs(Y)
  meshtyp = meshname[1:findlast("_", meshname).start-1]
  meshidx = meshname[findlast("_", meshname).start+1:end]
  max_cell_area = maximum((∫(1)Measure(Ω, 2))[Ω])
  h = sqrt(max_cell_area)
  simpl && (h = sqrt(4 * max_cell_area / sqrt(3)))

  function record_norms(λ, xh, dΩ)
    Uh, Kh1, Kh2, Ph1, Ph2, ph = xh
    Kh, Ph = vectors_to_tensor(Kh1, Kh2), vectors_to_tensor(Ph1, Ph2)
    ∇Ph = scalars_to_vector(∇ ⋅ Ph1, ∇ ⋅ Ph2)
    Ul2norm, Kl2norm = compute_l2_norm(Uh, dΩ), compute_l2_norm(Kh, dΩ)
    Pl2norm, pl2norm = compute_l2_norm(Ph, dΩ), compute_l2_norm(ph, dΩ)
    Phdivnorm = √(Pl2norm^2 + compute_l2_norm(∇Ph, dΩ)^2)

    write_line(filename,
      "$reffe_names,$meshtyp,$meshidx,$ncell,$h,$Udof,$Kdof,$Pdof,$pdof,$xdof,",
      "$λ,$Ul2norm,$Kl2norm,$Pl2norm,$Phdivnorm,$pl2norm")

    fields = [
      "displacement (u = $λ)" => Uh, "strain" => Kh, "stress" => Ph,
      "div stress" => ∇Ph, "pressure" => ph, "Jacobian (u = $λ)" => det(Kh + I2)
    ]
    writevtk(Ω, joinpath(DATA_DIR, "stretch2d_$(reffe_names)_results"), cellfields=fields)
  end

  ftol = 1e-9
  nls = NLSolver(show_trace=true, method=:newton, ftol=ftol, linesearch=HagerZhang(), iterations=60)
  solver = FESolver(nls)
  cache = nothing
  xh = zero(Y)

  ncell = num_cells(model)
  println("\n[$reffe_names] $ncell cells optimisation begins...")
  for λi in [0:0.1:0.5..., 0.55:0.05:1.0..., 1.02:0.02:1.5...]
    λ[1] = λi
    X = MultiFieldFESpace([TrialFESpace(H1U, Ū), CuK, CuK, DvP, DvP, L2p])
    xh = FEFunction(X, get_free_dof_values(xh))
    op = FEOperator(r, j, X, Y)
    println("λ = $λi optimisation starting...\n")
    xh, cache = solve!(xh, solver, op, cache)

    rmax = norm(assemble_vector(y -> r(xh, y), Y), Inf)
    if rmax > ftol
      println("λ = $λi optimisation did not converge, inf norm of residual: $rmax.")
      println_now("..................................\n")
      continue
    end

    λi in 0.5:0.5:2.0 && record_norms(λi, xh, dΩ⁺)
    println_now("λ = $λi optimisation finished\n..................................\n")
  end
end

mtype = "delaunay"
pairs = ["P1c1d̄0P̄0", "P2c2d̄1P̄1"]
pair_to_mindices = Dict(
  pairs[1] => [80, 66, 53, 43, 35, 29],
  pairs[2] => [78, 60, 48, 39, 31])

pair = pairs[parse(Int, ARGS[1])]
for midx in pair_to_mindices[pair]
  try
    run_stretch2d_option0(pair, "$(mtype)_$midx")
  catch e
    println_now("[$pair, $mtype size $midx] failed: $e")
  end
end
