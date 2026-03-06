using Gridap, Gridap.ReferenceFEs, Gridap.Polynomials
using Gridap.Arrays
import Gridap.Polynomials: Monomial
import Gridap.Fields: FieldGradientArray

## branch: "https://github.com/gridap/Gridap.jl.git"#nedelec-new-kinds


struct NedelecMixedPrebasis <: AbstractVector{Monomial}
  order::Int

  function NedelecMixedPrebasis(order::Integer=1)
    new(Int(order))
  end
end

function Base.size(a::NedelecMixedPrebasis)
  15
end

Base.getindex(a::NedelecMixedPrebasis, i::Integer) = Monomial()
Base.IndexStyle(::Type{<:NedelecMixedPrebasis}) = IndexLinear()

Polynomials.return_type(::NedelecMixedPrebasis) = VectorValue{3,Float64}
Polynomials.num_terms(a::NedelecMixedPrebasis) = length(a)
Polynomials.get_order(a::NedelecMixedPrebasis) = a.order

function Arrays.return_cache(
  f::NedelecMixedPrebasis, x::AbstractVector{<:Point})
  np = length(x)
  ndofs = num_terms(f)
  V = eltype(x)
  a = zeros(V, (np, ndofs))
  P = MonomialBasis{3}(VectorValue{3,Float64}, f.order, (e, order) -> sum(e) <= order)
  cP = return_cache(P, x)
  CachedArray(a), cP, P
end

function Arrays.evaluate!(cache, f::NedelecMixedPrebasis, x::AbstractVector{<:Point})
  ca, cP, P = cache
  np = length(x)
  ndofs = num_terms(f)
  ndofsP = length(P)
  setsize!(ca, (np, ndofs))
  Px = evaluate!(cP, P, x)
  a = ca.array
  for (ip, p) in enumerate(x)
    for j in 1:ndofsP
      a[ip, j] = Px[ip, j]
    end
    x1, x2, x3 = p
    x1x2, x1x3, x2x3 = x1 * x2, x1 * x3, x2 * x3
    x1x2x3 = x1 * x2x3
    a[ip, ndofsP+1] = VectorValue(x2x3 - x2 * x2x3 - x2x3 * x3, x1x2x3, x1x2x3)
    a[ip, ndofsP+2] = VectorValue(x1x2x3, x1x3 - x1 * x1x3 - x1x3 * x3, x1x2x3)
    a[ip, ndofsP+3] = VectorValue(x1x2x3, x1x2x3, x1x2 - x1 * x1x2 - x1x2 * x2)
  end
  a
end

function Arrays.return_cache(
  g::FieldGradientArray{1,<:NedelecMixedPrebasis},
  x::AbstractVector{<:Point})
  f, D = g.fa, 3
  np = length(x)
  ndofs = num_terms(f)
  xi = testitem(x)
  V = eltype(x)
  G = Gridap.Fields.gradient_type(V, xi)
  a = zeros(G, (np, ndofs))
  mb = MonomialBasis{D}(VectorValue{D,Float64}, f.order, (e, order) -> sum(e) <= order)
  P = Broadcasting(∇)(mb)
  cP = return_cache(P, x)
  CachedArray(a), cP, P
end

function Arrays.evaluate!(
  cache,
  g::FieldGradientArray{1,<:NedelecMixedPrebasis},
  x::AbstractVector{<:Point})
  ca, cP, P = cache
  f = g.fa
  np = length(x)
  ndofs = num_terms(f)
  setsize!(ca, (np, ndofs))
  a = ca.array
  fill!(a, zero(eltype(a)))
  ndofsP = length(P)
  Px = evaluate!(cP, P, x)
  for (ip, p) in enumerate(x)
    for j in 1:ndofsP
      a[ip, j] = Px[ip, j]
    end
    x1, x2, x3 = p
    x1x2, x1x3, x2x3 = x1 * x2, x1 * x3, x2 * x3
    x12, x22, x32 = x1 * x1, x2 * x2, x3 * x3
    zv = zero(x1)
    a[ip, ndofsP+1] = TensorValue(
      zv, x3 - 2x2x3 - x32, x2 - x22 - 2x2x3,
      x2x3, x1x3, x1x2,
      x2x3, x1x3, x1x2
    )
    a[ip, ndofsP+2] = TensorValue(
      x2x3, x1x3, x1x2,
      x3 - 2x1x3 - x32, zv, x1 - x12 - 2x1x3,
      x2x3, x1x3, x1x2
    )
    a[ip, ndofsP+3] = TensorValue(
      x2x3, x1x3, x1x2,
      x2x3, x1x3, x1x2,
      x2 - 2x1x2 - x22, x1 - x12 - 2x1x2, zv
    )
  end
  a
end

function _NedelecMixed_nodes_and_moments(::Type{et}, p::Polytope) where {et}
  D = num_dims(p)
  ft, pt = VectorValue{D,et}, Point{D,et}

  nf_nodes = [zeros(pt, 0) for _ in 1:num_faces(p)]
  nf_moments = [zeros(ft, 0, 0) for _ in 1:num_faces(p)]

  ecips, emoments = ReferenceFEs._Nedelec_edge_values(p, et, 1)
  erange = get_dimrange(p, 1)
  nf_nodes[erange] = ecips
  nf_moments[erange] = emoments

  ccips, cmoments = ReferenceFEs._Nedelec_cell_values(p, et, 2)
  crange = get_dimrange(p, D)
  nf_nodes[crange] = ccips
  nf_moments[crange] = cmoments

  nf_nodes, nf_moments
end

function NedelecMixedRefFE(::Type{et}) where {et}
  p = TET
  prebasis = NedelecMixedPrebasis()

  nf_nodes, nf_moments = _NedelecMixed_nodes_and_moments(et, p)
  face_dofs = ReferenceFEs._face_own_dofs_from_moments(nf_moments)
  dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

  ndofs = num_dofs(dof_basis)
  metadata = nothing

  reffe = GenericRefFE{Nedelec}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    CurlConformity(),
    metadata,
    face_dofs)

  reffe
end
