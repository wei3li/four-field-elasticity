import Gridap: Measure, CellField, VectorValue, TensorValue, ⊙, ∑, ∫, ∇
import Gridap.FESpaces: SingleFieldFESpace, get_fe_dof_basis
import GridapDistributed: DistributedCellField, DistributedMeasure
import FillArrays: Fill
using Gridap.ReferenceFEs


const I2 = TensorValue([1. 0.; 0. 1.])
const I3 = TensorValue([1. 0. 0.; 0. 1. 0.; 0. 0. 1.])

const _L2d = VectorValue([1., 0.])
const _R2d = VectorValue([0., 1.])

const _L3d = VectorValue([1., 0., 0.])
const _M3d = VectorValue([0., 1., 0.])
const _R3d = VectorValue([0., 0., 1.])

function _PROJECT_DIR()
  curdir = pwd()
  stopidx = findfirst("HyperelasticityPaper", curdir).stop
  curdir[1:stopidx]
end

const DATA_DIR = joinpath(_PROJECT_DIR(), "notebooks/data")
const MESH_DIR = joinpath(_PROJECT_DIR(), "meshes")

vectors_to_tensor(v1, v2) = v1 ⊗ _L2d + v2 ⊗ _R2d
vectors_to_tensor(v1, v2, v3) = v1 ⊗ _L3d + v2 ⊗ _M3d + v3 ⊗ _R3d

scalars_to_vector(s1, s2) = s1 * _L2d + s2 * _R2d
scalars_to_vector(s1, s2, s3) = s1 * _L3d + s2 * _M3d + s3 * _R3d

symbol_to_space = Dict(
  "c̄" => o -> ReferenceFE(nedelec1, Float64, o),
  "c" => o -> ReferenceFE(nedelec2, Float64, o),
  "d̄" => o -> ReferenceFE(raviart_thomas, Float64, o),
  "d" => o -> ReferenceFE(bdm, Float64, o),
  "P̄" => (t, o) -> ReferenceFE(lagrangian, t, o, space=:P),
  "P" => (t, o) -> ReferenceFE(lagrangian, t, o, space=:P),
  "Q̄" => (t, o) -> ReferenceFE(lagrangian, t, o),
  "Q" => (t, o) -> ReferenceFE(lagrangian, t, o)
)

symbol_to_conformity = Dict(
  "c̄" => :HCurl,
  "c" => :HCurl,
  "d̄" => :HDiv,
  "d" => :HDiv,
  "P̄" => :L2,
  "P" => :H1,
  "Q̄" => :L2,
  "Q" => :H1,
)

function _construct_cellfield(u, dΩ::Measure)
  isa(u, CellField) && (return u)
  CellField(u, dΩ.quad.trian)
end

function _construct_cellfield(u, dΩ::DistributedMeasure)
  isa(u, DistributedCellField) && (return u)
  CellField(u, dΩ.trian)
end

function init_membrane_model1x1()
  ptr = [1, 5]
  data = [1, 2, 3, 4]
  cell_vertex_lids = Gridap.Arrays.Table(data, ptr)
  node_coords = Vector{Point{2,Float64}}(undef, 4)
  node_coords[1] = Point{2,Float64}(0, 0)
  node_coords[2] = Point{2,Float64}(48, 44)
  node_coords[3] = Point{2,Float64}(0, 44)
  node_coords[4] = Point{2,Float64}(48, 60)

  polytope = QUAD
  cell_types = collect(Fill(1, length(cell_vertex_lids)))
  cell_reffes = [ReferenceFE(polytope, lagrangian, Float64, 1)]
  grid = Gridap.Geometry.UnstructuredGrid(node_coords,
    cell_vertex_lids,
    cell_reffes,
    cell_types,
    Gridap.Geometry.NonOriented())
  model = Gridap.Geometry.UnstructuredDiscreteModel(grid)
  labels = get_face_labeling(model)
  labels.d_to_dface_to_entity[1] .= [1, 2, 1, 2]
  labels.d_to_dface_to_entity[2] .= [2, 2, 1, 2]
  labels.d_to_dface_to_entity[3] .= 3
  add_tag!(labels, "dirichlet", [1])
  add_tag!(labels, "neumann", [2])
  add_tag!(labels, "interior", [3])

  model
end

function compute_l2_norm(u, dΩ)
  uh = _construct_cellfield(u, dΩ)
  √(∑(∫(uh ⊙ uh)dΩ))
end

function compute_err_norm3d(typ, uhs, u, dΩ)
  uh = length(uhs) == 1 ? uhs[1] : vectors_to_tensor(uhs...)
  l2err = compute_l2_norm(uh - u, dΩ)
  if 'P' ∉ typ && 'Q' ∉ typ
    uh1, uh2, uh3 = uhs
    ∇eh = if 'c' in typ
      ∇eh1 = ∇ × (uh1 - col1 ∘ u)
      ∇eh2 = ∇ × (uh2 - col2 ∘ u)
      ∇eh3 = ∇ × (uh3 - col3 ∘ u)
      vectors_to_tensor(∇eh1, ∇eh2, ∇eh3)
    else
      scalars_to_vector(∇ ⋅ uh1, ∇ ⋅ uh2, ∇ ⋅ uh3) - ∇ ⋅ u
    end
    √(l2err^2 + compute_l2_norm(∇eh, dΩ)^2)
  elseif '̄' ∈ typ
    l2err
  else
    ∇eh = ∇(uh - u)
    √(l2err^2 + compute_l2_norm(∇eh, dΩ)^2)
  end
end

function println_now(msg)
  println(msg)
  flush(stdout)
end

function parse_range(s::String, T=Int32)
  (:)(parse.(T, split(s, ':'))...)
end

function create_file_with_header(filename, headers...)
  diridx = findlast('/', filename)
  if diridx !== nothing
    dirpath = filename[1:diridx-1]
    !ispath(dirpath) && mkpath(dirpath)
  end
  if !isfile(filename)
    open(filename, "w+") do io
      for h in headers
        write(io, string(h))
      end
      write(io, "\n")
    end
  end
end

function write_line(filename, args...)
  open(filename, "a") do io
    for arg in args
      write(io, string(arg))
    end
    write(io, "\n")
  end
end

function factorise_integer(x::Integer)
  factors = Integer[]
  i = 2
  while i * i <= x
    if (x % i == 0)
      while (x % i == 0)
        push!(factors, i)
        x = x ÷ i
      end
    end
    i += 1
  end
  x != 1 && push!(factors, x)
  factors
end
