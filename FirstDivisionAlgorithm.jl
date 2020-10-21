# the cell state consists of
# [1] size, [2] age, [3] growth rate, [4] birth-volume, [5] division size, [6:end] molecules, .

module FirstDivisionAlgorithm

export FDAsim, cellModel, CellType, PopulationType

# some type abbreviations
const CellType = Array{Float64}
const PopulationType = Array{CellType}

include("ABM_model.jl") # model definitions are contained in here
include("ABM_extrande.jl")
include("ABM_fda.jl")

end
