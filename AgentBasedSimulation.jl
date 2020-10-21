### Main Simulator File

## Simulation parameters
# chemostat size (Inf=no chemostat)
chemo=1
# discretisation of histograms
dx=0.05
# time step of data collection
timeStep=2.
# transients
preStep=10.

# include the code
include("CHistograms.jl")
import .ConditionalHistograms
include("FirstDivisionAlgorithm.jl")
using .FirstDivisionAlgorithm

function main()

# initialise the model
model = cellModel()
# initialise state of cell population with one cell
# [1] size, [2] age, [3] growth rate, [4] birth size, [5] division size, [6] molecules,
state = CellType[]
push!(state,[1.,0.,1.,1.,model.getDivisionSize(1.),0.,0.])

# simulate some transients
println("Transients...")
FDAsim(model,state,preStep,chemostat=chemo)

# define histograms
hist_size     = ConditionalHistograms.Histo()
hist_mol      = ConditionalHistograms.Histo()
hist_mol_size = ConditionalHistograms.ConditionalHist(dx)

# run the simulation
println("Start simulation...")
global n=0
while true
  # simulation step
  FDAsim(model,state,timeStep,chemostat=chemo)

  # collect the data for analysis
  for cell in state
    size = [cell[1]]
    molNo = [cell[6]]
    # collect histograms
    hist_size.add(size)
    hist_mol.add(molNo)
    hist_mol_size.add(molNo,size)
  end

  # write stats to file periodically
  if n==10
    ConditionalHistograms.dumpHist(hist_size, "hist_size.dat")
    ConditionalHistograms.dumpHist(hist_mol,  "hist_protein.dat")
    ConditionalHistograms.dumpStats(hist_mol,  "stats_protein.dat")
    ConditionalHistograms.dumpStats(hist_mol_size, "hist_protein_size.dat")
    ConditionalHistograms.dumpStats(hist_mol_size, "stats_protein_size.dat")
    n=0
  end
  global n+=1

end # end while

end # end main

main()
