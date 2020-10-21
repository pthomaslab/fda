
# advances size/age of cells
function GrowSingleCell(cell::CellType, tau::Float64)
    cell[1] = cell[1] * exp(cell[3]*tau)
    cell[2] += tau
end

# computes the time interval from now to division
function getDivisionTime(cell::CellType)
  return log( cell[5]/cell[1] )/cell[3]
end

# computes the next dividing cell and its division time
function getNextDivision(inState::PopulationType)
    nextCell = -1; nextTau = 0.; idx = 1;
    for cell in inState
        tau = getDivisionTime(cell)
        if ( tau < nextTau || nextCell == -1 )
           nextTau = tau
           nextCell = idx
        end
        idx+=1
    end
    return [ nextTau, nextCell ]
end

# Binomial partitioning of molecules
function PartitionMolecules(molecules::Array{Float64}, p)
  return [rand(Binomial(Int(m),p)) for m in molecules]
end

# Partition a cell state into mother and daughter cells
function partitionCells(model::cellModel,cell::CellType)

    # partition cell size
    pDiv::Float64 = model.getTheta()
    daughterSize=cell[1]*pDiv;
    motherSize=cell[1]*(1 .-pDiv);
    # partition molecules
    mols = cell[6:end]
    daughterMols = PartitionMolecules(mols,pDiv)
    motherMols = mols - daughterMols;
    # assign new state vectors
    daughterState = [daughterSize; 0.; model.getGrowthRate(cell[3]); daughterSize; model.getDivisionSize(daughterSize); daughterMols]
    motherState   = [motherSize;   0.; model.getGrowthRate(cell[3]); motherSize;   model.getDivisionSize(motherSize); motherMols]

    return motherState, daughterState

end

# runs the First-Division Algorithm on a cell population
function FDAsim(model::cellModel, inState::PopulationType, tOut::Float64; chemostat=Inf)

time::Float64=0.
while time<tOut

#	Get division time and next cell dividing
	tau::Float64, cellIdx::Int = getNextDivision(inState)
#	Grow all cells and exit
	if (time+tau > tOut)
	   GrowCells(model, inState, time, tOut)
       break
    end
#  Grow all cells
   GrowCells(model, inState, time, time+tau)
#  Partition molecules
   motherState, daughterState = partitionCells(model, inState[cellIdx])
#  Give birth to daughter
   push!(inState, daughterState)
#  Divide mother
   inState[cellIdx] .= motherState
#  Update time
   time+=tau
#  Apply a chemostat
   if length(inState)>chemostat
	  splice!(inState, rand(1:length(inState))) # removes a random cell
   end

end # end while

end # end function FDAsim()
