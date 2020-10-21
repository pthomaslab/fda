using Distributions

# Extrande simulation algorithm of reactions in growing cells
function Extrande(model::cellModel,cell::Array{Float64}, tIn::Float64, tOut::Float64)
   time::Float64 = tIn
   while true

        # Compute the propensity bound (assuming propensities increase with cell size)
        propensityBound::Float64  = sum(model.propensities([cell[5]], cell[6:end]))
        # propose reaction time
        tau::Float64 = -log(rand())/propensityBound

        # exit in time
        if (time + tau) > tOut
           GrowSingleCell(cell, tOut - time)
           break
        end

        # grow cell until the next reaction occurs
        GrowSingleCell(cell, tau)
        # choose a reaction and execute it
        a = model.propensities(cell[1:5], cell[6:end])
        propensitySum::Float64 = sum(a)
        r::Float64 = rand()*propensityBound
        if propensitySum >= r
           reaction::Int = 1
           pSum::Float64 = a[1]
           while pSum < r
              reaction+=1
              pSum += a[reaction]
           end
           cell[6:end] += model.Stoichiometry[:,reaction]
        else
          # or do nothing (thinning step)
        end
        time += tau
    end

end

# apply Extrande to all cells
function GrowCells(model::cellModel, population::PopulationType, tIn, tOut)
   for cell in population
     Extrande(model, cell, tIn, tOut)
   end
end
