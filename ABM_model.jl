### Model definitions
using Distributions
mutable struct cellModel
    slope::Float64 # slope of linear growth model
    cv2p::Float64 # CV^2 of division errors
    cv2size::Float64 # CV^2 of cell size control errors
    noiseDist::ContinuousUnivariateDistribution # distribution
    pDist:: ContinuousUnivariateDistribution

    SpeciesNames::Array{String} # names of all species
    Stoichiometry::Any # stoichiometric matrix
    propensities::Function
    getDivisionSize::Function
    getTheta::Function
    getGrowthRate::Function

    function cellModel()
        this = new()

        # set slope of linear growth model
        this.slope = 1.0 #(adder=1, sizer=0.)
        # set the partition distribution
        this.cv2p =(0.01)^2
        this.pDist = Beta((1-this.cv2p)/(2*this.cv2p),(1-this.cv2p)/(2*this.cv2p))
        # sets the added size distribution
        this.cv2size = 0.075
        this.noiseDist = Gamma((2-this.slope)/this.cv2size,this.cv2size)

        # set Stoichiometry and Propensities
        this.SpeciesNames  = ["mRNA","Protein"]
        this.Stoichiometry = [ [1. , 0.] [-1. , 0.] [0. , 1.] ]
        this.propensities = function(extState::Array{Float64}, molState::Array{Float64})
           # rate constants
           k0::Float64 = 10.; kdm::Float64 = 9.; ks::Float64 = 100.;
           return [ k0*extState[1] ; kdm*molState[1] ; ks*molState[1] ]
        end

        # cell size control model
        this.getDivisionSize = function(birthSize::Float64)
            return this.cv2size>0. ? this.slope*birthSize + rand(this.noiseDist) : this.slope*birthSize + (2-this.slope)
        end
        # division errors
        this.getTheta = function()
            return this.cv2p>0. ? rand(this.pDist) : 0.5
        end
        # growth rate
        this.getGrowthRate = function(x)
            return x
        end

        return this
    end


end
