# some helper functions to collect conditional (and unconditional) histograms andstatistics
module ConditionalHistograms

export Histo, ConditionalHist, dumpHist, dumpStats

import Base.isless
methods(isless)

# vector comparison
function isless(x::Vector{Int64},y::Vector{Int64})
   for (i,val) in enumerate(y)
     if x[i]<val
       return true
     end
   end
   return false
end

# sorts a dictionary
function sorting(dict::AbstractDict)
  return sort(collect(zip(keys(dict),values(dict))))
end

# binning functions
function trim(x::Float64,dx)
   return Int(floor(x/dx))
end
function trim(x::Array{Float64},dx)
   return map((y)->trim(y,dx),x)
end

abstract type AbstractStatsContainer end

# An easy histogram container which also computes mean and covariance
mutable struct Histo <: AbstractStatsContainer
  scale::Float64
  count::BigInt
  data::Dict{Vector{Int64},Int64}
  add::Function
  mean::Vector{Float64}
  #cov::Matrix{Float64}
  cov::Array{Float64}
  getMean::Function
  getCovariance::Function
  getHistogram::Function

  function Histo(scale::Float64=1.)

      this = new()
      this.count = 0
      this.data  = Dict{Array{Int64},Int64}()
      this.scale = scale

      # add a value to the histogram
      this.add  = function(values::Array{Float64})

            # discretise values
            val = trim(values,this.scale)
            # update histogram
            if get(this.data,val,0)==0
              push!(this.data,val => 1)
            else
             this.data[val]+=1
            end

            # compute mean and mean square
            if this.count > 0
              for i in 1:length(values)
                this.mean[i]  = (this.mean[i]*this.count+values[i])/(this.count+1)
                for j in 1:i
                   this.cov[i,j] = (this.cov[i,j]*this.count+values[i]*values[j])/(this.count+1)
                end
              end
            else
              this.mean=values
              this.cov=convert(Array{Float64}, [length(values),length(values)])
              for i in 1:length(values)
                for j in 1:i
                   this.cov[i,j] = (values[i]*values[j])
                end
              end
            end
            # increase count
            this.count+=1

      end # end add

    # mean
    this.getMean = function()
       return this.mean
    end

    # covariance matrix
    this.getCovariance = function()
      co = Vector{Float64}()
      for i in 1:length(this.mean)
        for j in 1:i
          v = (this.cov[i,j] - this.mean[i]*this.mean[j])
          push!(co, v)
        end
      end
      return co
    end

    this.getHistogram = function()
       # implement normalization
    end


    return this
  end #function Histo()

end #end histo

# Container which only computes mean and covariance
mutable struct StatsContainer  <: AbstractStatsContainer
  count::BigInt
  add::Function
  mean::Vector{Float64}
  cov::Matrix{Float64}
  getMean::Function
  getCovariance::Function

  function StatsContainer()

      this = new()
      this.count = 0

      this.add  = function(values::Vector{Float64})

            # compute mean and mean square
            if this.count > 0
              for i in 1:length(values)
                this.mean[i]  = (this.mean[i]*this.count+values[i])/(this.count+1)
                for j in 1:i
                   this.cov[i,j] = (this.cov[i,j]*this.count+values[i]*values[j])/(this.count+1)
                end
              end
            else
              this.mean=val
              this.cov=Array{Float64}(undef,length(values),length(values))
              for i in 1:length(values)
                for j in 1:i
                   this.cov[i,j] = (values[i]*values[j])
                end
              end
            end

            # increase count
            this.count+=1

      end # end add

    #
    this.getMean = function()
       return this.mean
    end

    this.getCovariance = function()
      co = Vector{Float64}()
      for i in 1:length(this.mean)
        for j in 1:i
          v = (this.cov[i,j] - this.mean[i]*this.mean[j])
          push!(co, v)
        end
      end
      return co
    end

    return this
    end

end #end StatsContainer

abstract type AbstractConditionalStatsContainer end

# container that records conditional histograms p(X|Y)
mutable struct ConditionalHist <: AbstractConditionalStatsContainer
   scaleX::Float64
   scaleY::Float64
   data::Dict{Vector{Int64},Histo}
   add::Function

   function ConditionalHist(dy::Float64=1.,dx::Float64=1.)
     this = new()
     this.scaleY = dy
     this.scaleX = dx
     this.data  = Dict{Vector{Int64},Histo}()

     this.add = function(values::Array{Float64},keys)
       ukey = trim(keys,this.scaleY)
       if get(this.data,ukey,0)==0
         push!(this.data,ukey => Histo(this.scaleX))
       end
       this.data[ukey].add(values)
     end

     return this
   end

end # end ConditionalHist

# container that records conditional statistics (mean & covariance matrix)
mutable struct ConditionalStats <: AbstractConditionalStatsContainer
   data::Dict{Vector{Int64},StatsContainer}
   add::Function
   scaleY::Float64

   function ConditionalStats()
     this = new()
     this.data = Dict{Vector{Int64},StatsContainer}()

     this.add = function(values,keys)
       ukey = trim(keys,this.scaleY)
       if get(this.data,ukey,0)==0
         push!(this.data,ukey => StatsContainer())
       end
       this.data[ukey].add(values)
     end

     return this
   end
end # end ConditionalStats

function dumpStats(h::AbstractStatsContainer, fname::AbstractString)
    f=open(fname,"w")
    for m in h.getMean()
        write(f,"\t $m")
    end
    cov = h.getCovariance()
    for v in cov
        write(f,"\t $v")
    end
    write(f,"\n")
    close(f)
end

function dumpStats(h::AbstractConditionalStatsContainer, fname::AbstractString, shiftY::Float64=0.)
    f=open(fname,"w")
    for (key,val) in sorting(h.data)
      for k in key
        keyVal=h.scaleY*(k+shiftY)
        write(f,"$keyVal \t")
      end
      for m in val.getMean()
        write(f,"\t $m")
      end
      cov = val.getCovariance()
      for v in cov
        write(f,"\t $v")
      end
      write(f,"\n")
    end
    close(f)
end

function dumpHist(h::Histo, fname::AbstractString, shiftX::Float64=0.)
  f=open(fname,"w")
  for (state,prob) in sorting(h.data)
    for s in state
      keyVal=h.scale*(s+shiftX)
      write(f,"\t $keyVal")
    end
    p=prob/h.count/(h.scale^length(state))
    write(f,"\t $p \n")
  end
  close(f)
end

function dumpHist(h::ConditionalHist, fname::AbstractString, shiftX::Float64=0., shiftY::Float64=0.)

    f=open(fname,"w")
    for (key,val) in sorting(h.data)
      for k in key
        keyVal=h.scaleY*(k+shiftY)
        for (state,prob) in sorting(val.data)
          write(f,"$keyVal")
          for s in state
            st=h.scaleX*(s+shiftX)
            write(f,"\t $st")
          end
          norm=val.count*(h.scaleX^length(state))
          p=prob/norm
          write(f,"\t $p \n")
        end
      end
    end
    close(f)

end

end
