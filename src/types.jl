struct CutIterData #used for plotting performance of method
    time::Float64
    Cut::Int128
    obj::Float64
    bestbound::Float64
    boundgap::Float64
    status::Symbol
end

mutable struct Cut
    p::Float64
    ∇s::Array{Float64}
    status::Symbol
  end

  struct Dual
    α::Array{Float64}
    λ::Float64
    βl::Array{Float64}
    βu::Array{Float64}
    ρ::Array{Float64}
    w::Array{Float64}
    ofv::Float64
    status::Symbol
  end

  struct SparsePortfolioData
    μ::Array{Float64}
    Y::Array{Float64}
    d::Array{Float64}
    X::Array{Float64, 2}
    A::Array{Float64, 2}
    l::Array{Float64}
    u::Array{Float64}
    k::Int64
    n::Int64 # No stocks
    m::Int64 # No constraints
    f::Int64 # Rank of cov matrix
    min_investment::Array{Float64, 1}
    γ::Array{Float64}
  end

struct CardinalityConstrainedPortfolio
  indices::Array{Int}
  w::Array{Float64}
  λ::Float64
  α::Array{Float64}
  βl::Array{Float64}
  βu::Array{Float64}
  ρ::Array{Float64}
  boundGap::Float64
  boundGapSOCP::Float64
  numIters::Int128
  expectedReturn::Float64
  portfolioVariance::Float64
  solveTime::Float64
  nodeCount::Float64
end

struct PortfolioData #Old struct, not currently in use.
  μ::Array{Float64}
  X::Array{Float64, 2}
  normalization::Array{Float64, 2}
  X_normalized::Array{Float64, 2}
  min_investment::Array{Float64, 1}
end
