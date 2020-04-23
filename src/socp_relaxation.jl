include("inner_problem.jl")


# qcqp version of relaxation: it's faster than the socp version
# please see the raw.jl file for the dual problem to this relaxation
function portfolios_socp(thePortfolio::SparsePortfolioData)

  n =thePortfolio.f
  m =thePortfolio.m
  f =thePortfolio.n 

  socpProblem = Model(solver=MosekSolver(MSK_DPAR_INTPNT_QO_TOL_PFEAS=1e-5, MSK_DPAR_INTPNT_QO_TOL_DFEAS=1e-5, MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0))

  @variable(socpProblem, α[1:n])
  @variable(socpProblem, λ)
  @variable(socpProblem, t>=0)
  @variable(socpProblem, v[1:f]>=0.0)
  @variable(socpProblem, w[1:f])
  @variable(socpProblem, βl[1:m]>=0.0)
  @variable(socpProblem, βu[1:m]>=0.0)
  @variable(socpProblem, ρ[1:n]>=0.0)

  @constraint(socpProblem, w.>=thePortfolio.X'*α+thePortfolio.A'*(βl-βu)+λ*ones(f)+ρ-thePortfolio.d)
  @constraint(socpProblem, knorm[i=1:f], v[i]+t+ρ[i]*thePortfolio.min_investment[i]>=(thePortfolio.γ[i]/2.0)*w[i]^2)

  @objective(socpProblem, Max, (-0.5*α'*α+thePortfolio.Y'*α+λ+βl'*thePortfolio.l-βu'*thePortfolio.u-sum(v)-thePortfolio.k*t)') #+π'*thePortfolio.min_investment

  status=solve(socpProblem)

  return Dual(getvalue(α), getvalue(λ), getvalue(βl), getvalue(βu), getvalue(ρ), getvalue(w), getobjectivevalue(socpProblem), status), getobjectivevalue(socpProblem)
end

# SOCP version of relaxation. Not in use, as this version is currently slower.
function portfolios_socp2(thePortfolio::SparsePortfolioData)

  n = size(thePortfolio.X, 1)
  m = size(thePortfolio.A, 1)
  f = size(thePortfolio.X, 2)

  socpProblem = Model(solver=MosekSolver( MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-8, MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-8))

  @variable(socpProblem, α[1:n])
  @variable(socpProblem, λ)
  @variable(socpProblem, t>=0)
  @variable(socpProblem, v[1:f]>=0.0) 
  @variable(socpProblem, w[1:f])
  @variable(socpProblem, βl[1:m]>=0.0)
  @variable(socpProblem, βu[1:m]>=0.0)
  @variable(socpProblem, ρ[1:n]>=0.0)
  @variable(socpProblem, τ>=0.0)

  @constraint(socpProblem, w.>=thePortfolio.X'*α+thePortfolio.A'*(βl-βu)+λ*ones(f)+ρ-thePortfolio.d)
  @constraint(socpProblem, knorm[i=1:f], v[i]+t+ρ[i]*thePortfolio.min_investment[i]+1.0>=norm([sqrt(2.0*thePortfolio.γ[i])*w[i];v[i]+t+ρ[i]*thePortfolio.min_investment[i]-1.0]))
  @constraint(socpProblem, τ+1.0>=norm([2.0*α;τ-1.0]))

  @objective(socpProblem, Max, (-0.5*τ+thePortfolio.Y'*α+λ+βl'*thePortfolio.l-βu'*thePortfolio.u-sum(v)-thePortfolio.k*t)') #+π'*thePortfolio.min_investment

  status=solve(socpProblem)

  return Dual(getvalue(α), getvalue(λ), getvalue(βl), getvalue(βu), getvalue(ρ), getvalue(w), status, getobjectivevalue(socpProblem)), getobjectivevalue(socpProblem)
end
