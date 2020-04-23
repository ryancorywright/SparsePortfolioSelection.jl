function cplex_raw_bigM(thePortfolio::SparsePortfolioData, ΔT_max::Float64=3600, gap::Float64=1e-4, debugLevel::Int=0)

  master_stock_problem_cplex=Model(solver=CplexSolver(CPX_PARAM_MIQCPSTRAT=1, CPX_PARAM_TILIM=ΔT_max, CPX_PARAM_EPGAP=gap, CPX_PARAM_THREADS=1))

  @variable(master_stock_problem_cplex, z[1:thePortfolio.n], Bin)
  @variable(master_stock_problem_cplex, x[1:thePortfolio.n]>=0)
  @constraint(master_stock_problem_cplex, sum(z[i] for i=1:thePortfolio.n)<=thePortfolio.k)
  @constraint(master_stock_problem_cplex, sparsity[i=1:thePortfolio.n], x[i]<=z[i])
  @constraint(master_stock_problem_cplex, minInvestment[i=1:thePortfolio.n], x[i]>=z[i]*thePortfolio.min_investment[i])
  @constraint(master_stock_problem_cplex, sumsToOne, sum(x[i] for i=1:thePortfolio.n)==1.0)
  @constraint(master_stock_problem_cplex, thePortfolio.A*x.>=thePortfolio.l)
  @constraint(master_stock_problem_cplex, thePortfolio.A*x.<=thePortfolio.u)
  @objective(master_stock_problem_cplex, Min, (0.5*sum((1.0/(thePortfolio.γ[i])*x[i]^2 for i=1:thePortfolio.n)) +0.5*x'*thePortfolio.X'*thePortfolio.X*x-thePortfolio.μ'*x+0.5*thePortfolio.Y'*thePortfolio.Y)')

  solve(master_stock_problem_cplex)

  optimalStocks=getvalue(x)

  ofv = getobjectivevalue(master_stock_problem_cplex)

  optimalIndicies=findall(a->a>0.5, getvalue(z))

  return CardinalityConstrainedPortfolio(optimalIndicies, optimalStocks[optimalIndicies], 0.0, zeros(thePortfolio.n), zeros(thePortfolio.m), zeros(thePortfolio.m), zeros(thePortfolio.m), getobjgap(master_stock_problem_cplex), getobjgap(master_stock_problem_cplex), 0, thePortfolio.μ'*optimalStocks, optimalStocks'*thePortfolio.X'*thePortfolio.X*optimalStocks, getsolvetime(master_stock_problem_cplex), getnodecount(master_stock_problem_cplex))

end

function cplex_raw_MISOCP(thePortfolio::SparsePortfolioData, ΔT_max::Float64=3600, gap::Float64=1e-4, debugLevel::Int=0)

  misocp_cplex=Model(solver=CplexSolver(CPX_PARAM_MIQCPSTRAT=1, CPX_PARAM_TILIM=ΔT_max, CPX_PARAM_EPGAP=gap, CPX_PARAM_THREADS=1))

  @variable(misocp_cplex, z[1:thePortfolio.n], Bin) # defines whether a stock is selected (big-M)
  @variable(misocp_cplex, x[1:thePortfolio.n]>=0.0)    # defines the quantity selected for each stock (short selling is permitted)
  @variable(misocp_cplex, θ[1:thePortfolio.n]>=0.0)
  @variable(misocp_cplex, τ>=0.0)

  @constraint(misocp_cplex, socpconstraint[i=1:thePortfolio.n], θ[i]+z[i]>=norm([2.0*x[i]; θ[i]-z[i]]))
  @constraint(misocp_cplex, sum(z[i] for i=1:thePortfolio.n)<=thePortfolio.k)
  @constraint(misocp_cplex, thePortfolio.A*x.>=thePortfolio.l)
  @constraint(misocp_cplex, thePortfolio.A*x.<=thePortfolio.u)
  @constraint(misocp_cplex, minInvestment[i=1:thePortfolio.n], x[i]>=z[i]*thePortfolio.min_investment[i])
  @constraint(misocp_cplex, sumsToOne, sum(x[i] for i=1:thePortfolio.n)==1.0)
  @constraint(misocp_cplex, τ+1.0>=norm([2.0*thePortfolio.X*x;τ-1.0]))
  @objective(misocp_cplex, Min, (0.5*sum(θ[i]/thePortfolio.γ[i] for i=1:thePortfolio.n)+0.5*τ-thePortfolio.μ'*x+0.5*thePortfolio.Y'*thePortfolio.Y)')

  solve(misocp_cplex)

  optimalStocks=getvalue(x)

  ofv = getobjectivevalue(misocp_cplex)
  optimalIndicies=findall(a->a>0.5, getvalue(z))
  numIters::Int128=0

  return CardinalityConstrainedPortfolio(optimalIndicies, optimalStocks[optimalIndicies], 0.0, zeros(thePortfolio.n), zeros(thePortfolio.m), zeros(thePortfolio.m),  zeros(thePortfolio.m), getobjgap(misocp_cplex), getobjgap(misocp_cplex), numIters, thePortfolio.μ'*optimalStocks, optimalStocks'*thePortfolio.X'*thePortfolio.X*optimalStocks, getsolvetime(misocp_cplex), getnodecount(misocp_cplex))

end

# This function is currently used for diagnostic purposes only (i.e. not part of the main code base)
function cplex_raw_MISOCP_subset(thePortfolio::SparsePortfolioData, z::Array{Float64})

  socp_mosek=Model(solver=MosekSolver(MSK_DPAR_INTPNT_QO_TOL_PFEAS=1e-6, MSK_DPAR_INTPNT_QO_TOL_DFEAS=1e-6, MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0))

  @variable(socp_mosek, x[1:thePortfolio.n]>=0.0)    # defines the quantity selected for each stock (short selling is permitted)
  @variable(socp_mosek, θ[1:thePortfolio.n]>=0.0)
  @variable(socp_mosek, τ>=0.0)

  @constraint(socp_mosek, socpconstraint[i=1:thePortfolio.n], θ[i]+z[i]>=norm([2.0*x[i]; θ[i]-z[i]]))
  @constraint(socp_mosek, sum(z[i] for i=1:thePortfolio.n)<=thePortfolio.k)
  @constraint(socp_mosek, thePortfolio.A*x.>=thePortfolio.l)
  @constraint(socp_mosek, thePortfolio.A*x.<=thePortfolio.u)
  @constraint(socp_mosek, minInvestment[i=1:thePortfolio.n], x[i]>=z[i]*thePortfolio.min_investment[i])
  @constraint(socp_mosek, sumsToOne, sum(x[i] for i=1:thePortfolio.n)==1.0)
  @constraint(socp_mosek, τ+1.0>=norm([2.0*thePortfolio.X*x;τ-1.0]))
  @objective(socp_mosek, Min, (0.5*sum(θ[i]/thePortfolio.γ[i] for i=1:thePortfolio.n)+0.5*τ-thePortfolio.μ'*x+0.5*thePortfolio.Y'*thePortfolio.Y)')

  solve(socp_mosek)

  ofv = getobjectivevalue(socp_mosek)
end



function cplex_raw_MISOCP_subset2(thePortfolio::SparsePortfolioData, z::Array{Float64})

  socp_mosek=Model(solver=MosekSolver(MSK_DPAR_INTPNT_QO_TOL_PFEAS=1e-6, MSK_DPAR_INTPNT_QO_TOL_DFEAS=1e-6, MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0))

  @variable(socp_mosek, x[1:thePortfolio.n]>=0.0)    # defines the quantity selected for each stock (short selling is permitted)
  @variable(socp_mosek, θ[1:thePortfolio.n]>=0.0)
  @variable(socp_mosek, τ>=0.0)

  @constraint(socp_mosek, socpconstraint[i=1:thePortfolio.n], θ[i]+z[i]>=norm([2.0*x[i]; θ[i]-z[i]]))
  @constraint(socp_mosek, sum(z[i] for i=1:thePortfolio.n)<=thePortfolio.k)
  @constraint(socp_mosek, thePortfolio.A*x.>=thePortfolio.l)
  @constraint(socp_mosek, thePortfolio.A*x.<=thePortfolio.u)
  @constraint(socp_mosek, minInvestment[i=1:thePortfolio.n], x[i]>=z[i]*thePortfolio.min_investment[i])
  @constraint(socp_mosek, sumsToOne, sum(x[i] for i=1:thePortfolio.n)==1.0)
  @constraint(socp_mosek, τ+1.0>=norm([2.0*thePortfolio.X*x;τ-1.0]))
  @objective(socp_mosek, Min, (0.5*sum(θ[i]/thePortfolio.γ[i] for i=1:thePortfolio.n)+0.5*τ-thePortfolio.μ'*x+0.5*thePortfolio.Y'*thePortfolio.Y)')

  theStatus=solve(socp_mosek)

  ofv = getobjectivevalue(socp_mosek)
  return theStatus, getvalue(x)
end



function cplex_MISOCP_relaxation(thePortfolio::SparsePortfolioData, ΔT_max::Float64=3600)

  misocp_cplex=Model(solver=CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_TILIM=ΔT_max))

  @variable(misocp_cplex, z[1:thePortfolio.n]>=0.0) # defines whether a stock is selected (big-M)
  @constraint(misocp_cplex, z.<=ones(thePortfolio.n))
  @variable(misocp_cplex, x[1:thePortfolio.n]>=0.0)    # defines the quantity selected for each stock (short selling is permitted)
  @variable(misocp_cplex, θ[1:thePortfolio.n]>=0.0)
  @variable(misocp_cplex, τ>=0.0)

  @constraint(misocp_cplex, socpconstraint[i=1:thePortfolio.n], θ[i]+z[i]>=norm([2.0*x[i]; θ[i]-z[i]]))
  @constraint(misocp_cplex, sum(z[i] for i=1:thePortfolio.n)<=thePortfolio.k)
  @constraint(misocp_cplex, thePortfolio.A*x.>=thePortfolio.l)
  @constraint(misocp_cplex, thePortfolio.A*x.<=thePortfolio.u)
  @constraint(misocp_cplex, minInvestment[i=1:thePortfolio.n], x[i]>=z[i]*thePortfolio.min_investment[i])
  @constraint(misocp_cplex, sumsToOne, sum(x[i] for i=1:thePortfolio.n)==1.0)
  @constraint(misocp_cplex, τ+1.0>=norm([2.0*thePortfolio.X*x;τ-1.0]))
  @objective(misocp_cplex, Min, (0.5*sum(θ[i]/thePortfolio.γ[i] for i=1:thePortfolio.n)+0.5*τ-thePortfolio.μ'*x+0.5*thePortfolio.Y'*thePortfolio.Y)')

  solve(misocp_cplex)

  optimalStocks=getvalue(x)

  ofv = getobjectivevalue(misocp_cplex)
  optimalIndicies=findall(a->a>0.5, getvalue(z))
  numIters::Int128=0

  return getvalue(z)

end
