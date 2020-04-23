
function portfolios_hillclimb(thePortfolio::SparsePortfolioData, indices_new::Array{Int64,1})

  L=1e1 #Note that this parameter should be tuned by the user if the warm-start isn't performing well
  λ=0.0
  iter=0
  indices = []
  α = zeros(thePortfolio.f)
  βl= zeros(size(thePortfolio.A, 1))
  βu= copy(βl)
  ρ=zeros(thePortfolio.n)

  SoRCoeff=1.7 # Improves performance of the method.

  while (iter < 20 || indices != indices_new && iter<=50)
    iter = iter+1
    indices = indices_new
    # Maximize over α for a given s
    dual_vars= inner_dual(thePortfolio, indices)
    ρ_full=zeros(thePortfolio.n)
    ρ_full[indices]=copy(dual_vars.ρ)
    # Averaging
    α = (α*(iter-SoRCoeff) + SoRCoeff*dual_vars.α)/iter
    λ = (λ*(iter-SoRCoeff) + SoRCoeff*dual_vars.λ)/iter
    βl= (βl*(iter-SoRCoeff) + SoRCoeff*dual_vars.βl)/iter
    βu= (βu*(iter-SoRCoeff) + SoRCoeff*dual_vars.βu)/iter
    ρ = (ρ*(iter-SoRCoeff) + SoRCoeff*ρ_full)/iter
    # Maximize over s for a given set of dual variables
    # Using a discrete first-order heuristic (see Bertsekas 1999, Bertsimas and King 2015, among others).

    p = -(thePortfolio.γ'/2.0).*(α'thePortfolio.X+(βl-βu)'*thePortfolio.A+λ*fill(1.0, thePortfolio.n)'+ρ'-thePortfolio.d').^2+(thePortfolio.min_investment.*ρ)' #technically speaking we should be using w, since w>=..., but this is a heuristic, so this doesn't affect correctness/performance
    w = (thePortfolio.γ[indices]).*(thePortfolio.X[:, indices]'*α+thePortfolio.A[:, indices]'*(βl-βu)+λ*fill(1.0, size(indices, 1))+ρ[indices]-thePortfolio.d[indices])
    x=fill(0.0,thePortfolio.n)
    x[indices]= w'
    indices_new = sort(sortperm((abs.(x'-(1/L)*p)[1:end]), rev=true)[1:thePortfolio.k]) # drops additional values past k here (e.g. if called by heuristiccallback and cplex wants k+1 nonzero)
    # Clean solution if infeasible due to min_investment constraints.
    while sum(thePortfolio.min_investment[indices_new])>1.0
      index=rand(indices_new)
      indices_new=filter(e->e≠index, indices_new)
    end
  end

  dual_vars= inner_dual(thePortfolio, indices_new)
  return indices_new, thePortfolio.γ'.*dual_vars.w
end
