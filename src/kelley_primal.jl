function getKelleyPrimalCuts(thePortfolio::SparsePortfolioData, useCopyOfVariables::Bool, minReturnConstraint::Bool, theStabilizationPoint::Array{Float64}, kelleyPrimalEpochs::Int64)

    rootModel = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
    @variable(rootModel, sRoot[1:thePortfolio.n]>=0.0)
    @constraint(rootModel, sRoot.<=ones(thePortfolio.n))
    @variable(rootModel, tRoot>=-1e12)

    # Objective
    @objective(rootModel, Min, tRoot)

    # Constraints
    @constraint(rootModel, sum(sRoot)<=thePortfolio.k)
    @constraint(rootModel, sum(sRoot)>=1)

    # Max investment constraint
    if useCopyOfVariables
      @constraint(rootModel, sum(sRoot[i]*(thePortfolio.A[1,i].>thePortfolio.l[1]) for i=1:thePortfolio.n)>=1.0) #ensures feasibility of problem
      @variable(rootModel,   xRoot[1:thePortfolio.n]>=0.0)
      @constraint(rootModel, thePortfolio.A*xRoot.>=thePortfolio.l)
      @constraint(rootModel, thePortfolio.A*xRoot.<=thePortfolio.u)
      @constraint(rootModel, minInvestment[i=1:thePortfolio.n], xRoot[i]>=sRoot[i]*thePortfolio.min_investment[i])
      @constraint(rootModel, bigMMaster[i=1:thePortfolio.n], xRoot[i]<=thePortfolio.u[i+1].*sRoot[i])
      @constraint(rootModel, sumsToOne, sum(xRoot[i] for i=1:thePortfolio.n)==1.0)
    end

    if minReturnConstraint #This line assumes that minimum return constraint comes first whenever it is included
      @constraint(rootModel, sum(sRoot[i]*(thePortfolio.A[1,i].>thePortfolio.l[1]) for i=1:thePortfolio.n)>=1.0) #ensures feasibility of problem
    end

    UB = Inf; LB = -Inf
    rootStabilizationTrick = :inOut
    stabilizationPoint=theStabilizationPoint
    rootCutsSense=1
    ε= 1e-10
    λ = (rootStabilizationTrick == :inOut) ? .1 : 1.
    δ = (rootStabilizationTrick == :inOut || rootStabilizationTrick == :twoEps) ? 2*ε : 0.
    rootCutsLim=kelleyPrimalEpochs
    rootCutCount = 0
    oaRootCutCount=0
    consecutiveNonImprov_1 = 0
    consecutiveNonImprov_2 = 0
    theCutPool=Cut[]
    for epoch in 1:kelleyPrimalEpochs
      #@show epoch, LB, UB, oaRootCutCount
      solve(rootModel)
      zstar = min.(max.(getvalue(sRoot), 0.), 1.) #cope with numerical errors
      stabilizationPoint += zstar; stabilizationPoint /= 2

      if LB >= getobjectivevalue(rootModel) - eps()
        if consecutiveNonImprov_1 == 5
          consecutiveNonImprov_2 += 1
        else
          consecutiveNonImprov_1 += 1
        end
      else
        if consecutiveNonImprov_1 < 5
          consecutiveNonImprov_1 = 0
        end
        consecutiveNonImprov_2 = 0
      end
      LB = max(LB, getobjectivevalue(rootModel))

      if consecutiveNonImprov_1 == 5
        λ = 1
      elseif consecutiveNonImprov_2 == 5
        δ = 0.
      end

      z0 = λ*zstar + (1-λ)*stabilizationPoint .+ δ
      z0=min.(ones(thePortfolio.n), max.(z0, zeros(thePortfolio.n)))

      Cut  = portfolios_objective2(thePortfolio, z0)
      addOACut=true
      if Cut.status==:Optimal
        UB=min(UB, Cut.p)
      else
        z0 = zstar
        z0=min.(ones(thePortfolio.n), max.(z0, zeros(thePortfolio.n)))
        @suppress Cut  = portfolios_objective2(thePortfolio, z0)
        addOACut=false
      end


      if Cut.status==:Optimal
        @constraint(rootModel, tRoot >= Cut.p + dot(Cut.∇s,sRoot-z0))
        if ( (rootCutsSense > 0)&(epoch <= rootCutsLim) )||( (rootCutsSense < 0)&(epoch >= kelleyPrimalEpochs-rootCutsLim+1) ) & addOACut
            Cut.p+=dot(Cut.∇s,-z0)
            push!(theCutPool, Cut)
        end
      else
        @constraint(rootModel, sum(z0[i]*(1.0-sRoot[i])+sRoot[i]*(1.0-z0[i]) for i=1:thePortfolio.n)>=1.0)
        stabilizationPoint=theStabilizationPoint
      end

      if abs(UB-LB)/abs(UB) <= 1e-4 || consecutiveNonImprov_2 >= 10
        #println("Root node gap: ", abs(UB-LB))
      end
    end
    #println("Kelley Primal: Generated $kelleyPrimalEpochs cuts")

  return theCutPool
end
