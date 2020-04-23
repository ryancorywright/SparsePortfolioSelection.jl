# The below line sometimes needs to be executed on the engaging cluster
#ENV["CPLEX_STUDIO_BINARIES"]="/home/software/cplex/128/cplex/bin/x86-64_linux"

function cutting_planes_portfolios(thePortfolio::SparsePortfolioData, ΔT_max::Float64=600.0, gap::Float64=1e-4;
    mip_solver::MathProgBase.AbstractMathProgSolver=CplexSolver(CPX_PARAM_EPGAP=gap, CPX_PARAM_TILIM=ΔT_max), numRandomRestarts::Int=5, TrackCuts::Int=0, useSOCPLB::Bool=false, minReturnConstraint=false, useHeuristic=true, useWarmStart=true, useCopyOfVariables=false, usingKelleyPrimal::Bool=false, maxCutCallbacks::Int64=0)

  numIters::Int128=0
  bestbound::Float64=0.0

  miop = Model(solver=mip_solver)
  CurrentTime=time()

  # Optimization variables
  @variable(miop, s[1:thePortfolio.n], Bin)
  @variable(miop, t>=-1e12)

  # Objective
  @objective(miop, Min, t)

  # Constraints
  @constraint(miop, sum(s)<=thePortfolio.k)
  @constraint(miop, sum(s)>=1)

  # Min return constraint
  if minReturnConstraint #This line assumes that minimum return constraint comes first whenever it is included
    @constraint(miop, sum(s[i]*(thePortfolio.A[1,i].>thePortfolio.l[1]) for i=1:thePortfolio.n)>=1.0) #ensures feasibility of problem
  end

  # Max investment constraint
  if useCopyOfVariables 
    @variable(miop,   x_master[1:thePortfolio.n]>=0.0)
    @constraint(miop, thePortfolio.A*x_master.>=thePortfolio.l)
    @constraint(miop, thePortfolio.A*x_master.<=thePortfolio.u)
    @constraint(miop, minInvestment[i=1:thePortfolio.n], x_master[i]>=s[i]*thePortfolio.min_investment[i])
    @constraint(miop, bigMMaster[i=1:thePortfolio.n], x_master[i]<=thePortfolio.u[i+1]*s[i]) #todo: code up less adhoc way of enforcing the u_i part.
    @constraint(miop, sumsToOne, sum(x_master[i] for i=1:thePortfolio.n)==1.0)
  end

  # Buy-in constraint to retain feasibility
  @constraint(miop, sum(s[i]*thePortfolio.min_investment[i] for i=1:thePortfolio.n)<=1.0)

  lowerBound=-1e12

  if useSOCPLB
    dual_vars, socplb =portfolios_socp(thePortfolio)
    @show socplb
    lowerBound=socplb
  end

  bests0=zeros(thePortfolio.n)
  theCut=Cut(-1e12, zeros(thePortfolio.n), :Optimal)
  if useWarmStart
    @suppress bests0=getWarmStart(thePortfolio, numRandomRestarts)
    theCut =portfolios_objective(thePortfolio, bests0)
  end

  if sum(bests0)<=thePortfolio.k  && theCut.status==:Optimal && useWarmStart # If we have general constraints, warm-start may inject an infeasible z as a warm-start, which can give incorrect output.
    setvalue(s, bests0)
  end
  if theCut.status==:Optimal && useWarmStart
    @constraint(miop, t>=theCut.p + dot(theCut.∇s, s-bests0))
  end

  bestUB=1e12
  bestCplexUB=1e12
  heuristicUB=1e12
  if useWarmStart
    bestUB=theCut.p
    bestCplexUB=theCut.p
    heuristicUB=theCut.p
  end

  globalbests0=zeros(thePortfolio.n)

  socpgap=(bestUB-lowerBound)/abs(bestUB)
  if useSOCPLB && theCut.status==:Optimal && useWarmStart
    println("***********************")
    println("*** SOCP Gap is: ",round.(socpgap*10000)/100, "% *****")
    println("***********************")
  end
  bcdata = CutIterData[]

  ############################################
  ####  Kelley Cut loop here       ###########
  ############################################
  # Optimization variables

   kelleyCutPool=Cut[]
   if usingKelleyPrimal
     zSOCP=cplex_MISOCP_relaxation(thePortfolio, ΔT_max)
     stabilizationPoint=zSOCP
     kelleyCutPool=getKelleyPrimalCuts(thePortfolio, useCopyOfVariables, minReturnConstraint, stabilizationPoint, 200)
     if length(kelleyCutPool) > 0
       for c in kelleyCutPool
         @constraint(miop, t >= c.p + dot(c.∇s, s))
       end
     end
   end

############################################
#### End Kelley Cut loop here    ###########
############################################

  function outer_approximation(cb)
    numIters=numIters+1

    # if (bestCplexUB-getobjbound(miop))<1e-4 && useSOCPLB # socp bound is tight, so use it.
    #       #@lazyconstraint(cb, t>= socplb)
    # end

    s0 = min.(max.(getvalue(s), 0.), 1.) #Round s to be integral
    s0 = 1.0.*round.(Int, s0)
    @suppress Cut = portfolios_objective(thePortfolio, s0)
    #theStatus, bestx= cplex_raw_MISOCP_subset(thePortfolio, s0) # this should match Cut.p, so can use for debugging
    indices0=findall(s->s>0.5, s0)

    if useHeuristic
      if Cut.p<1.05*bestCplexUB
        indices0 = portfolios_hillclimb(thePortfolio, indices0)[1]
        s00 = zeros(thePortfolio.n)
        s00[indices0].=1
        Cut2 = portfolios_objective(thePortfolio, s00)
        if Cut2.p<bestCplexUB
          globalbests0=copy(s00)
          bestUB=Cut2.p
        end
        @lazyconstraint(cb, t>=Cut2.p +dot(Cut2.∇s, s-s00))
      end
    end
    if Cut.status==:Optimal
        @lazyconstraint(cb, t>=Cut.p +dot(Cut.∇s, s-s0))
    else
      @lazyconstraint(cb, sum(s0[i]*(1.0-s[i])+s[i]*(1.0-s0[i]) for i=1:thePortfolio.n)>=1.0) # cut off this solution, as it is infeasible
      println("Warning: Generated Feasibility Cut!")
    end

    bestCplexUB=min(bestCplexUB, Cut.p)
    bestUB=min(bestUB, Cut.p)
    bestbound = MathProgBase.cbgetbestbound(cb)
    if TrackCuts>0
      push!(bcdata, CutIterData(time()-CurrentTime,numIters,bestUB,bestbound, (bestUB-bestbound)/bestUB, Cut.status))
    end
    Cut=nothing
  end
  addlazycallback(miop, outer_approximation, fractional=false)

  numCutCallbacks=0
  function fractional_cut(cb)
    if numCutCallbacks < maxCutCallbacks
      s0 = getvalue(s)
      stabilizationPoint=s0
      userCutPool=getKelleyPrimalCuts(thePortfolio, useCopyOfVariables, minReturnConstraint, stabilizationPoint, 10)
      if length(userCutPool) > 0
        for c in userCutPool
          @usercut(cb, t >= c.p + dot(c.∇s, s))
        end
      end
      numCutCallbacks+=1
    end
  end
  addcutcallback(miop, fractional_cut)

  solve(miop)

  if TrackCuts>0
    open("branchandcuttrack.csv","a") do fp
        println(fp, "time,cut,obj,bestbound, boundgap")
        for bb in bcdata
            println(fp, bb.time, ",", bb.Cut, ",",
                        bb.obj, ",", bb.bestbound, ",", bb.boundgap)
        end
          println(fp, " , , , , ")
    end
  end

  println("Number of lazy constraints used is:", numIters)
  if abs((heuristicUB-getobjectivevalue(miop))/getobjectivevalue(miop))<gap
    println("Heuristic warm-start was globally optimal")
  else
    println("Cutting-planes improved on the warm-start by ", round(10000*(heuristicUB-getobjectivevalue(miop))/getobjectivevalue(miop))/100, "%")
  end

  indices=findall(s->s>0.5, getvalue(s))

  U = thePortfolio.X[:, indices]
  f2=size(U,2)
  w=fill(0.0, f2)
  dual_vars= inner_dual(thePortfolio, indices)
  w = thePortfolio.γ[indices].*dual_vars.w
  return CardinalityConstrainedPortfolio(indices, w, dual_vars.λ, dual_vars.α, dual_vars.βl, dual_vars.βu, dual_vars.ρ, (getobjectivevalue(miop)-getobjbound(miop))/abs(getobjectivevalue(miop)), (getobjectivevalue(miop)-lowerBound)/abs(getobjectivevalue(miop)), numIters, (thePortfolio.μ[indices]'*w)[1], (w'*U'*U*w)[1], getsolvetime(miop), getnodecount(miop))
end

function portfolios_objective(thePortfolio::SparsePortfolioData, s::Array{Float64,1})

  indices = findall(s->s>0.5, s)

  k = size(indices,1)
  dual_vars=inner_dual(thePortfolio, indices)
  p=dual_vars.ofv
  ρ_full=zeros(thePortfolio.n)
  ρ_full[indices].=dual_vars.ρ
  w= thePortfolio.X'*dual_vars.α+dual_vars.λ*fill(1.0, thePortfolio.n)+thePortfolio.A'*(dual_vars.βl-dual_vars.βu)+ρ_full-thePortfolio.d
  for i=1:thePortfolio.n
    if abs(s[i])<1e-8
      if w[i]<=-1e-8
        w[i]=0.0
      end
    end
  end
  ∇s = -thePortfolio.γ.*(w).^2/(2.0)+ ρ_full.*thePortfolio.min_investment

  return Cut(p[1], ∇s, dual_vars.status)
end

function portfolios_objective2(thePortfolio::SparsePortfolioData, s0::Array{Float64,1})

  indices = findall(s->s>1e-6, s0)

  k = size(indices,1)
  dual_vars=inner_dual2(thePortfolio, indices, s0)
  p=dual_vars.ofv
  ρ_full=zeros(thePortfolio.n)
  ρ_full[indices].=dual_vars.ρ
  w= thePortfolio.X'*dual_vars.α+dual_vars.λ*fill(1.0, thePortfolio.n)+thePortfolio.A'*(dual_vars.βl-dual_vars.βu)+ρ_full-thePortfolio.d
  for i=1:thePortfolio.n
    if abs(s0[i])<1e-8 # Pareto-optimal cut selection
      if w[i]<=-1e-8
        w[i]=0.0
      end
    end
  end
  ∇s = -thePortfolio.γ.*(w).^2/(2.0)+ ρ_full.*thePortfolio.min_investment

  return Cut(p[1], ∇s, dual_vars.status)
end

function getWarmStart(thePortfolio::SparsePortfolioData, numRandomRestarts::Int)
    bestp0=1e14
    bests0=zeros(thePortfolio.n)
    for i=1:numRandomRestarts
      indices0 = portfolios_hillclimb(thePortfolio, sort(shuffle(1:thePortfolio.n)[1:thePortfolio.k]))[1] #start at random indices
      while sum(thePortfolio.min_investment[indices0])>1.0 #iterate is infeasible, so cleaning solution.
        index=rand(indices0)
        indices0=filter(e->e≠index, indices0)
      end
      s0 = zeros(thePortfolio.n)
      s0[indices0].=1
      Cut = portfolios_objective(thePortfolio, s0)
      #@show cplex_raw_MISOCP_subset(thePortfolio, s0) #good for debugging
      if Cut.p<bestp0
        bestp0=Cut.p
        bestIndices0=indices0
        bests0=zeros(thePortfolio.n)
        bests0[bestIndices0].=1
      end
    end
    return bests0
end
