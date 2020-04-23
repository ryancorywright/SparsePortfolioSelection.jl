include("src/core.jl")
using JLD
function createWilshire5000Data(rank, theK)
    wilshire5000Data=load("data/Wilshire5000Rank"*string(rank)*".jld", "Wilshire5000Rank"*string(rank))

    n=size(wilshire5000Data.X,2)
    μ=wilshire5000Data.μ*2769/(11*12)
    X=wilshire5000Data.X*sqrt((2769/(11*12)))
    Σ= X'*X

    min_investment=zeros(n)

    γ=(1.0/sqrt(n))*ones(n)

    # Compute r_min and r_max, in order to add minimum return constraint.

    # r_min
    m=Model(solver=MosekSolver(MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0));
    @variable(m, x[1:n]>=0.0);
    @constraint(m, sum(x)==1.0);
    @objective(m, Min,  ((1.0/(2.0*γ[1]))*x'*x +0.5*x'*Σ*x)');
    solve(m)
    r_min=μ'*getvalue(x)
    # r_max
    m2=Model(solver=MosekSolver(MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0))
    @variable(m2, x2[1:n]>=0.0)
    @constraint(m2, sum(x2)==1.0)
    @objective(m2, Max,  (μ'*x2-(1.0/(2.0*γ[1]))*x2'*x2)')
    solve(m2)
    r_max=μ'*getvalue(x2)
    r_bar=r_min+0.3*(r_max-r_min)

    A=μ'
    l=[r_bar]
    u=[1e12]

    μ=zeros(n)
    #########################################################################

    k=theK

    # A=zeros(1,n)
    # l=[0.0]
    # u=[0.0]
    f=size(X, 1)
    m=size(A, 1)
    Y=(X*X')\(X*μ)
    d=X'*Y-μ

    return SparsePortfolioData(μ, Y, d, X, A, l, u, k, n, m, f, min_investment, γ)
end

for theK in [10, 50, 100, 200]
        for rank in [100, 200, 500, 1000]
        thePortfolioData=createWilshire5000Data(rank, theK)

        sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=false, minReturnConstraint=true, useHeuristic=false, useWarmStart=true, useCopyOfVariables=false, usingKelleyPrimal=true, maxCutCallbacks=0);
        @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
        println("Variance is:", sparse_Portfolio.portfolioVariance)
        println("Return is:", sparse_Portfolio.expectedReturn)
        println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
        open("results/Table9ResultsOutput_withminreturn.csv","a") do fp
          println(fp, "Instance, cardinality, regularizer, method, solvetime, boundgap, boundgapSOCP, nodecount, number of iterations")
          println(fp, "rank"*string(rank), ",", thePortfolioData.k, ",", "1/sqrt(n)", ",", "OA", ",",   sparse_Portfolio.solveTime, ",", sparse_Portfolio.boundGap, ",", sparse_Portfolio.boundGapSOCP, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
        end


        sparse_Portfolio_Cplex=cplex_raw_MISOCP(thePortfolioData, 600.0, 1e-4, 1)
        #@show @test(abs.(sum(sparse_Portfolio_Cplex.w)-1.0)<1e-4)
        println("Variance is:", sparse_Portfolio_Cplex.portfolioVariance)
        println("Return is:", sparse_Portfolio_Cplex.expectedReturn)
        println("Reg term is:", 0.5*sum(sparse_Portfolio_Cplex.w[i]^2/thePortfolioData.γ[sparse_Portfolio_Cplex.indices[i]] for i=1:size(sparse_Portfolio_Cplex.indices,1)))

        open("results/Table9ResultsOutput_withminreturn.csv","a") do fp
          println(fp, "Instance, cardinality, regularizer, method, solvetime, boundgap, boundgapSOCP, nodecount, number of iterations")
          println(fp, "rank"*string(rank), ",", thePortfolioData.k, ",",
                            "1/sqrt(n)", ",", "CPLEX MISOCP", ",",  sparse_Portfolio_Cplex.solveTime, ",", sparse_Portfolio_Cplex.boundGap, ",", sparse_Portfolio_Cplex.boundGapSOCP, ",", sparse_Portfolio_Cplex.nodeCount, ",", sparse_Portfolio_Cplex.numIters)
        end
    end
end
