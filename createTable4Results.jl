#ENV["CPLEX_STUDIO_BINARIES"]="/home/software/cplex/128/cplex/bin/x86-64_linux"
using DelimitedFiles, LinearAlgebra, Random
include("src/core.jl")
function getProblemDataBeasley_minreturn(pathName)
    testData=readdlm(pathName)
    n=testData[1,1]
    μ=testData[2:n+1, 1]
    μ=convert(Array{Float64}, μ)
    σ=testData[2:n+1, 2]
    Σ=zeros(n,n)
    for iter=n+2:size(testData,1)
        # This loop inserts diagonal elements twice. This is ok as benign, and a "fix" would require more lines of code.
        Σ[testData[iter, 1], testData[iter, 2]]=σ[testData[iter, 1]]*σ[testData[iter, 2]]*testData[iter, 3]
        Σ[testData[iter, 2], testData[iter, 1]]=σ[testData[iter, 1]]*σ[testData[iter, 2]]*testData[iter, 3]
    end
    X=sqrt(Σ)
    k=5 #you will need to manually change this to 10 or 20 to produce the other rows in the table
    γ=(100/sqrt(n))*ones(n)

    # Compute r_min and r_max, in order to add minimum return constraint.

    # r_min
    m=Model(solver=MosekSolver(MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0));
    @variable(m, x[1:n]>=0.0);
    @constraint(m, sum(x)==1.0);
    @objective(m, Min,  ((1.0/(2.0*γ[1]))*x'*x +0.5*x'*Σ*x)'); #Each value of \gamma is the same here, so ok
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
    min_investment=zeros(n)
    f=size(X, 2)
    m=size(A, 1)
    # Constraining minimum return instead, so these entries all become 0.
    μ=zeros(n)
    Y=(X*X')\(X*μ)
    d=X'*Y-μ
    d=round.(d, digits=8) #avoids round-off errors when Σ is full rank

    return SparsePortfolioData(μ, Y, d, X, A, l, u, k, n, m, f, min_investment, γ)
end


for portfolioNumber in 1:5
     pathName = "/Users/Ryancw2/Dropbox (MIT)/Sparse Portfolios Project/OptimalMarkowitzPortfolios.jl/data/port"
     pathName=pathName*string(portfolioNumber)*".txt"
     thePortfolioData=getProblemDataBeasley_minreturn(pathName)

     sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=false, minReturnConstraint=true, useHeuristic=false);
     @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
     println("Variance is:", sparse_Portfolio.portfolioVariance)
     println("Return is:", sparse_Portfolio.expectedReturn)
     println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
     open("results/Table4ResultsOutput.csv","a") do fp
     println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount, iters")
     println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                       "100/sqrt(n)", ",", "OA", ",",  sparse_Portfolio.solveTime, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
     end

     sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=false, minReturnConstraint=true, useHeuristic=false, usingKelleyPrimal=true, maxCutCallbacks=0);
     @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
     println("Variance is:", sparse_Portfolio.portfolioVariance)
     println("Return is:", sparse_Portfolio.expectedReturn)
     println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
     open("results/Table4ResultsOutput_b.csv","a") do fp
     println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount, iters")
     println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                       "100/sqrt(n)", ",", "OA-Kelley (50 cuts)", ",",  sparse_Portfolio.solveTime, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
     end

     # Test against Big-M for correctness
     sparse_Portfolio_Cplex_bigM=cplex_raw_bigM(thePortfolioData, 3600.0, 1e-8, 1)
     @show @test(abs.(sum(sparse_Portfolio_Cplex_bigM.w)-1.0)<1e-4)
     println("Variance is:", sparse_Portfolio_Cplex_bigM.portfolioVariance)
     println("Return is:", sparse_Portfolio_Cplex_bigM.expectedReturn)
     println("Reg term is:", 0.5*sum(sparse_Portfolio_Cplex_bigM.w[i]^2/thePortfolioData.γ[sparse_Portfolio_Cplex_bigM.indices[i]] for i=1:size(sparse_Portfolio_Cplex_bigM.indices,1)))
     open("results/Table4ResultsOutput.csv","a") do fp
     println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount")
     println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                     "100/sqrt(n)", ",", "CPLEX Big-M", ",",  sparse_Portfolio_Cplex_bigM.solveTime, ",", sparse_Portfolio_Cplex_bigM.nodeCount)
    end

    # # #  # Test against MISOCP for correctness
    sparse_Portfolio_Cplex=cplex_raw_MISOCP(thePortfolioData, 3600.0, 1e-8, 1)
    @show @test(abs.(sum(sparse_Portfolio_Cplex.w)-1.0)<1e-4)
    println("Variance is:", sparse_Portfolio_Cplex.portfolioVariance)
    println("Return is:", sparse_Portfolio_Cplex.expectedReturn)
    println("Reg term is:", 0.5*sum(sparse_Portfolio_Cplex.w[i]^2/thePortfolioData.γ[sparse_Portfolio_Cplex.indices[i]] for i=1:size(sparse_Portfolio_Cplex.indices,1)))
    open("results/Table4ResultsOutput.csv","a") do fp
    println(fp, "Instance, cardinality, regularizer, method, solvetime")
    println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                        "100/sqrt(n)", ",", "CPLEX MISOCP", ",",  sparse_Portfolio_Cplex.solveTime, ",", sparse_Portfolio_Cplex.nodeCount)
    end
end
