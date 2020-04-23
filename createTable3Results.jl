#ENV["CPLEX_STUDIO_BINARIES"]="/home/software/cplex/128/cplex/bin/x86-64_linux"
using DelimitedFiles, LinearAlgebra
include("src/core.jl")

function getProblemDataBeasley(pathName)
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
    Y=(X*X')\(X*μ)
    d=X'*Y-μ
    d=round.(d, digits=8) #avoids round-off errors when Σ is full rank
    A=zeros(1,size(X,2))
    l=zeros(size(A,1))
    u=zeros(size(A,1))
    min_investment=zeros(n)
    f=size(X, 2)
    m=size(A, 1)
    k=5
    γ=(100/sqrt(n))*ones(n)

    return SparsePortfolioData(μ, Y, d, X, A, l, u, k, n, m, f, min_investment, γ)
end


for portfolioNumber=1:5
    pathName = "/Users/Ryancw2/Dropbox (MIT)/Sparse Portfolios Project/OptimalMarkowitzPortfolios.jl/data/port"
    pathName=pathName*string(portfolioNumber)*".txt"
    thePortfolioData=getProblemDataBeasley(pathName)

    sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=false, useHeuristic=false);
    @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
    println("Variance is:", sparse_Portfolio.portfolioVariance)
    println("Return is:", sparse_Portfolio.expectedReturn)
    println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
    open("results/Table3ResultsOutput.csv","a") do fp
    println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount, iters")
    println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                     "100/sqrt(n)", ",", "OA", ",",  sparse_Portfolio.solveTime, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
    end


    # # Test against Big-M for correctness
    sparse_Portfolio_Cplex_bigM=cplex_raw_bigM(thePortfolioData, 300.0, 1e-8, 1)
    @show @test(abs.(sum(sparse_Portfolio_Cplex_bigM.w)-1.0)<1e-4)
    println("Variance is:", sparse_Portfolio_Cplex_bigM.portfolioVariance)
    println("Return is:", sparse_Portfolio_Cplex_bigM.expectedReturn)
    println("Reg term is:", 0.5*sum(sparse_Portfolio_Cplex_bigM.w[i]^2/thePortfolioData.γ[sparse_Portfolio_Cplex_bigM.indices[i]] for i=1:size(sparse_Portfolio_Cplex_bigM.indices,1)))
    open("results/Table3ResultsOutput.csv","a") do fp
    println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount")
    println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                     "100/sqrt(n)", ",", "CPLEX Big-M", ",",  sparse_Portfolio_Cplex_bigM.solveTime, ",", sparse_Portfolio_Cplex_bigM.nodeCount)
    end

    # # Test against MISOCP for correctness
    sparse_Portfolio_Cplex=cplex_raw_MISOCP(thePortfolioData, 300.0, 1e-8, 1)
    @show @test(abs.(sum(sparse_Portfolio_Cplex.w)-1.0)<1e-4)
    println("Variance is:", sparse_Portfolio_Cplex.portfolioVariance)
    println("Return is:", sparse_Portfolio_Cplex.expectedReturn)
    println("Reg term is:", 0.5*sum(sparse_Portfolio_Cplex.w[i]^2/thePortfolioData.γ[sparse_Portfolio_Cplex.indices[i]] for i=1:size(sparse_Portfolio_Cplex.indices,1)))
    open("results/Table3ResultsOutput.csv","a") do fp
    println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount")
    println(fp, "beasley"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                    "100/sqrt(n)", ",", "CPLEX MISOCP", ",",  sparse_Portfolio_Cplex.solveTime, ",", sparse_Portfolio_Cplex.nodeCount)
    end
end
