#ENV["CPLEX_STUDIO_BINARIES"]="/home/software/cplex/128/cplex/bin/x86-64_linux"
#ENV["MOSEKBINDIR"]="/home/software/mosek/9/tools/platform/linux64x86/bin"
#ENV["MOSEKLM_LICENSE_FILE"]="/home/software/mosek/mosek.lic"
include("src/core.jl")
function getProblemDataFG_minreturn(pathName, pathName2, pathName3, probSize, type, instance, theK)
    # Size: 200, 300, 400
    # Type: Pard
    # Instance: 1-10, converted to equivalent letter.
    alphabet=["a","b","c","d","e","f","g","h","i","j"]

    # Adjust file path appropriately
    if probSize==200
        if type==:Pard
            pathName=pathName*string(probSize)*"/pard"*string(probSize)*"_"*alphabet[instance]
        end
    else
        if type==:Pard
            pathName=pathName*string(probSize)*"/pard"*string(probSize)*"_"*alphabet[instance]
        end
    end

    # Read in .txt file to find n, μ
    returnData=readdlm(pathName*".txt")
    n=returnData[1,1]
    μ=returnData[2:n+1, 1]
    # Read in .rho file to find r_min
    minReturn=readdlm(pathName*".rho")
    r_bar=minReturn[1,1]
    # Read in .bds file to find min_investment, u
    bounds=readdlm(pathName*".bds")
    min_investment=bounds[1:n, 1]
    ub=bounds[1:n, 2]
    # Read in .mat to find sigma
    cov=readdlm(pathName*".mat")
    Σ=cov[2:n+1, 1:n]
    Σ=convert(Array{Float64,2}, Σ)

    # Read in diagonal vector precomputed by Zheng et al
    if probSize==200
        if type==:Pard
            pathName2=pathName2*"/pard"*string(probSize)*"_"*alphabet[instance]
        end
    else
        if type==:Pard
            pathName2=pathName2*"/pard"*string(probSize)*"_"*alphabet[instance]
        end
    end

    if probSize==200
        if type==:Pard
            pathName3=pathName3*"/pard"*string(probSize)*"_"*alphabet[instance]
        end
    else
        if type==:Pard
            pathName3=pathName3*"/pard"*string(probSize)*"_"*alphabet[instance]
        end
    end
    γ=(1000/n)*ones(n)

    # Purpose of these lines was to attempt to take convex combination of diagonal matrices, as proposed in Zheng et. al. (2014). Unfortunately, just using the Frangioni/Gentile is faster
    dominance=readdlm(pathName2*".diag")
    theν=dominance[2:n+1,1]
    dominance=readdlm(pathName3*".diag")
    theν=0.7*theν+0.3*dominance[2:n+1,1]


    for i=1:n
        γ[i]=((1.0/γ[i])+theν[i])^(-1) #adjust γ to have more diagonal dominance, Σ to be smaller
    end
    Σ=Σ-Diagonal(theν)
    X=sqrt(Σ)
    k=theK

    A=vcat(μ', Matrix(1.0I, n,n)) # Can also explicitly create an upper bound on x, if creating a diagonal matrix is too slow/cumbersome.
    l=vcat([r_bar], zeros(n))
    u=vcat([1e12], ub)
    f=size(X, 2)
    m=size(A, 1)
    μ=zeros(n) #Since constraining μ rather than putting in objective
    Y=zeros(f)
    d=zeros(n)

    return SparsePortfolioData(μ, Y, d, X, A, l, u, k, n, m, f, min_investment, γ)
end

for portfolioNumber in [8]

    for theK in [10]
     thePath=pwd()
     pathName = thePath*"/data/MV/size"
     pathName2 = thePath*"/data/diagonals/s/"
     pathName3 = thePath*"/data/diagonals/s/"
     thePortfolioData=getProblemDataFG_minreturn(pathName,pathName2, pathName3, 200, :Pard, portfolioNumber, theK)

     sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=true, minReturnConstraint=true, useHeuristic=false, useWarmStart=false, useCopyOfVariables=true);
     @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
     println("Variance is:", sparse_Portfolio.portfolioVariance)
     println("Return is:", sparse_Portfolio.expectedReturn)
     println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
     open("results/Table6ResultsOutput_noddterm.csv","a") do fp
     println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount, iters")
     println(fp, "pard400_"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                       "1000/(n)", ",", "OA", ",",  sparse_Portfolio.solveTime, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
     end

     sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=true, minReturnConstraint=true, useHeuristic=false, useWarmStart=false, useCopyOfVariables=true, usingKelleyPrimal=true, maxCutCallbacks=0);
     @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
     println("Variance is:", sparse_Portfolio.portfolioVariance)
     println("Return is:", sparse_Portfolio.expectedReturn)
     println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
     open("results/Table6ResultsOutput_noddterm.csv","a") do fp
     println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount, iters")
     println(fp, "pard400_"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                       "1000/(n)", ",", "OA Kelley", ",",  sparse_Portfolio.solveTime, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
     end

     sparse_Portfolio=cutting_planes_portfolios(thePortfolioData; useSOCPLB=true, minReturnConstraint=true, useHeuristic=false, useWarmStart=false, useCopyOfVariables=true, usingKelleyPrimal=true, maxCutCallbacks=200);
     @show @test(abs.(sum(sparse_Portfolio.w)-1.0)<1e-4)
     println("Variance is:", sparse_Portfolio.portfolioVariance)
     println("Return is:", sparse_Portfolio.expectedReturn)
     println("Reg term is:", 0.5*sum(sparse_Portfolio.w[i]^2/thePortfolioData.γ[sparse_Portfolio.indices[i]] for i=1:size(sparse_Portfolio.indices,1)))
     open("results/Table6ResultsOutput_noddterm.csv","a") do fp
     println(fp, "Instance, cardinality, regularizer, method, solvetime, nodecount, iters")
     println(fp, "pard400_"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                       "1000/(n)", ",", "OA Kelley 50", ",",  sparse_Portfolio.solveTime, ",", sparse_Portfolio.nodeCount, ",", sparse_Portfolio.numIters)
     end

    sparse_Portfolio_Cplex=cplex_raw_MISOCP(thePortfolioData, 600.0, 1e-8, portfolioNumber)
    @show @test(abs.(sum(sparse_Portfolio_Cplex.w)-1.0)<1e-4)
    println("Variance is:", sparse_Portfolio_Cplex.portfolioVariance)
    println("Return is:", sparse_Portfolio_Cplex.expectedReturn)
    println("Reg term is:", 0.5*sum(sparse_Portfolio_Cplex.w[i]^2/thePortfolioData.γ[sparse_Portfolio_Cplex.indices[i]] for i=1:size(sparse_Portfolio_Cplex.indices,1)))
    open("results/Table6ResultsOutput_noddterm.csv","a") do fp
    println(fp, "Instance, cardinality, regularizer, method, solvetime")
    println(fp, "pard400_"*string(portfolioNumber), ",", thePortfolioData.k, ",",
                        "1000/sqrt(n)", ",", "CPLEX MISOCP", ",",  sparse_Portfolio_Cplex.solveTime, ",", sparse_Portfolio_Cplex.nodeCount)
        end
    end
end
