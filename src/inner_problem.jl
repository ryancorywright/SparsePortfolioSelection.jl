function inner_dual(thePortfolio::SparsePortfolioData, indices::Array{Int64})


    n=thePortfolio.f
    m=thePortfolio.m
    f=size(indices,1)

    innerQP = Model(solver=MosekSolver(MSK_DPAR_INTPNT_QO_TOL_PFEAS=1e-6, MSK_DPAR_INTPNT_QO_TOL_DFEAS=1e-6, MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0)) #
    #

    @variable(innerQP, α[1:n])
    @variable(innerQP, λ)
    @variable(innerQP, w[1:f])
    @variable(innerQP, ρ[1:f]>=0.0)
    @variable(innerQP, βl[1:m]>=0.0)
    @variable(innerQP, βu[1:m]>=0.0)

    @constraint(innerQP, w.>=(thePortfolio.X[:, indices]'*α+thePortfolio.A[:, indices]'*(βl-βu)+λ*ones(f)+ρ-thePortfolio.d[indices]))
    @objective(innerQP, Max, (-0.5*α'*α-sum((thePortfolio.γ[indices]/2.0).*w.*w)+thePortfolio.Y'*α+λ+βl'*thePortfolio.l-βu'*thePortfolio.u+ρ'*thePortfolio.min_investment[indices])')

    status=solve(innerQP; suppress_warnings=true)
    return Dual(getvalue(α), getvalue(λ), getvalue(βl), getvalue(βu), getvalue(ρ), getvalue(w), getobjectivevalue(innerQP), status)
end

function inner_dual2(thePortfolio::SparsePortfolioData, indices::Array{Int64}, s::Array{Float64})

    n=thePortfolio.f
    m=thePortfolio.m
    f=size(indices,1)

    innerQP = Model(solver=MosekSolver(MSK_DPAR_INTPNT_QO_TOL_PFEAS=1e-6, MSK_DPAR_INTPNT_QO_TOL_DFEAS=1e-6, MSK_IPAR_LOG=0, MSK_IPAR_MAX_NUM_WARNINGS=0)) #
    #

    @variable(innerQP, α[1:n])
    @variable(innerQP, λ)
    @variable(innerQP, w[1:f])
    @variable(innerQP, ρ[1:f]>=0.0)
    @variable(innerQP, βl[1:m]>=0.0)
    @variable(innerQP, βu[1:m]>=0.0)

    @constraint(innerQP, w.>=(thePortfolio.X[:, indices]'*α+thePortfolio.A[:, indices]'*(βl-βu)+λ*ones(f)+ρ-thePortfolio.d[indices]))
    @objective(innerQP, Max, (-0.5*α'*α-sum((thePortfolio.γ[indices]/2.0).*s[indices].*w.*w)+thePortfolio.Y'*α+λ+βl'*thePortfolio.l-βu'*thePortfolio.u+ρ'*(s[indices].*thePortfolio.min_investment[indices]))')

    status=solve(innerQP; suppress_warnings=true)
    return Dual(getvalue(α), getvalue(λ), getvalue(βl), getvalue(βu), getvalue(ρ), getvalue(w), getobjectivevalue(innerQP), status)
end
