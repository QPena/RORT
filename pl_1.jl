include("instances_io.jl")

using JuMP
using GLPKMathProgInterface



function problem_1(inst::Data)

    m = Model(solver=GLPKSolverMIP())

    arks = [(i,j) for i in 1:inst.n for j in 1:inst.n if inst.adj[i][j]]


    @variable(m, x[(i,j) in arks], Bin)
    @variable(m, pi[1:inst.n])

    @constraint(m, cons[(i,j) in arks], pi[j] - pi[i] - inst.d[i][j] * x[(i,j)] <= inst.c[i][j])
    @constraint(m, cardinality, sum(x[(i,j)] for (i,j) in arks) <= inst.k)

    @objective(m, Max, pi[inst.n] - pi[1])

    solve(m)
    #print(getobjectivevalue(m))

    return getobjectivevalue(m)

end