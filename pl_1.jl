include("instances_io.jl")

using JuMP
using GLPKMathProgInterface


@doc """
Premier problème

On résout directement un programme linéaire représentant le problème d'interdictions de plus courts chemins
"""
function problem_1(inst::Data)

    m = Model(solver=GLPKSolverMIP())

    # On récupère la liste des arcs
    arks = [(i,j) for i in 1:inst.n for j in 1:inst.n if inst.adj[i][j]]

    # On crée les variables associées à chaque arc et à chaque noeud
    @variable(m, x[(i,j) in arks], Bin)
    @variable(m, pi[1:inst.n])

    # On ajoute les contraintes
    @constraint(m, cons[(i,j) in arks], pi[j] - pi[i] - inst.d[i][j] * x[(i,j)] <= inst.c[i][j])
    # Contrainte de cardinalité
    @constraint(m, cardinality, sum(x[(i,j)] for (i,j) in arks) <= inst.k)

    # Fonction objectif
    @objective(m, Max, pi[inst.n] - pi[1])

    # On résout le PL
    solve(m)

    # On renvoie la valeur du PL
    return getobjectivevalue(m)

end