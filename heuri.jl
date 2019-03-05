include("instances_io.jl")

using JuMP
using GLPKMathProgInterface

# Modes : * 1 = random
#         * 2 = Plus grande distance
#         * 3 = Plus grande pénalité
#         * 4 = Plus grande distance + pénalité
#         * 5 = Plus petite distance
#         * 6 = Plus petite pénalité
#         * 7 = Plus petite distance + pénalité
modes = [1,2,3,4,5,6,7]


function heuristicSolve!(file)
    inf = 0
    sup = Inf
    best_inf = -1
    best_sup = -1
    for mode in modes
        inst = generate(file)
        m_inf, m_sup = heuri!(inst, mode)
        inf = max(inf, m_inf)
        if inf == m_inf
            best_inf = mode
        end
        sup = min(sup, m_sup)
        if sup == m_sup
            best_sup = mode
        end
    end
    return inf, sup, best_inf, best_sup
end

function heuri!(inst::Data, mode::Int)
    A_star = []
    obj = -1
    arks = [(i,j) for i in 1:inst.n for j in 1:inst.n if inst.adj[i][j]]
    while length(A_star) <= inst.k
        #println("\n Itération ", length(A_star))
        m = Model(solver=GLPKSolverMIP())

        @variable(m, y[(i,j) in arks] >= 0)
        @variable(m, z >= 0)

        @constraint(m, penalty[(i,j) in arks], inst.d[i][j] * y[(i,j)] <= z)
        @constraint(m, flow[v in 1:inst.n], sum(y[(u,v)] for u in 1:inst.n if (u,v) in arks) - sum(y[(v,u)] for u in 1:inst.n if(v,u) in arks) ==  (v == inst.n ? 1 : v == 1 ? -1 : 0)) 

        @objective(m, Min, sum(inst.c[u][v] * y[(u,v)] for (u,v) in arks) + inst.k * z)
    
        solve(m)

        if length(A_star) < inst.k
            sature = Array{Tuple{Int64,Int64}}([])
            for (u,v) in arks
                if (inst.d[u][v] * getvalue(y[(u,v)]) - getvalue(z))^2 <= 0.0000001
                    push!(sature, (u,v))
                end
            end
            if mode == 1
            ### Choose random ark
                (u,v) = sature[rand(1:length(sature))]
            
            elseif mode == 2
            ### Choose highest c
                sort!(sature, by=x->-inst.c[x[1]][x[2]]) # tri dans l'ordre croissant
                (u,v) = sature[1]

            elseif mode == 3
            ### Choose highest d
                sort!(sature, by=x->-inst.d[x[1]][x[2]]) # tri dans l'ordre croissant
                (u,v) = sature[1]

            elseif mode == 4
            ### Choose highest c+d
                sort!(sature, by=x->-(inst.c[x[1]][x[2]] + inst.d[x[1]][x[2]])) # tri dans l'ordre croissant
                (u,v) = sature[1]

            elseif mode == 5
            ### Choose lowest c
                sort!(sature, by=x->inst.c[x[1]][x[2]]) # tri dans l'ordre croissant
                (u,v) = sature[1]

            elseif mode == 6
            ### Choose lowest d
                sort!(sature, by=x->inst.d[x[1]][x[2]]) # tri dans l'ordre croissant
                (u,v) = sature[1]

            elseif mode == 7
            ### Choose lowest c + d
                sort!(sature, by=x->inst.c[x[1]][x[2]] + inst.d[x[1]][x[2]]) # tri dans l'ordre croissant
                (u,v) = sature[1]
            end

            push!(A_star, (u,v))
            #println(u, ",", v)
            inst.c[u][v] += inst.d[u][v]
            inst.d[u][v] = 0
        else
            obj = getobjectivevalue(m)
            break
        end
    end


    s1 = Model(solver=GLPKSolverMIP())

    @variable(s1, x[(i,j) in arks], Bin)

    @constraint(s1, depart,sum(x[(1,j)] for j in 1:inst.n if (1,j) in arks) == 1)
    @constraint(s1, arrivee,sum(x[(i,inst.n)] for i in 1:inst.n if (i,inst.n) in arks) == 1)
    @constraint(s1, flot[v in 2:inst.n-1], sum(x[(i,v)] for i in 1:inst.n if (i,v) in arks) - sum(x[(v,j)] for j in 1:inst.n if (v,j) in arks) == 0)


    @objective(s1, Min, sum(inst.c[i][j]*x[(i,j)] for (i,j) in arks))

    solve(s1)

    inf = getobjectivevalue(s1)

    #println("Borne inf : ", inf)
    #println("Borne sup : ", obj)
    #println("Gap :       ", 100*(obj - getobjectivevalue(s1))/obj, "%")

    return inf, obj

end