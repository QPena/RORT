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

function heuri!(inst::Data, mode::Int)
    A_star = []
    obj = -1
    while length(A_star) <= inst.k
        #println("\n Itération ", length(A_star))
        m = Model(solver=GLPKSolverMIP())

        arks = [(i,j) for i in 1:inst.n for j in 1:inst.n if inst.adj[i][j]]

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
            ### Choose lowest c+d
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

    return obj

end