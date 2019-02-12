include("instances_io.jl")

using JuMP
using GLPKMathProgInterface

function heuri!(inst::Data)
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
            ### Choose random ark
            #(u,v) = sature[rand(1:length(sature))]
            
            ### Choose highest c
            #sort!(sature, by=x->-inst.c[x[1]][x[2]]) # tri dans l'ordre croissant
            #(u,v) = sature[1]

            ### Choose highest d
            sort!(sature, by=x->-inst.d[x[1]][x[2]]) # tri dans l'ordre croissant
            (u,v) = sature[1]

            ### Choose lowest c+d
            #sort!(sature, by=x->inst.c[x[1]][x[2]] + inst.d[x[1]][x[2]]) # tri dans l'ordre croissant
            #(u,v) = sature[1]

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