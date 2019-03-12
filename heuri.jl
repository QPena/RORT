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


@doc """
Fonction chapeau : résout le problème contenu dans le fichier file avec les différentes heuristiques de choix
""" ->
function heuristicSolve(file)
    inf = 0 # meilleure borne inférieure
    sup = Inf # meilleure borne supérieure
    best_inf = -1 # mode qui donne la meilleure borne inférieure
    best_sup = -1 # mode qui donne la meilleure borne supérieure
    # On résout l'heuristique pour chaque heuristique de choix
    for mode in modes
        # On lit l'instance que l'on veut résoudre
        # Note : on la re-génère à chaque fois, car la fonction de résolution modifie l'instance
        # On pourrait aussi effectuer une deepcopy de l'instance initiale (mais pas sûr que ce soit plus rapide)
        inst = generate(file)
        # On résout l'instance avec le mode sélectionné
        m_inf, m_sup = heuri!(inst, mode)
        # On met à jour la meilleure borne inférieure
        inf = max(inf, m_inf)
        if inf == m_inf
            best_inf = mode
        end
        # On met à jour la meilleure borne inférieure
        sup = min(sup, m_sup)
        if sup == m_sup
            best_sup = mode
        end
    end
    return inf, sup, best_inf, best_sup
end

@doc """
Fonction résolution : résout l'instance inst avec le mode sélectionné
""" ->
function heuri!(inst::Data, mode::Int)
    A_star = [] # les arcs sélectionnés
    sup = -1
    # Liste des arcs du graphe
    arks = [(i,j) for i in 1:inst.n for j in 1:inst.n if inst.adj[i][j]]
    while length(A_star) <= inst.k
        # Création du modèle
        m = Model(solver=GLPKSolverMIP())

        # Variables du problème dual
        @variable(m, y[(i,j) in arks] >= 0)
        @variable(m, z >= 0)

        # Contrainte de pénalité
        @constraint(m, penalty[(i,j) in arks], inst.d[i][j] * y[(i,j)] <= z)
        # Contraintes de flot
        @constraint(m, flow[v in 1:inst.n], sum(y[(u,v)] for u in 1:inst.n if (u,v) in arks) - sum(y[(v,u)] for u in 1:inst.n if(v,u) in arks) ==  (v == inst.n ? 1 : v == 1 ? -1 : 0)) 

        # Objectif
        @objective(m, Min, sum(inst.c[u][v] * y[(u,v)] for (u,v) in arks) + inst.k * z)
    
        # On résout le dual
        solve(m)

        # Si l'on a pas déjà sélectionné k arcs
        if length(A_star) < inst.k
            # On récupère tous les arcs qui saturent la contrainte de pénalité
            sature = Array{Tuple{Int64,Int64}}([])
            for (u,v) in arks
                # Note : le codage des flottants amènent parfois des écarts infimes,
                # on utilise donc ce petit stratagème pour s'assurer de bien récupérer tous les arcs
                if (inst.d[u][v] * getvalue(y[(u,v)]) - getvalue(z))^2 <= 0.0000001
                    push!(sature, (u,v))
                end
            end
            # On sélectionne un arc selon le mode retenu
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

            # On ajoute l'arc sélectionné à la liste des arcs
            push!(A_star, (u,v))

            # On met à jour la distance de l'arc sélectionné
            inst.c[u][v] += inst.d[u][v]
            inst.d[u][v] = 0
        # Si on a déjà sélectionné k arcs, la valeur du dual est notre borne supérieure
        else
            sup = getobjectivevalue(m)
            break
        end
    end

    # On résout le programme linéaire permettant de déterminer le plus court chemin pénalisant les arcs sélectionnés
    # Les distances des arcs pénalisés ont déjà été mises à jour avec la pénalité
    s1 = Model(solver=GLPKSolverMIP())

    # Variables binaires si l'arc (i,j) est dans le plus court chemin
    @variable(s1, x[(i,j) in arks], Bin)

    # Contraintes de flot
    @constraint(s1, depart,sum(x[(1,j)] for j in 1:inst.n if (1,j) in arks) == 1)
    @constraint(s1, arrivee,sum(x[(i,inst.n)] for i in 1:inst.n if (i,inst.n) in arks) == 1)
    @constraint(s1, flot[v in 2:inst.n-1], sum(x[(i,v)] for i in 1:inst.n if (i,v) in arks) - sum(x[(v,j)] for j in 1:inst.n if (v,j) in arks) == 0)

    # Objectif : plus court chemin
    @objective(s1, Min, sum(inst.c[i][j]*x[(i,j)] for (i,j) in arks))

    # On résout le programme linéaire
    solve(s1)
    # On obtient la borne inférieure du problème
    inf = getobjectivevalue(s1)


    return inf, sup

end