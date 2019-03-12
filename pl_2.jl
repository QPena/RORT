include("instances_io.jl")

using JuMP
using GLPKMathProgInterface
#inst=generate(5,10,7,20,20)

function problem_2(inst::Data)
	Y=[]

	arks = [(i,j) for i in 1:inst.n for j in 1:inst.n if inst.adj[i][j]]

	#initialisation de Y par calcul du pcc
	s1 = Model(solver=GLPKSolverMIP())

	@variable(s1, x[(i,j) in arks], Bin)

	@constraint(s1, depart,sum(x[(1,j)] for j in 1:inst.n if (1,j) in arks) == 1)
	@constraint(s1, arrivee,sum(x[(i,inst.n)] for i in 1:inst.n if (i,inst.n) in arks) == 1)
	@constraint(s1, flot[v in 2:inst.n-1], sum(x[(i,v)] for i in 1:inst.n if (i,v) in arks) - sum(x[(v,j)] for j in 1:inst.n if (v,j) in arks) == 0)


	@objective(s1, Min, sum(inst.c[i][j]*x[(i,j)] for (i,j) in arks))

	solve(s1)

	Y0=[]
	X0=getvalue(x)
	for (i,j) in arks
		if X0[(i,j)]==1
			push!(Y0,(i,j))
		end
	end
	#on rajoute le pcc calculé dans Y
	push!(Y,Y0)

	#initialisation de la valeur du problême maître
	Z=1000
	#initialisation de la valeur du plus court chemin calculée par le problême esclave
	valeur=0
	#initialisation du nombre de coupes
	K=0
	X=[]

	#initialisation du problême maître
	m = Model(solver=GLPKSolverMIP())
		
	@variable(m, y[(i,j) in arks], Bin)

	@variable(m, z>=0)

	@constraint(m, cardinality, sum(y[(i,j)] for (i,j) in arks) <= inst.k)
	@objective(m, Max, z)

	while Z>valeur
	#tant que la valeur du problême maître est strictement supérieure à la valeur du plus court chemin, la boucle continue
	#problème maître
		
		#mise à jour du problême maître par le rajout de la coupe calculée lors de la dernière itération
	    @constraint(m ,sum(inst.c[i][j] + inst.d[i][j]*y[(i,j)] for (i,j) in Y[K]) - z >= 0)

	    #On résout le problême maître
		solve(m)
		X=getvalue(y)
		#mise à jour de Z
		Z=getobjectivevalue(m)

		#mise à jour du problême esclave avec les arcs pénalisés par le problême maître
		@objective(s1, Min, sum((inst.c[i][j] + inst.d[i][j]*X[(i,j)])*x[(i,j)] for (i,j) in arks))
		#On résout le problême esclave
		solve(s1)
		valeur=getobjectivevalue(s1)

		if valeur < Z
			Yi=[]
			Xi=getvalue(x)
			for (i,j) in arks
				if Xi[(i,j)]==1
					push!(Yi,(i,j))
				end

			end
	        push!(Y,Yi)
	 #on rajoute le plus court chemin calculé par le problême maître dans Y
	 		#K=nombre de coupes
			K+=1

		else 
		#le problême est résolu
	        #on retourne la valeur du problème
	        return Z
			break
				
		end
			

	end
end
