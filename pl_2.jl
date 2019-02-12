include("instances_io.jl")

using JuMP
using GLPKMathProgInterface
inst=generate(10,10,5,20,20)


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
print("initialisation de Y")
println(getvalue(x))
println(getobjectivevalue(s1))
Y0=[]
X0=getvalue(x)
for (i,j) in arks
	if X0[(i,j)]==1
		push!(Y0,(i,j))
	end
end
push!(Y,Y0)

Z=1000
valeur=0
K=1
X=[]

m = Model(solver=GLPKSolverMIP())
	
@variable(m, y[(i,j) in arks], Bin)

@variable(m, z>=0)

@constraint(m, cardinality, sum(y[(i,j)] for (i,j) in arks) <= inst.k)
@objective(m, Max, z)

while Z>valeur

#problème maître
	
    @constraint(m ,sum(inst.c[i][j] + inst.d[i][j]*y[(i,j)] for (i,j) in Y[K]) - z >= 0)


	

	solve(m)
	X=getvalue(y)
	Z=getobjectivevalue(m)
	println("résolution problème maître")
	println("Z= ",Z)

	#s = Model(solver=GLPKSolverMIP())


	#@variable(s, x[(i,j) in arks], Bin)

#	@constraint(s, depart,sum(x[(1,j)] for j in 1:inst.n if (1,j) in arks) == 1)
#	@constraint(s, arrivee,sum(x[(i,inst.n)] for i in 1:inst.n if (i,inst.n) in arks) == 1)
#	@constraint(s, flot[v in 2:inst.n-1], sum(x[(i,v)] for i in 1:inst.n if (i,v) in arks) - sum(x[(v,j)] for j in 1:inst.n if (v,j) in arks) == 0)
    
	@objective(s1, Min, sum((inst.c[i][j] + inst.d[i][j]*X[(i,j)])*x[(i,j)] for (i,j) in arks))

	solve(s1)
	valeur=getobjectivevalue(s1)
    println("valeur problème esclave= ",valeur)


	println("résolution problème esclave")

	if valeur<Z
		Yi=[]
		Xi=getvalue(x)
		for (i,j) in arks
			if Xi[(i,j)]==1
				push!(Yi,(i,j))
			end

		end
        push!(Y,Yi)
 #       println("Yi= ",Yi)

		global K+=1
#		println("itération ",K)
#        println(Y[K])
	else 
		println("z= ",Z,"x= ",X)
		break
			
	end
		

end

