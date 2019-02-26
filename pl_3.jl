include("instances_io.jl")

using JuMP
using Cbc
@doc """
Troisème problème

Ne pas : @variable(m, -1 <) x)
Mais plutôt : @variable(m, x>=-1)

Ne pas Pkg.add('GLPK')
Mais plutôt Pkg.add('GLPKMathProgInterface') (double quote)
@variable(m, x[[[(i,j) for i in 1:10
                       for j in 1:10
                       if i <= j]])

On calcule le plus coourt chemin de s ) t ne passant pas par les arcs uv tels que xuv = 1
On note ici un chemin y indicé par les arêtes ; liste d'arêtes ? Ou matrice booléenne de même taille que adj (homomorphe) ?
""" ->
function shortest_path(graph::Data, x::Array{Array{Bool,1},1})
  m = Model(solver=CbcSolver())
  @variable(m, y[u in 1:graph.n, v in 1:graph.n], Bin)

  @constraint(m, sum(graph.adj[1][v]*y[1,v] for v in 1:graph.n) == 1)
  @constraint(m, sum(graph.adj[v][1]*y[v,1] for v in 1:graph.n) == 0)
  @constraint(m, sum(graph.adj[v][graph.n]*y[v,graph.n] for v in 1:graph.n) == 1)
  @constraint(m, sum(graph.adj[graph.n][v]*y[graph.n,v] for v in 1:graph.n) == 0)

  @constraint(m, cons[u in 2:(graph.n-1)], sum(graph.adj[u][w]*y[u,w] - graph.adj[w][u]*y[w,u] for w in 1:graph.n) == 0)

  @constraint(m, xory[u = 1:graph.n, v = 1:graph.n],  x[u][v] + y[u,v] <= graph.adj[u][v])

  @objective(m, Min, sum(graph.adj[u][v]*y[u,v]*graph.c[u][v] for u in 1:graph.n for v in 1:graph.n))

  status = solve(m, relaxation=false)
  return m, status, y
end








@doc """
A partir de la variable y et de la taille n, renvoie la liste des noeuds du chemin contenu dans y
""" ->
function get_path(y, n)
  path = [1]
  node = 1
  Y = []
  print("Path is : ")
  while node != n
    for i in 1:n
      if getvalue(y[node,i]) == 1
        push!(path, i)
	push!(Y, (node,i))
	print(" ", node, " ")
	node = i
	break
      end
    end
  end
  return Y
end


function get_x(x, n)
  final_x = []
  for u in 1:n
    for v in 1:n
      if getvalue(x[u, v]) == 1
        push!(final_x, (u,v))
      end
    end
  end
  return final_x
end






@doc """
Fonction chapeau : résout le problème master, appelle les sous problèmes
""" ->
function master(graph::Data)
  # Hyp faite : duv = infini
  Y = []
  longest_path_length = -Inf
  x0 = [[false for i in 1:graph.n] for j in 1:graph.n]
  m_0, status, y = shortest_path(graph, x0)
  longest_x = x0
  if status != :Optimal
    println("Error : first path not found :")
    return
  end
  path = get_path(y, graph.n)
  push!(Y, path)
  println("Initialization done.")
  while status == :Optimal
    println("Solving master...")
    m = Model(solver=CbcSolver())

    @variable(m, x[1:graph.n, 1:graph.n], Bin)

    @constraint(m, cons[path in Y], sum(x[u, v] for (u,v) in path) >= 1)
    @constraint(m, sum(x[u, v] for u in 1:graph.n for v in 1:graph.n) <= graph.k)
    @constraint(m, xoradj[u in 1:graph.n, v in 1:graph.n], x[u, v] <= graph.adj[u][v])

    status = solve(m, relaxation=false)
    if status != :Optimal
      println("Master problem unfeasible")
      break
    end
    x_val = [[getvalue(x[u, v]) == 1.0 for v in 1:graph.n] for u in 1:graph.n]
    println("Master solved with status ", status)


    println("Solving slave...")
    model, slave_status, slave_path =  shortest_path(graph, x_val)
    if slave_status != :Optimal
      println("Error : sub problem is unfeasible for x = ", x_val)
      break
    end
    if getobjectivevalue(model) > longest_path_length
      longest_x = x
      longest_path_length = getobjectivevalue(model)
    end
    push!(Y, get_path(slave_path, graph.n))
  end
  print("Final solution is ", get_x(longest_x, graph.n))
end
