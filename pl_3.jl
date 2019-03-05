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
function shortest_path(graph::Data, x, arks)
  m = Model(solver=CbcSolver())
  @variable(m, y[(u,v) in arks], Bin)

  @constraint(m, sum(y[(u,v)] for (u,v) in arks if u == 1) == 1)
  @constraint(m, sum(y[(u,v)] for (u,v) in arks if v == 1) == 0)
  @constraint(m, sum(y[(u,v)] for (u,v) in arks if v == graph.n) == 1)
  @constraint(m, sum(y[(u,v)] for (u,v) in arks if u == graph.n) == 0)

  @constraint(m, cons[u in 2:(graph.n-1)], sum(y[(u,w)] for w in 1:graph.n if (u,w) in arks) - sum(y[(w,u)] for w in 1:graph.n if (w,u) in arks)== 0)

  @constraint(m, xory[(u,v) in arks],  x[(u,v)] + y[(u,v)] <= 1)

  @objective(m, Min, sum(y[(u,v)]*graph.c[u][v] for (u,v) in arks))

  status = solve(m, relaxation=false)
  return m, y
end








@doc """
A partir de la variable y et de la taille n, renvoie la liste des noeuds du chemin contenu dans y
""" ->
function get_path(y, n, arks)
  path = [1]
  node = 1
  Y = []
  print("Path is : ")
  while node != n
    for (u,v) in arks
      if u == node && getvalue(y[(u,v)]) == 1
        push!(path, u)
	push!(Y, (u,v))
	print(" ", node, " ")
	node = v
	break
      end
    end
  end
  return Y
end


function get_x(x, n, arks)
  final_x = []
  for (u,v) in arks
    if getvalue(x[(u, v)]) == 1
      push!(final_x, (u,v))
    end
  end
  return final_x
end

function get_x_from_value(x, n, arks)
  final_x = []
  for (u,v) in arks
    if x[(u, v)] == 1
      push!(final_x, (u,v))
    end
  end
  return final_x
end





@doc """
Fonction chapeau : résout le problème master, appelle les sous problèmes
""" ->
function master(graph::Data)
  # Hyp faite : duv = infini
  #Y = []
  arks = [(i,j) for i in 1:graph.n for j in 1:graph.n if graph.adj[i][j]]
  longest_path_length = -Inf
  x0 = Dict((i,j) => false for (i,j) in arks)
  m_0, y = shortest_path(graph, x0, arks)
  status = solve(m_0, relaxation=false)
  longest_x = x0
  if status != :Optimal
    println("Error : first path not found :")
    return
  end
  Y_path = get_path(y, graph.n, arks)
  #push!(Y, path)

  m = Model(solver=CbcSolver())
  @variable(m, x[(u,v) in arks], Bin)
  @constraint(m, sum(x[(u, v)] for (u,v) in Y_path) >= 1)
  @constraint(m, sum(x[(u, v)] for (u,v) in arks) <= graph.k)
  @constraint(m, xoradj[(u,v) in arks], x[(u, v)] <= graph.adj[u][v])

  println("Initialization done.")
  while status == :Optimal
    println("Solving master...")
    # m = Model(solver=CbcSolver())
    # @variable(m, x[(u,v) in arks], Bin)
    # @constraint(m, cons[path in Y], sum(x[(u, v)] for (u,v) in Y_path) >= 1)
    # @constraint(m, sum(x[(u, v)] for (u,v) in arks) <= graph.k)
    # @constraint(m, xoradj[(u,v) in arks], x[(u, v)] <= graph.adj[u][v])

    status = solve(m, relaxation=false)
    if status != :Optimal
      println("Master problem unfeasible")
      break
    end
    x_val = getvalue(x) # [getvalue(x[(u, v)]) == 1.0 for (u,v) in arks]
    println("Master solved with status ", status, " and x = ", get_x(x, graph.n, arks))


    println("Solving slave...")
    @constraint(m_0, [(u,v) in arks],  getvalue(x[(u,v)]) + y[(u,v)] <= 1)

    model, slave_path =  shortest_path(graph, x_val, arks)
    # model = m_0
    # slave_path = y
    slave_status = solve(model, relaxation=false)
    
    if slave_status != :Optimal
      println("Error : sub problem is unfeasible for x = ", x_val)
      break
    end
    if getobjectivevalue(model) > longest_path_length
      longest_x = getvalue(x)
      longest_path_length = getobjectivevalue(model)
    end
    # push!(Y, get_path(slave_path, graph.n, arks))
    Y_path = get_path(slave_path, graph.n, arks)
    @constraint(m, sum(x[(u, v)] for (u,v) in Y_path) >= 1)
  end
  print("Final solution is ", get_x_from_value(longest_x, graph.n, arks), " with value ", longest_path_length)
end
