@doc """
Type représentant une entrée du problème d'interdiction de plus court chemins.

Il contient 6 attributs:
- `n` est le nombre de noeuds du graphe
- `m` est le nombre d'arcs du graphe
- `k` est le nombre d'arcs pénalisés par l'attaquant
- `adj` est un tableau de tableau de `n` x `n` booléens, `adj[i][j]` est vrai si (i, j) est un arc du graphe
- `c` est un tableau de tableau de `n` x `n` entiers, `c[i][j]` est le coût de l'arc (i, j).
- `d` est un tableau de tableau de `n` x `n` entiers, `d[i][j]` est le coût de l'arc (i, j).

Il n'est pas nécessaire de construire vous même les entrées, vous avez 3 fonctions pour cela:
- `generate(l, c, k, maxc, maxd)` génère un graphe de type grille
- `generate(l, c, k, maxc, maxd, filename)` génère un graphe de type grille et écrit le tout dans le fichier filename. 
- `generate(filename) renvoie l'instance générée dans le fichier filename
""" ->
mutable struct Data
  n::Int # nbre sommets
  m::Int # nbre arcs
  k::Int # nbre arcs penalises
  adj::Array{Array{Bool}} # matrice adjacence
  c::Array{Array{Int}} # matrice couts
  d::Array{Array{Int}} # matrice couts additionels (penalite)
end

Data() = Data(0, 0, 0, [], [], []) # Constructeur par défaut


@doc """
generate(l::Int, c::Int, k::Int, maxc::Int, maxd::Int)
c >= 2
l >= 2

Fonction de génération d'instance.


Cette fonction génère et renvoie une instance du problème d'interdiction de plus courts chemins. Le graphe est une quasi-grille avec `l` * `c` + 2 noeuds : une source et un puits, et une grille de `l` lignes et `c` colonnes. La source est reliée à tous les noeuds de la premiève colonne. La dernière colonne est reliée au puits. Chaque noeud d'une colonne est reliée à au plus trois noeuds de la colonne suivante: le noeud v sur la même ligne, ainsi que les deux noeuds au dessus et en dessous de v, si ces noeuds existent. Enfin, sur toutes les colonnes, exceptées la première et la dernière, chaque noeud est relié au noeud de la ligne d'au dessus et de la ligne d'en dessous, si ce noeud existe.

Par exemple, dans l'exemple suivant où `l` = `c` = 3:\n
...1..4..7\n
s..2..5..8..t\n
...3..6..9

- Le noeud 1 est relié aux noeuds 4 et 5.
- Le noeud 2 est relié aux noeuds 4, 5, et 6.
- Le noeud 4 est relié aux noeuds 5, 7 et 8.
- Le noeud 5 est relié aux noeuds 4, 6, 7, 8 et 9.


L'attaquant peut pénaliser `k` arcs. Le coût initial d'un arc est choisi aléatoirement entre 1 et `maxc`. La pénalité d'un arc est choisie entre 1 et `maxd`. Si `maxd` est négatif alors la pénalité est infinie (elle est égale à `maxc` * le nombre d'arêtes + 1, supérieur au plus long des chemins de s à t dans le graphe initial).
""" ->
function generate(l::Int, c::Int, k::Int, maxc::Int, maxd::Int)
  sv = Data()
  sv.n = l * c + 2
  sv.m = (c-2)*(5*l-4)+5*l-2
  sv.k = k
    
  # Initialisation des matrices d'adjacence, de coût et de pénalité
  for i in 1:sv.n
    push!(sv.adj, [])
    push!(sv.c, [])
    push!(sv.d, [])
    for j in 1:sv.n
      push!(sv.adj[i], false) # on remplit la ligne de "Faux"
      push!(sv.c[i], 0) # on remplit la ligne de zéros
      push!(sv.d[i], 0) # on remplit la ligne de zéros
    end
  end
 
  arcs = []
  # Liaisons entre la source et la première colonne
  for i in 1:l
    push!(arcs, (1, i + 1))
  end  
  
  # Liaisons entre la dernière colonne et le puits
  for i in 1:l
    push!(arcs, (sv.n - i, sv.n))
  end  
  
  # Liaisons entre deux colonnes successives
  for column in 1:(c - 1)
    for line in 1:l

      u = (column - 1) * l + 1 + line
      push!(arcs, (u, u + l))
      if line > 1
        push!(arcs, (u, u + l - 1))
      end
      if line < l 
        push!(arcs, (u, u + l + 1))
      end

    end
  end  

  # Liaisons entre les noeuds d'une même colonne
  for column in 2:(c - 1)
    for line in 1:(l-1)
      u = (column - 1) * l + 1 + line
      push!(arcs, (u, u + 1))
      push!(arcs, (u + 1, u))
    end
  end

  # On rempli les matrics avec les arcs sélectionnés.
  for (i, j) in arcs
    sv.adj[i][j] = true
    sv.c[i][j] = rand(1:maxc)
    sv.d[i][j] = (maxd > 0) ? rand(1:maxd) : (maxc * sv.m + 1)
  end

  return sv
end


@doc """
generate(l::Int, c::Int, k::Int, maxc::Int, maxd::Int, filename::String)
c >= 2
l >= 2

Utilise la fonction `generate(l, c, k, maxc, maxd)` pour construire une instance puis écrit cette instance dans le fichier dont le chemin est `filename`. Enfin, elle renvoie l'instance générée.

Le format du fichier est le suivant:
- Une première ligne contenant les 3 entiers positifs n, m et k, correspondant respectivement au nombre de noeuds, d'arcs et d'arcs pénalisés par l'attaquant.
- m lignes contenant 4 entiers positifs u v c et d où u et v sont entre 1 et n (inclus) signifiant qu'il existe un arc entre u et v de coût c et de pénalité d.

""" ->
function generate(l::Int, c::Int, k::Int, maxc::Int, maxd::Int, filename::String)
  sv = generate(l, c, k, maxc, maxd)
  open(filename, "w") do file
    write(file, "$(sv.n) $(sv.m) $(sv.k)\n") 
    for i in 1:sv.n
      for j in 1:sv.n
        if sv.adj[i][j]
          write(file, "$i $j $(sv.c[i][j]) $(sv.d[i][j])\n")
        end
      end
    end
  
  end
  return sv
end


@doc """ 
generate(filename::String)

Lecture fichier de données.

Cette fonction lit un fichier dont le chemin est `filename` et renvoie un objet de type Data contenant les informations écrites dans le fichier.
Le fichier doit être au format généré par la fonction `generate(l::Int, c::Int, k::Int, maxc::Int, maxd::Int, filename::String)`
- Une première ligne contenant 3 entiers positifs n, m et k
- m lignes contenant 4 entiers positifs u v c et d où u et v sont entre 1 et n (inclus)
""" ->
function generate(filename::String)
  
  open(filename) do file
    lines = readlines(file)  # fichier de l'instance à resoudre
    
    sv = Data()

    # Nombre de noeuds et d'arêtes, et nombre d'arcs pénalisés
    line = lines[1]
    line_decompose=split(line)
    sv.n = parse(Int64, line_decompose[1])
    sv.m = parse(Int64, line_decompose[2])
    sv.k = parse(Int64, line_decompose[3])

    # Initialisation des matrices d'adjacence, de coût et de pénalité
    for i in 1:sv.n
      push!(sv.adj, [])
      push!(sv.c, [])
      push!(sv.d, [])
      for j in 1:sv.n
        push!(sv.adj[i], false) # on remplit la ligne de "Faux"
        push!(sv.c[i], 0) # on remplit la ligne de zéros
        push!(sv.d[i], 0) # on remplit la ligne de zéros
      end
    end

    # Lectures des informations concernant l'adjecence, les coûts et les pénalités
    for i in 2:(sv.m + 1)
      line = lines[i]
      line_decompose = split(line)
      som_deb = parse(Int64, line_decompose[1])
      som_fin = parse(Int64, line_decompose[2])
      sv.adj[som_deb][som_fin] = true
      sv.c[som_deb][som_fin] = parse(Int64,line_decompose[3])
      sv.d[som_deb][som_fin] = parse(Int64,line_decompose[4])
    end
    return sv
  end

end
