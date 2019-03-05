include("heuri.jl")
include("pl_1.jl")
include("pl_2.jl")
include("pl_3.jl")

workingdir = "C:\\Users\\penaq\\Desktop\\MPRO\\RORT\\Instances"

function main(pl1=false, pl2=false, pl3=false, heuri=false, filetemplate = "", modes = [1,2,3,4,5,6,7])
	for file in readdir(workingdir)
		if startswith(file, filetemplate)
			println("\n",file)
			inst = generate(workingdir * "\\" * file)
			if pl1
				print("PL 1 : ")
				@time pb1 = problem_1(inst)
				println(pb1)
			end

			if pl2
				print("PL 2 : ")
				@time pb2 = problem_2(inst)
				println(pb2)
			end

			if pl3
				print("PL 3 : ")
				@time pb3 = master(inst);
				println(pb3)
			end

			if heuri
				println("Heuristique : ")
				@time inf, sup, minf, msup = heuristicSolve!(workingdir * "\\" * file)
				println("Globale borne inf : ", inf, " (mode ", minf, ")")
				println("Globale borne sup : ", sup, " (mode ", msup, ")")
				println("Best gap          : ", 100*(sup-inf)/sup, "%")
				#for mode in modes
				#	inst = generate(workingdir * "\\" * file)
				#	@time heursol = heuri!(inst, mode)
				#	println("Mode ", mode, " : ", heursol)
				#end
			end
		end
	end

end


function generate_instances(l,c,k,maxc,maxd, N)
	shift = 0
	print_maxd = maxd>0 ? string(maxd) : "inf"
	for i = 1:N
		filename = string(l) * "_" * string(c) * "_" * string(k) * "_" * string(maxc) * "_" *print_maxd*"_"*string(i+shift)*".dat"
		while filename in readdir(workingdir)
			shift += 1
			filename = string(l) * "_" * string(c) * "_" * string(k) * "_" * string(maxc) * "_" *print_maxd*"_"*string(i+shift)*".dat"
		end
		generate(l,c,k,maxc,maxd,"C:\\Users\\penaq\\Desktop\\MPRO\\RORT\\Instances\\" * filename )
		println("Generated ", filename)
	end
end
# pb3 : 157 en 6s