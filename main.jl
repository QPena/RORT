include("heuri.jl")
include("pl_1.jl")

function main()
	inst = generate("C:\\Users\\penaq\\Desktop\\MPRO\\RORT\\10_10_5_50_50_1.dat")
	#println(problem_1(inst))
	println(heuri!(inst))

end

