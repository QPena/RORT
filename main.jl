include("heuri.jl")
include("pl_1.jl")
include("pl_3.jl")

function main()
	inst = generate("C:\\Users\\penaq\\Desktop\\MPRO\\RORT\\5_4_3_10_10_1.dat")
	#@time pb1 = problem_1(inst)
	#@time pb3 = master(inst);
	@time heuri = heuri!(inst)

end

