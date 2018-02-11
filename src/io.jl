function get_erf_integral()
    data = readdlm("assets/erf_integral.dat")
    erf_table = data[:,2]
end

function write_sim_parameters()
    writedlm("output/parameters.csv", [@sprintf("%.5e", dtm) @sprintf("%d",nz) @sprintf("%.5e", dt_print) @sprintf("%.5e", tm_quiet) @sprintf("%.5e", L)], ",")
end

function open_files()
    	c1::Int = 0
    	c2::Int = 0
    	c3::Int = 0

    	for i = 1:(nspec+1)*3+2
    		if i==1
    			filename = "r.csv"
    		elseif i==(nspec+1)*3+2
    			filename = "efield.csv"
    		elseif (i>1) & (i<=(nspec+1)+1)
    			c1 += 1
    			filename = "vel" * string(c1) * ".csv"
    		elseif (i>(nspec+1)+1) & (i<=2*(nspec+1)+1)
                c2 += 1
    			filename = "den" * string(c2) * ".csv"
    		elseif (i>2*(nspec+1)+1) & (i<=3*(nspec+1)+1)
    			c3 += 1
    			filename = "temp" * string(c3) * ".csv"
    		end
    		open("output/" * filename, "w+")
    	end
end


function write_data()
	c1::Int = 0
	c2::Int = 0
	c3::Int = 0
    var = Array{Float64}(nz)

	for i = 1:(nspec+1)*3+2
		if i==1
			var = r
            println(r)
            filename = "r.csv"
		elseif i==(nspec+1)*3+2
			var = Efield
            filename = "efield.csv"
		elseif (i>1) & (i<=(nspec+1)+1)
			c1 += 1
			var = u[:,c1]
            filename = "vel" * string(c1) * ".csv"
		elseif (i>(nspec+1)+1) & (i<=2*(nspec+1)+1)
			c2 += 1
			if c2<=nspec
				var = 1.e-6 * rho[:,c2]	/ mi[c2]
			else #electrons
				var = 1.e-6 * rho[:,c2]	/ me
			end
            filename = "den" * string(c2) * ".csv"
		elseif (i>2*(nspec+1)+1) & (i<=3*(nspec+1)+1)
			c3 += 1
			var = T[:,c3] / qe_C
            filename = "temp" * string(c3) * ".csv"
		end
		for k = 1:nz
			write("output/" * filename, var[k])
		end
	end

end
