function get_erf_integral()
    data = readdlm("assets/erf_integral.dat")
    erf_table = data[:,2]
end

function write_sim_parameters()
    filename = "output/parameters.csv"
    writedlm(filename, [@sprintf("%.5e", dtm) @sprintf("%d",nz) @sprintf("%.5e", dt_print) @sprintf("%.5e", tm_quiet) @sprintf("%.5e", L)], ",")
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
            f = open("output/" * filename, "a")
    		@eval $(Symbol("f_$i")) = $f
    	end
end


function write_data(hydro)
	c1::Int = 0
	c2::Int = 0
	c3::Int = 0
    var = Array{Float64}(nz)

	for i = 1:(nspec+1)*3+2
		if i==1
			var = r
            filename = "r.csv"
		elseif i==(nspec+1)*3+2
			var = hydro.Efield
            filename = "efield.csv"
		elseif (i>1) & (i<=(nspec+1)+1)
			c1 += 1
			var = hydro.u[:,c1]
            filename = "vel" * string(c1) * ".csv"
		elseif (i>(nspec+1)+1) & (i<=2*(nspec+1)+1)
			c2 += 1
			if c2<=nspec
				var = 1.e-6 * hydro.rho[:,c2]	/ mi[c2]
			else #electrons
				var = 1.e-6 * hydro.rho[:,c2]	/ me
			end
            filename = "den" * string(c2) * ".csv"
		elseif (i>2*(nspec+1)+1) & (i<=3*(nspec+1)+1)
			c3 += 1
			var = hydro.T[:,c3] / qe_C
            filename = "temp" * string(c3) * ".csv"
		end
        f = @eval $(Symbol("f_$i"))
		for k = 1:nz
			writedlm(f, [@sprintf("%.5e", var[k])])
		end
	end
end

function close_files()
    for i = 1:(nspec+1)*3+2
        f = @eval $(Symbol("f_$i"))
        close(f)
    end
end




#---------------------------------------------------------
function write_all_data(hydro)
	c1::Int = 0
	c2::Int = 0
	c3::Int = 0
    var = Array{Float64}(nz)

    filename = "./output/U1D.csv"
    var = hydro.U1D
    writedlm(filename, var, ",")

    filename = "./output/U1D_p.csv"
    var = hydro.U1D_p
    writedlm(filename, var, ",")

    filename = "./output/U1D_c.csv"
    var = hydro.U1D_c
    writedlm(filename, var, ",")

    filename = "./output/F1D.csv"
    var = hydro.F1D
    writedlm(filename, var, ",")

    filename = "./output/G1D.csv"
    var = hydro.G1D
    writedlm(filename, var, ",")

    filename = "./output/R1D.csv"
    var = hydro.R1D
    writedlm(filename, var, ",")

    filename = "./output/C1D.csv"
    var = hydro.C1D
    writedlm(filename, var, ",")

    filename = "./output/T.csv"
    var = hydro.T
    writedlm(filename, var, ",")

end
