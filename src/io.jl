function get_erf_integral()
    data = readdlm("assets/erf_integral.dat")
    erf_table = data[:,2]
end

function write_sim_parameters()
    writedlm("output/parameters.csv", [@sprintf("%.5e", dtm) @sprintf("%d",nz) @sprintf("%.5e", dt_print) @sprintf("%.5e", tm_quiet) @sprintf("%.5e", L)], ",")
end
