push!(LOAD_PATH, "./init")
push!(LOAD_PATH, "./src")

using Constants
using Get_Sim_Params
using Allocate_Arrays

include("../init/process_params.jl")
const dr, mi_g, mi, qi, nquiet, cs, CFL, dtm,
        lm, nsmooth, vsmooth, tsmooth = process_params()

include("../init/init_sim.jl")
const r = init_spacegrid(dr)
const eps_visc = init_art_viscosity(geom)
if !restart
    init_variables()
end
do_smoothing()
init_fluxes()

include("io.jl")
get_erf_integral()
write_sim_parameters()
open_files()
write_data()

println("C'est parti!")

close_files()
