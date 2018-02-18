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
erf_table = get_erf_integral()
write_sim_parameters()
open_files()
write_data()

include("utils.jl")
include("boundary.jl")
include("numerics.jl")
include("collisions.jl")

nq = 0
j = 0
U1D_p, U1D_c = init_predictor_corrector()

while tm <= maxTime
    j += 1
    println(j)
    is_print = is_print_time(tm, dt_print)
    tm += dtm

    nq, time1, time2 = is_quiet_time!(time1, time2, tm, tm_quiet, geom, nz, nq)

    apply_BC!()
    #----- predictor -----------
    #---------------------------
    predictor!(nq, j)
    update_variables!(U1D_p)
    calculate_collisions!(erf_table)
end

println("---Simulation ended---")

close_files()
