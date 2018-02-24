push!(LOAD_PATH, "./init")
push!(LOAD_PATH, "./src")

using Constants
using Get_Sim_Params

include("../init/process_params.jl")
const dr, mi_g, mi, qi, nquiet, cs, CFL, dtm,
        lm, nsmooth, vsmooth, tsmooth = process_params()

include("types.jl")
hydro = hydro_vars(nz, nspec, neqi)

include("../init/init_sim.jl")
const r = init_spacegrid(dr)
const eps_visc = init_art_viscosity(geom)
if !restart
    hydro = init_variables(hydro)
end
hydro = do_smoothing(hydro)
hydro = init_fluxes(hydro)

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
    println("=================================")
    println("ne = ", ne)
    calculate_collisions!(erf_table)
    source_terms!(nq)

    #----- corrector -----------
    #---------------------------
    corrector!(nq,j)

    #----- final step ----------
    U1D[2:nq-1,1:neqi+1] = 0.5 * (U1D[2:nq-1,1:neqi+1] + U1D_c[2:nq-1,1:neqi+1]) +
            eps_visc[1:nq-2,1:neqi+1] .* (U1D[3:nq,1:neqi+1] - 2 * U1D[2:nq-1,1:neqi+1] + U1D[1:nq-2,1:neqi+1])

    update_variables!(U1D)
    calculate_collisions!(erf_table)
    source_terms!(nq)

    if is_print
        write_data()
        is_print = false
    end
end

println("---Simulation ended---")

close_files()
