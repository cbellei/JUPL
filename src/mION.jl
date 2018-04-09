push!(LOAD_PATH, "./init")
push!(LOAD_PATH, "./src")
using Glob

include("io.jl")
output_files = glob("output/*.csv")
for file in output_files
    rm(file)
end
jpg_files = glob("output/*.jpg")
for file in jpg_files
    rm(file)
end

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
init_predictor_corrector(hydro)

erf_table = get_erf_integral()
write_sim_parameters()
open_files()
write_data(hydro)
close_files()

include("utils.jl")
include("boundary.jl")
include("numerics.jl")
include("collisions.jl")

nq = 0
j = 0
while tm <= maxTime
    j += 1
    tm += dtm
    is_print = is_print_time(tm, dt_print)

    nq, time1, time2 = is_quiet_time!(time1, time2, tm, tm_quiet, geom, nz, nq, dt_print, dtm)

    apply_BC!(hydro)

    #----- predictor -----------
    #---------------------------
    predictor!(hydro, nq, j)
    update_variables!(hydro.U1D_p, hydro)
    hydro, k_DT, ke, Qextra = calculate_collisions!(hydro, erf_table)
    source_terms!(hydro, k_DT, ke, Qextra, nq)
      # #----start debug
      # if j==1
      #     open_files()
      #     write_data(hydro)
      #     close_files()
      #     write_all_data(hydro)
      #     hydro, k_DT, ke, Qextra = calculate_collisions!(hydro, erf_table)
      #     systemerror(0)
      # end
      # #----end debug

    #----- corrector -----------
    #---------------------------
    corrector!(hydro, nq, j)

    #----- final step ----------
    hydro.U1D[2:nq-1,1:neqi+1] = 0.5 * (hydro.U1D[2:nq-1,1:neqi+1] + hydro.U1D_c[2:nq-1,1:neqi+1]) +
        eps_visc[1:nq-2,1:neqi+1] .*
        (hydro.U1D[3:nq,1:neqi+1] - 2 * hydro.U1D[2:nq-1,1:neqi+1] + hydro.U1D[1:nq-2,1:neqi+1])

    update_variables!(hydro.U1D, hydro)
    hydro, k_DT, ke, Qextra = calculate_collisions!(hydro, erf_table)
    source_terms!(hydro, k_DT, ke, Qextra, nq)

    if is_print
        open_files()
        write_data(hydro)
        close_files()
        is_print = false
    end

end

println("---Simulation ended---")

close_files()
