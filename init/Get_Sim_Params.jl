module Get_Sim_Params
    using Constants
    import YAML

    path = ARGS[2] #"./runs/run01/run01.yml"
    data = YAML.load(open(path))

    export geom, nz, nregions, L, smoothing
    export coeff_CFL, maxTime, smoothing, dr_smooth, tm_quiet
    export dt_print, efield_switch, electron_switch
    export heattransfer_switch, flimit, electron_heat_conduction
    export friction_switch, Lab_min, check_frequency, smax, ion_viscosity
    export bremsstrahlung, restart, nedges
    export nspec, Ai, Zi, r0, den0, temp0, vel0
    export rmin, neqi, eps_visc_max, eps_compress
    export smooth_vars, smooth_type
    export boundary, phi

    key = "grid"
    const geom = data[key]["geom"]
    const nz = data[key]["nz"]
    const nregions = data[key]["nregions"]
    const L = data[key]["L"]
    const boundary = data[key]["boundary"]

    key = "output"
    const dt_print = data[key]["dt_print"]

    key = "numerics"
    const restart = data[key]["restart"]
    const tm_quiet = data[key]["tm_quiet"]
    const coeff_CFL = data[key]["coeff_CFL"]
    const maxTime = data[key]["maxTime"]
    const smoothing = data[key]["smoothing"]
    const dr_smooth = data[key]["dr_smooth"]
    const eps_visc_max = data[key]["eps_visc_max"]
    const eps_compress = data[key]["eps_compress"]
    const nedges = data[key]["nedges"]
    const smooth_vars = data[key]["smooth_vars"]
    const phi = data[key]["phi"]

    key = "physics"
    subkey = "electric_field"
    const efield_switch = data[key][subkey]["efield_switch"]
    const electron_switch = data[key][subkey]["electron_switch"]
    subkey = "heat_conduction"
    const heattransfer_switch = data[key][subkey]["heattransfer_switch"]
    const flimit = data[key][subkey]["flimit"]
    const electron_heat_conduction = data[key][subkey]["electron_heat_conduction"]
    subkey = "collisions"
    const friction_switch = data[key][subkey]["friction_switch"]
    const Lab_min = data[key][subkey]["Lab_min"]
    const check_frequency = data[key][subkey]["check_frequency"]
    const smax = data[key][subkey]["smax"]
    const ion_viscosity = data[key][subkey]["ion_viscosity"]
    subkey = "radiation"
    const bremsstrahlung = data[key][subkey]["bremsstrahlung"]

    key = "species"
    const nspec = data[key]["nspec"]
    ai = zeros(Float64, nspec)
    zi = zeros(Float64, nspec)
    r = zeros(Float64, nspec, nregions + 1)
    den = zeros(Float64, nspec, nregions)
    temp = zeros(Float64, nspec, nregions)
    vel = zeros(Float64, nspec, nregions)
    smooth_type = Array{String}(nedges, nspec)
    for ispec = 1:nspec
        ai[ispec] = data[key]["species"*string(ispec)]["Ai"]
        zi[ispec] = data[key]["species"*string(ispec)]["Zi"]
        r_species = data[key]["species"*string(ispec)]["r"]
        for iedge = 1:nregions + 1
            if isa(r_species[iedge], Float64)
                r[ispec, iedge] = r_species[iedge]
            else
                r[ispec, iedge] = L
            end
        end
        den_species = data[key]["species"*string(ispec)]["den"]
        temp_species = qe_C * data[key]["species"*string(ispec)]["temp"]
        vel_species = data[key]["species"*string(ispec)]["vel"]
        for ireg = 1:nregions
            den[ispec, ireg] = den_species[ireg]
            temp[ispec, ireg] = temp_species[ireg]
            vel[ispec, ireg] = vel_species[ireg]
        end
        smooth_type_species = data[key]["species"*string(ispec)]["smooth_type"]
        for iedge = 1:nedges
            smooth_type[iedge, ispec] = smooth_type_species
        end

    end
    const Ai = ai
    const Zi = zi
    const r0 = r
    const den0 = den
    const temp0 = temp
    const vel0 = vel

    if geom=="slab"
        const rmin = 0.
    else #spherical
        const rmin = 2.5e-6
    end

    const neqi = 3*nspec
end
