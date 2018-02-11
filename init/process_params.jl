function process_params()

    mi = zeros(Float64, nspec)
    mi_g = zeros(Float64, nspec)
    qi = zeros(Float64, nspec)

    dr = L / nz
    if geom=="spherical"
        dr = -dr #i.e., dr is negative
    end

    for ispec = 1:nspec
        mi_g[ispec] = Ai[ispec] * mp_g #mass of ion in grams
        mi[ispec] = 1.e-3 * mi_g[ispec]
        qi[ispec] = Zi[ispec] * qe #charge of species in statcoulomb
    end

    nquiet = floor(1.18e-2 / abs(dr))
    cs = sqrt(2 * g * maximum(temp0) / mean(mi)) #it's an "average" speed of sound
    CFL = abs(dr) / cs
    dtm = coeff_CFL * CFL #may need to be varied to avoid NaN when tm > tm_coll
    lm = dtm / dr

    nsmooth::Int32 = floor(abs(dr_smooth/dr)) #only used if smoothing==.true.
    vsmooth = nsmooth  #smoothing for velocity
    tsmooth = nsmooth  #smoothing for temperature


    return dr, mi_g, mi, qi, nquiet, cs, CFL, dtm, lm,
            nsmooth, vsmooth, tsmooth
end
