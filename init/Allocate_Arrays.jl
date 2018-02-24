
module Allocate_Arrays
    using Get_Sim_Params

    export xi, q_FS, q_diff, xi12, nu0, ke, du, vth_e
    export Te, Efield, ne, Qextra, p, u, rho, T
    export Q_DT, ni, mu_ion, xiab, L_ab, k_DT, nu_DT, U1D, F1D
    export R1D, U1D_c, U1D_p, G1D, C1D
    export dr_shell, r_shell
    export T0, N0, V0

    xi = zeros(Float64, nz)
    q_FS = zeros(Float64, nz)
    q_diff = zeros(Float64, nz)
    xi12 = zeros(Float64, nz)
    nu0 = zeros(Float64, nz)
    ke = zeros(Float64, nz)
    du = zeros(Float64, nz)
    vth_e = zeros(Float64, nz)
    Te = zeros(Float64, nz)
    Efield = zeros(Float64, nz)
    ne = zeros(Float64, nz)
    Qextra = zeros(Float64, nz)

    p = zeros(Float64, nz, nspec+1)
    u = zeros(Float64, nz, nspec+1)
    rho = zeros(Float64, nz, nspec+1)
    T = zeros(Float64, nz, nspec+1)
    Q_DT = zeros(Float64, nz, nspec+1)

    ni = zeros(Float64, nz, nspec)
    mu_ion = zeros(Float64, nz, nspec)

    xiab = zeros(Float64, nz, nspec, nspec)
    L_ab = zeros(Float64, nz, nspec, nspec)
    k_DT = zeros(Float64, nz, nspec+1, nspec+1)
    nu_DT = zeros(Float64, nz, nspec+1, nspec+1)

    U1D = zeros(Float64, nz, neqi+1)
    F1D = zeros(Float64, nz, neqi+1)
    R1D = zeros(Float64, nz, neqi+1)
    U1D_p = zeros(Float64, nz, neqi+1)
    U1D_c = zeros(Float64, nz, neqi+1)
    G1D = zeros(Float64, nz, neqi+1)

    C1D = zeros(Float64, nz, neqi)

    r_shell = zeros(Float64, nspec)
    dr_shell = zeros(Float64, nspec)

    r0 = zeros(Float64, nspec, nregions+1)
    den0 = zeros(Float64, nspec, nregions+1)
    temp0 = zeros(Float64, nspec, nregions+1)
    vel0 = zeros(Float64, nspec, nregions+1)

    T0 = zeros(Float64, nz, nspec)
    N0 = zeros(Float64, nz, nspec)
    V0 = zeros(Float64, nz, nspec)

end
