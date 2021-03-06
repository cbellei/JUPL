mutable struct hydro_vars
    U1D::Array{Float64,2}
    U1D_p::Array{Float64,2}
    U1D_c::Array{Float64,2}
    F1D::Array{Float64,2}
    G1D::Array{Float64,2}
    R1D::Array{Float64,2}
    C1D::Array{Float64,2}
    rho::Array{Float64,2}
    u::Array{Float64,2}
    T::Array{Float64,2}
    p::Array{Float64,2}
    Efield::Array{Float64,1}
    ni::Array{Float64,2}
    ne::Array{Float64,1}

    function hydro_vars(nz::Int64, nspec::Int64, neqi::Int64)
        U1D = zeros(Float64, nz, neqi+1)
        U1D_p = zeros(Float64, nz, neqi+1)
        U1D_c = zeros(Float64, nz, neqi+1)
        F1D = zeros(Float64, nz, neqi+1)
        G1D = zeros(Float64, nz, neqi+1)
        R1D = zeros(Float64, nz, neqi + 1)
        C1D = zeros(Float64, nz, neqi)
        p = zeros(Float64, nz, nspec+1)
        u = zeros(Float64, nz, nspec+1)
        rho = zeros(Float64, nz, nspec+1)
        T = zeros(Float64, nz, nspec+1)
        Efield = zeros(Float64, nz)
        ni = zeros(Float64, nz, nspec)
        ne = zeros(Float64, nz)

        new(
            U1D,
            U1D_p,
            U1D_c,
            F1D,
            G1D,
            R1D,
            C1D,
            p,
            u,
            rho,
            T,
            Efield,
            ni,
            ne
        )
    end
end
