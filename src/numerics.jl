using Constants

function predictor!(hydro, nq, jj)

    nq::Int = nq
	lm_ = dtm_ / dr

	#take care of ions
	for i = 1:neqi
		hydro.U1D_p[2:nq-1,i] = hydro.U1D[2:nq-1,i] - lm_ * (hydro.F1D[3:nq,i] - hydro.F1D[2:nq-1,i])
				+ dtm_ * (phi * hydro.G1D[3:nq,i] + (1-phi) * hydro.G1D[2:nq-1,i])
		if (geom=="spherical") #add correction for spherical geometry
			hydro.U1D_p[2:nq-1,i] = hydro.U1D_p[2:nq-1,i] - 2. * dtm_ / r[2:nq-1]
						*(phi * ( hydro.F1D[3:nq,i] - hydro.C1D[3:nq,i] ) &
							+ (1-phi) * ( hydro.F1D[2:nq-1,i] - hydro.C1D[2:nq-1,i]))
		end
	end

	#take care of electrons
	if mod(jj,2) == 1
		hydro.U1D_p[2:nq-1,neqi+1] = hydro.U1D[2:nq-1,neqi+1] -
                lm_ * (hydro.F1D[3:nq,neqi+1] - hydro.F1D[2:nq-1,neqi+1]) +
				dtm_ * (phi * hydro.G1D[3:nq,neqi+1] + (1-phi) * hydro.G1D[2:nq-1,neqi+1])
		if geom=="spherical" #add correction for spherical geometry
			hydro.U1D_p[2:nq-1,neqi+1] = hydro.U1D_p[2:nq-1,neqi+1] - 2. * dtm_ / r[2:nq-1] *
						(phi * hydro.F1D[3:nq,neqi+1] + (1-phi) * hydro.F1D[2:nq-1,neqi+1])
		end
	else
		hydro.U1D_p[2:nq-1,neqi+1] = hydro.U1D[2:nq-1,neqi+1] -
                lm_ * (hydro.F1D[2:nq-1,neqi+1] - hydro.F1D[1:nq-2,neqi+1]) +
				dtm_ * ((1-phi) * hydro.G1D[2:nq-1,neqi+1] + phi * hydro.G1D[1:nq-2,neqi+1])
		if (geom=="spherical") #add correction for spherical geometry
			hydro.U1D_p[2:nq-1,neqi+1] = hydro.U1D_p[2:nq-1,neqi+1] - 2. * dtm_ / r[2:nq-1] *
						((1-phi) * hydro.F1D[2:nq-1,neqi+1] + phi * hydro.F1D[1:nq-2,neqi+1])
		end
	end
    return hydro
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function corrector!(hydro, nq, jj)

    nq::Int = nq
	lm_ = dtm_ / dr

	for i = 1:neqi
		hydro.U1D_c[2:nq-1,i] = hydro.U1D_p[2:nq-1,i] -
                lm_ * (hydro.F1D[2:nq-1,i] - hydro.F1D[1:nq-2,i]) +
				dtm_ * ((1-phi) * hydro.G1D[2:nq-1,i] + phi * hydro.G1D[1:nq-2,i])
		if geom=="spherical"
			hydro.U1D_c[2:nq-1,i] = hydro.U1D_c[2:nq-1,i] - 2. * dtm_ ./ r[2:nq-1]  .*
						((1-phi) * (hydro.F1D[2:nq-1,i] - hydro.C1D[2:nq-1,i]) +
						phi * (hydro.F1D[1:nq-2,i] - hydro.C1D[1:nq-2,i]))
		end
	end

	#take care of electrons
	if mod(jj,2) == 1
		hydro.U1D_c[2:nq-1,neqi+1] = hydro.U1D_p[2:nq-1,neqi+1] -
                lm_ * (hydro.F1D[2:nq-1,neqi+1] - hydro.F1D[1:nq-2,neqi+1])  +
				dtm_ * ((1-phi) * hydro.G1D[2:nq-1,neqi+1] + phi * hydro.G1D[1:nq-2,neqi+1])
		if geom=="spherical" #add correction for spherical geometry
			hydro.U1D_c[2:nq-1,neqi+1] = hydro.U1D_c[2:nq-1,neqi+1] - 2. * dtm_ ./ r[2:nq-1]  .*
						((1-phi) * hydro.F1D[2:nq-1,neqi+1] + phi * hydro.F1D[1:nq-2,neqi+1])
		end
	else
		hydro.U1D_c[2:nq-1,neqi+1] = hydro.U1D_p[2:nq-1,neqi+1] -
                lm_ *  (hydro.F1D[3:nq,neqi+1] - hydro.F1D[2:nq-1,neqi+1]) +
				dtm_ * (phi * hydro.G1D[3:nq,neqi+1] + (1-phi) * hydro.G1D[2:nq-1,neqi+1])
		if geom=="spherical" #add correction for spherical geometry
			hydro.U1D_c[2:nq-1,neqi+1] = hydro.U1D_c[2:nq-1,neqi+1] - 2. * dtm_ ./ r[2:nq-1] .*
						(phi * hydro.F1D[3:nq,neqi+1] + (1-phi) * hydro.F1D[2:nq-1,neqi+1])
		end
	end

    return hydro
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function update_variables!(UU, hydro)
	# need to initialize electron density and speed to zero
	hydro.rho[:,nspec+1] = 0.
	hydro.ne = zeros(Float64, nz)
	hydro.u[:,nspec+1] = 0.

	for i = 1:nspec
		hydro.rho[:,i] = UU[:,3*(i-1)+1]
		hydro.u[:,i] = UU[:,3*(i-1)+2] ./ hydro.rho[:,i]
		hydro.T[:,i] = (g-1) * mi[i] * (UU[:,3*(i-1)+3] ./ hydro.rho[:,i] - 0.5*hydro.u[:,i].^2)
		hydro.p[:,i] = hydro.rho[:,i] .* hydro.T[:,i] / mi[i]

		hydro.F1D[:,3*(i-1)+1] = hydro.rho[:,i] .* hydro.u[:,i]
		hydro.F1D[:,3*(i-1)+2] = hydro.rho[:,i] .* hydro.u[:,i].^2 + hydro.p[:,i]
		hydro.F1D[:,3*(i-1)+3] = hydro.rho[:,i] .* hydro.u[:,i] .*
            UU[:,3*(i-1)+3] ./ hydro.rho[:,i] + hydro.p[:,i] .* hydro.u[:,i]

		if geom=="spherical"  # define correcting matrix
			hydro.C1D[:,3*(i-1)+1] = 0.
			hydro.C1D[:,3*(i-1)+2] = copy(hydro.p[:,i])
			hydro.C1D[:,3*(i-1)+3] = 0.
		end

	    # electrons: quasi-neutrality + zero-current equation
	    hydro.rho[:,nspec+1] = hydro.rho[:,nspec+1] + me * Zi[i] * hydro.rho[:,i] / mi[i]  # quasi-neutrality
	    hydro.u[:,nspec+1] = hydro.u[:,nspec+1] +
                Zi[i] * hydro.rho[:,i] .* hydro.u[:,i] / mi[i]  # zero-current condition

		hydro.ni[:,i] = hydro.rho[:,i] / mi[i]
		hydro.ne = hydro.ne + Zi[i] * hydro.ni[:,i]
	end
	hydro.u[:,nspec+1] = hydro.u[:,nspec+1] ./ hydro.ne # correct with ne

	# electron temperature and pressure
	hydro.T[:,nspec+1] = (g-1) * me * (UU[:,neqi+1] ./ hydro.rho[:,nspec+1] - 0.5*hydro.u[:,nspec+1].^2 )
	hydro.p[:,nspec+1] = hydro.rho[:,nspec+1] .* hydro.T[:,nspec+1] / me
	hydro.F1D[:,neqi+1] = hydro.u[:,nspec+1] .* UU[:,neqi+1]  + hydro.p[:,nspec+1] .* hydro.u[:,nspec+1]

	# limit max temperature by Tmax and avoid negative temperatures
	# for j = 1:nspec+1
	# 	T[:,j] = min( T[:,j], Tmax  )
	#     T[:,j] = max( T[:,j], 0.  ) # avoid negative temperature
	# end

    return hydro
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function source_terms!(hydro, k_DT, ke, Qextra, nq)

    Q_DT = zeros(Float64, nz, nspec + 1)
    hydro.G1D[:, :] = 0.0
    q_diff = zeros(Float64, nz)
    q_FS = zeros(Float64, nz)

	if heattransfer_switch & electron_switch
		for i = 1:nspec+1
			for j = 1:nspec+1
				if j != i
					Q_DT[:,i] = Q_DT[:,i] + k_DT[:,i,j] .* (hydro.T[:,j] - hydro.T[:,i])
				end
			end
		end
		Q_DT[:,nspec+1] = Q_DT[:,nspec+1] + Qextra  #electrons
	elseif heattransfer_switch & ~electron_switch   #only ions
		for i = 1:nspec
			for j = 1:nspec
				if j != i
					Q_DT[:,i] = Q_DT[:,i] + k_DT[:,i,j] .* (hydro.T[:,j] - hydro.T[:,i])
				end
			end
		end
	end

	if efield_switch & electron_switch
		hydro.Efield[2:nq-1] = - 0.5 * (hydro.p[3:nq,nspec+1] - hydro.p[1:nq-2,nspec+1]) /
		                  dr ./ (qe_C * hydro.ne[2:nq-1] )
	else
		hydro.Efield[:] = 0.
	end

	for i = 1:nspec
		hydro.G1D[:,3*(i-1)+2] = Zi[i] * qe_C * hydro.ni[:,i] .* hydro.Efield  #source term
		hydro.G1D[:,3*(i-1)+2] = hydro.G1D[:,3*(i-1)+2] + hydro.R1D[:,3*(i-1)+2]
		hydro.G1D[:,3*(i-1)+3] = Q_DT[:,i] + Zi[i] * qe_C * hydro.ni[:,i] .* hydro.Efield .* hydro.u[:,i]
		hydro.G1D[:,3*(i-1)+3] = hydro.G1D[:,3*(i-1)+3] + hydro.R1D[:,3*(i-1)+3]
	end

	q_diff[2:nq-1] = 1. / dr * ke[2:nq-1] .* (hydro.T[3:nq,nspec+1] - hydro.T[2:nq-1,nspec+1])

	vth_e = sqrt.(hydro.T[:,nspec+1] / me)
	q_FS[2:nq-1] = hydro.rho[3:nq,nspec+1] / me .* hydro.T[3:nq,nspec+1] .* vth_e[3:nq]  #ensures that it is forward-differencing later
	q_diff[:] = sign.(q_diff) .* min.(flimit * q_FS, abs.(q_diff))

	hydro.G1D[:,neqi+1] = Q_DT[:,nspec+1] - qe_C * hydro.ne .* hydro.Efield .* hydro.u[:,nspec+1]

	if electron_heat_conduction
		hydro.G1D[2:nq-1,neqi+1] = hydro.G1D[2:nq-1,neqi+1] + 1. / dr * (q_diff[2:nq-1] - q_diff[1:nq-2])
		if geom=="spherical"  #add correction due to spherical geometry
			hydro.G1D[:,neqi+1] = hydro.G1D[:,neqi+1] + 2. * q_diff ./ r
		end
	end

	#Lindl, ch. 3 eq. (26) - rho in g/cm3, T in keV
	if bremsstrahlung
		hydro.G1D[:,neqi+1] = hydro.G1D[:,neqi+1] - 3.0e16 * (1.e3 * me * hydro.ne).^2 .*
			             (1.e-3 * hydro.T[:,nspec+1] / qe_C)
	end

    return hydro
end
