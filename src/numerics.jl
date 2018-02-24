using Constants

function predictor!(nq, jj)

    nq::Int = nq
	lm_ = dtm_ / dr

	#take care of ions
	for i = 1:neqi
		U1D_p[2:nq-1,i] = U1D[2:nq-1,i] - lm_ * ( F1D[3:nq,i] - F1D[2:nq-1,i] )
				+ dtm_ * ( phi * G1D[3:nq,i] + (1-phi) * G1D[2:nq-1,i] )
		if (geom=="spherical") #add correction for spherical geometry
			U1D_p[2:nq-1,i] = U1D_p[2:nq-1,i] - 2. * dtm_ / r[2:nq-1]
						*( phi * ( F1D[3:nq,i] - C1D[3:nq,i] ) &
							+ (1-phi) * ( F1D[2:nq-1,i] - C1D[2:nq-1,i] )  )
		end
	end

	#take care of electrons
	if mod(jj,2) == 1
		U1D_p[2:nq-1,neqi+1] = U1D[2:nq-1,neqi+1] - lm_ * ( F1D[3:nq,neqi+1] - F1D[2:nq-1,neqi+1] )
				+ dtm_ * ( phi * G1D[3:nq,neqi+1] + (1-phi) * G1D[2:nq-1,neqi+1] )
		if geom=="spherical" #add correction for spherical geometry
			U1D_p[2:nq-1,neqi+1] = U1D_p[2:nq-1,neqi+1] - 2. * dtm_ / r[2:nq-1]
						*( phi * F1D[3:nq,neqi+1] + (1-phi) * F1D[2:nq-1,neqi+1]   )
		end
	else
		U1D_p[2:nq-1,neqi+1] = U1D[2:nq-1,neqi+1] - lm_ * ( F1D[2:nq-1,neqi+1] - F1D[1:nq-2,neqi+1] )
				+ dtm_ * ( (1-phi) * G1D[2:nq-1,neqi+1] + phi * G1D[1:nq-2,neqi+1] )
		if (geom=="spherical") #add correction for spherical geometry
			U1D_p[2:nq-1,neqi+1] = U1D_p[2:nq-1,neqi+1] - 2. * dtm_ / r[2:nq-1]
						*( (1-phi) * F1D[2:nq-1,neqi+1] + phi * F1D[1:nq-2,neqi+1]   )
		end
	end
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function corrector!(nq, jj)

    nq::Int = nq
	lm_ = dtm_ / dr

	for i = 1:neqi
		U1D_c[2:nq-1,i] = U1D_p[2:nq-1,i] - lm_ *  (F1D[2:nq-1,i] - F1D[1:nq-2,i]) +
				dtm_ * ((1-phi) * G1D[2:nq-1,i] + phi * G1D[1:nq-2,i])
		if geom=="spherical"
			U1D_c[2:nq-1,i] = U1D_c[2:nq-1,i] - 2. * dtm_ ./ r[2:nq-1]  .*
						((1-phi) * (F1D[2:nq-1,i] - C1D[2:nq-1,i]) +
								phi * (F1D[1:nq-2,i] - C1D[1:nq-2,i]))
		end
	end

	#take care of electrons
	if mod(jj,2) == 1
		U1D_c[2:nq-1,neqi+1] = U1D_p[2:nq-1,neqi+1] - lm_ * (F1D[2:nq-1,neqi+1] - F1D[1:nq-2,neqi+1])  +
				dtm_ * ((1-phi) * G1D[2:nq-1,neqi+1] + phi * G1D[1:nq-2,neqi+1])
		if geom=="spherical" #add correction for spherical geometry
			U1D_c[2:nq-1,neqi+1] = U1D_c[2:nq-1,neqi+1] - 2. * dtm_ ./ r[2:nq-1]  .*
						((1-phi) * F1D[2:nq-1,neqi+1] + phi * F1D[1:nq-2,neqi+1])
		end
	else
		U1D_c[2:nq-1,neqi+1] = U1D_p[2:nq-1,neqi+1] - lm_ *  (F1D[3:nq,neqi+1] - F1D[2:nq-1,neqi+1]) +
				dtm_ * (phi * G1D[3:nq,neqi+1] + (1-phi) * G1D[2:nq-1,neqi+1])
		if geom=="spherical" #add correction for spherical geometry
			U1D_c[2:nq-1,neqi+1] = U1D_c[2:nq-1,neqi+1] - 2. * dtm_ ./ r[2:nq-1] .*
						(phi * F1D[3:nq,neqi+1] + (1-phi) * F1D[2:nq-1,neqi+1])
		end
	end
end

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function update_variables!(UU)
	# need to initialize electron density and speed to zero
	rho[:,nspec+1] = 0.
	ne = zeros(Float64, nz)
	u[:,nspec+1] = 0.

	for i = 1:nspec
		rho[:,i] = UU[:,3*(i-1)+1]
		u[:,i] = UU[:,3*(i-1)+2] ./ rho[:,i]
		T[:,i] = (g-1) * mi[i] * ( UU[:,3*(i-1)+3] ./ rho[:,i] - 0.5*u[:,i].^2 )
		p[:,i] = rho[:,i] .* T[:,i] / mi[i]

		F1D[:,3*(i-1)+1] = rho[:,i] .* u[:,i]
		F1D[:,3*(i-1)+2] = rho[:,i] .* u[:,i].^2 + p[:,i]
		F1D[:,3*(i-1)+3] = rho[:,i] .* u[:,i] .* UU[:,3*(i-1)+3] ./ rho[:,i] + p[:,i] .* u[:,i]

		if geom=="spherical"  # define correcting matrix
			C1D[:,3*(i-1)+1] = 0.
			C1D[:,3*(i-1)+2] = p[:,i]
			C1D[:,3*(i-1)+3] = 0.
		end

	    # electrons: quasi-neutrality + zero-current equation
	    rho[:,nspec+1] = rho[:,nspec+1] + me * Zi[i] * rho[:,i] / mi[i]  # quasi-neutrality
	    u[:,nspec+1] = u[:,nspec+1] + Zi[i] * rho[:,i] .* u[:,i] / mi[i]  # zero-current condition

		ni[:,i] = rho[:,i] / mi[i]
		ne = ne + Zi[i] * ni[:,i]
	end
	u[:,nspec+1] = u[:,nspec+1] ./ ne # correct with ne

	# electron temperature and pressure
	T[:,nspec+1] = (g-1) * me * ( UU[:,neqi+1] ./ rho[:,nspec+1] - 0.5*u[:,nspec+1].^2 )
	p[:,nspec+1] = rho[:,nspec+1] .* T[:,nspec+1] / me
	F1D[:,neqi+1] = u[:,nspec+1] .* UU[:,neqi+1]  + p[:,nspec+1] .* u[:,nspec+1]


    println("ne = ", ne)

	# limit max temperature by Tmax and avoid negative temperatures
	# for j = 1:nspec+1
	# 	T[:,j] = min( T[:,j], Tmax  )
	#     T[:,j] = max( T[:,j], 0.  ) # avoid negative temperature
	# end

end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
function source_terms!(nq)

    q_visc = zeros(Float64, nz)
    phi_visc = zeros(Float64, nz)

	if heattransfer_switch & electron_switch
		for i = 1:nspec+1
			for j = 1:nspec+1
				if j  != i
					Q_DT[:,i] = Q_DT[:,i] + k_DT[:,i,j] .* (T[:,j] -  T[:,i])
				end
			end
		end
		Q_DT[:,nspec+1] = Q_DT[:,nspec+1] + Qextra  #electrons
	elseif heattransfer_switch & ~electron_switch   #only ions
		for i = 1:nspec
			for j = 1:nspec
				if j != i
					Q_DT[:,i] = Q_DT[:,i] + k_DT[:,i,j] * ( T[:,j] -  T[:,i])
				end
			end
		end
	end

	if efield_switch & electron_switch
		Efield[2:nq-1] = - 0.5 * ( p[3:nq,nspec+1] - p[1:nq-2,nspec+1] ) /
		                  dr ./ (qe_C * ne[2:nq-1] )
	else
		Efield[:] = 0.
	end

	for i = 1:nspec
		G1D[:,3*(i-1)+2] = Zi[i] * qe_C * ni[:,i] .* Efield  #source term
		G1D[:,3*(i-1)+2] = G1D[:,3*(i-1)+2] + R1D[:,3*(i-1)+2]
		G1D[:,3*(i-1)+3] = Q_DT[:,i] + Zi[i] * qe_C * ni[:,i] .* Efield .* u[:,i]
		G1D[:,3*(i-1)+3] = G1D[:,3*(i-1)+3] + R1D[:,3*(i-1)+3]
	end

	if ion_viscosity

		 #find position of r = xmax_visc
 #		for i = 1, nz
 #			if (r[i] <= xmax_visc) then  #we found where r~=xmax_visc
 #				nn = i
 #				exit
 #			end
 #		end

		nn = maxloc(rho[:,1],1)

		for j = viscous_species:viscous_species #nspec
			G1D[nn+1:nq-2,3*(j-1)+3] = G1D[nn+1:nq-2,3*(j-1)+3] +
					0.25 * mu_ion[nn+1:nq-2,j] .* (u[nn+2:nq-1,j] - u[nn:nq-3,j]).^2 / dr^2 +
					mu_ion[nn+1:nq-2,j] .* (u[nn+1:nq-2,j] ./r[nn+1:nq-2]).^2	-
					0.5 * 2. * mu_ion[nn+1:nq-2,j] .* u[nn+1:nq-2,j] ./
                        r[nn+2:nq-1] .*(u[nn+2:nq-1,j] - u[nn:nq-3,j]) / dr
		end
	end



	q_diff[2:nq-1] = 1. / dr * ke[2:nq-1] .* (T[3:nq,nspec+1] - T[2:nq-1,nspec+1])

	vth_e = sqrt.(T[:,nspec+1] / me)
	q_FS[2:nq-1] = rho[3:nq,nspec+1] / me .* T[3:nq,nspec+1] .* vth_e[3:nq]  #ensures that it is forward-differencing later
	q_diff[:] = sign.(q_diff) .* min.(flimit * q_FS, abs.(q_diff))

	G1D[:,neqi+1] = Q_DT[:,nspec+1] - qe_C * ne .* Efield .* u[:,nspec+1]

	if electron_heat_conduction
		G1D[2:nq-1,neqi+1] = G1D[2:nq-1,neqi+1] + 1. / dr * (q_diff[2:nq-1] - q_diff[1:nq-2])
		if geom=="spherical"  #add correction due to spherical geometry
			G1D[:,neqi+1] = G1D[:,neqi+1] + 2. * q_diff ./ r
		end
	end

	 #Lindl, ch. 3 eq. (26) - rho in g/cm3, T in keV
	if bremsstrahlung
		G1D[:,neqi+1] = G1D[:,neqi+1] - 3.0e16 * (1.e3 * me * ne)^2 *
			             (1.e-3 * T[:,nspec+1] / qe_C)
	end
end
