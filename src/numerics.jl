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
	if mod(jj,2)==1
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

	# limit max temperature by Tmax and avoid negative temperatures
	# for j = 1:nspec+1
	# 	T[:,j] = min( T[:,j], Tmax  )
	#     T[:,j] = max( T[:,j], 0.  ) # avoid negative temperature
	# end

end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
