using Constants

time1 = true
time2 = true
tm = 0.0

if ion_viscosity
    dtm_ = 0.5 * dtm
    eps_visc_max = 1./3 * eps_visc_max #artificial viscosity
else
    dtm_ = dtm
end

function init_spacegrid(dr)
    r = zeros(Float64, nz)

    if(geom=="spherical")
        r[1] = L + rmin
        r0[:,nregions+1] = r0[:,nregions+1] + rmin #correct also r0 array
    else
        r[1] = 0.
    end

    for k = 2:nz
        r[k] = r[k-1] + dr #note that dr is negative
    end

    if(geom=="spherical")
        r[nz] = rmin
    end

    return r
end


function init_art_viscosity(geom)
    eps_visc = zeros(Float64, nz, neqi+1)

    if geom=="slab"
		for k = 1:nz
		   eps_visc[k,1:neqi+1] = eps_visc_max
		end
	else #spherical
		y1 = log10(eps_visc_max)
		y2 = log10(eps_visc_max / eps_compress)
		x1 = log10(r[nz])
		x2 = log10(r[1])
		b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
		a  = ( y1 - b ) / x1
		for k = 1:nz
			eps_visc[k,1:neqi+1] = r[k]^a * 10^b
		end
	end

    return eps_visc
end


function init_variables(hydro)

    T0 = zeros(Float64, nz, nspec)
    N0 = zeros(Float64, nz, nspec)
    V0 = zeros(Float64, nz, nspec)

	for j = 1:nspec
		for m = 1:nregions
			for k = 1:nz
				if (r[k] >= r0[j,m]) & (r[k] <= r0[j,m+1])
					T0[k,j] = temp0[j,m]
					N0[k,j] = den0[j,m]
					V0[k,j] = vel0[j,m]
				end
			end
		end
	end

    for k = 1:nz
    	for j = 1:nspec
			hydro.rho[k,j] = N0[k,j] * mi[j]
			hydro.u[k,j] = V0[k,j]
			hydro.p[k,j] = hydro.rho[k,j] / mi[j] * T0[k,j]
			hydro.T[k,j] = T0[k,j]
		end
    end

    return hydro
end

function init_predictor_corrector(hydro)
	hydro.U1D_p = copy(hydro.U1D)
	hydro.U1D_c = copy(hydro.U1D)
    return hydro
end


#----------------------------------------------------
function do_smoothing(hydro)

	kgrad = Array{Int32}(2)
    nA, nB, uA, uB, TA, TB = [Array{Float64}(2) for i = 1:6]
    sgn1 = [0, 1]
    sgn2 = [1, 0]
	ε = 1.e-10
	vel = Array{Float64}(nz)
	mm = Array{Int64}(nedges)

	if smoothing
		#smooth out interface between nz0 and nz0+nsmooth
		#----------------------------------------------------
		if smooth_vars[1]==1
			println("density smoothing function")
			println("--------------------------")
			#for each ion species, find position of max and min gradients, then smooth out density
			for j = 1:nspec
				maxgrad = 0.
				mingrad = 1.e20
				for k = 1:nz-1
					grad = hydro.rho[k,j] / hydro.rho[k+1,j]
					if grad > maxgrad
						maxgrad = grad
						kgrad[1] = k
						nA[1] = hydro.rho[k,j] / mi[j]
						nB[1] = hydro.rho[k+1,j] / mi[j]
					end
					if grad < mingrad
						mingrad = grad
						kgrad[2] = k
						nA[2] = hydro.rho[k,j] / mi[j]
						nB[2] = hydro.rho[k+1,j] / mi[j]
					end
				end

				for k = 1:nedges
					if smooth_type[k,j]=="max"
						mm[k] = 1
					elseif smooth_type[k,j]=="min"
						mm[k] = 2
					end
				end

				for jj = 1:nedges
					m = mm[jj]
					for k = kgrad[m]-sgn1[m]*nsmooth : kgrad[m]+sgn2[m]*nsmooth
							#now smooth out interface
							y1 = log10(nA[m])
							y2 = log10(nB[m])
							x1 = log10(r[kgrad[m] - sgn1[m]*nsmooth]) #note r is in micron
							x2 = log10(r[kgrad[m] + sgn2[m]*nsmooth]) #note r is in micron
							b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
							a  = ( y1 - b ) / x1

							hydro.rho[k,j] = 10^(a*log10(r[k]) + b) * mi[j]
					end
				end
			end
		end

		if smooth_vars[2]==1
			println("velocity smoothing function")
			println("--------------------------")
			#for each ion species, find position of max and min gradients, then smooth out density
            ok = [0, 0]
			for j = 1:nspec
				maxgrad = 0.
				mingrad = 1.e20
				vel = max(ε, abs(hydro.u[:,j]))
				for k = 1:nz-1
					grad = vel[k] / vel[k+1]
					if grad > maxgrad
						maxgrad = grad
						kgrad[1] = k
						uA[1] = vel[k]
						uB[1] = vel[k+1]
						if (k+tsmooth<nz) & (k-tsmooth>0)
							ok[1] = 1
						end
					end
					if grad < mingrad
						mingrad = grad
						kgrad[2] = k
						uA[2] = vel[k]
						uB[2] = vel[k+1]
						if (k+tsmooth < nz) & (k-tsmooth > 0)
							ok[2] = 1
						end
					end
				end

				for k = 1:nedges
					if smooth_type[k,j]=="max"
						mm[k]=1
					elseif smooth_type[k,j]=="min"
						mm[k]=2
					end
				end

				for jj = 1:nedges
					m = mm[jj]
					if ok[m].ne.1
                        continue #because it is a bad point...
                    end

					for k = kgrad[m]-sgn1[m]*vsmooth : kgrad[m]+sgn2[m]*vsmooth
							#now smooth out interface
							y1 = log10(uA[m])
							y2 = log10(uB[m])
							x1 = log10(r[kgrad[m] - sgn1[m]*vsmooth]) #note r is in micron
							x2 = log10(r[kgrad[m] + sgn2[m]*vsmooth]) #note r is in micron
							b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
							a  = ( y1 - b ) / x1

							hydro.u[k,j] = - (10^( a*log10(r[k]) + b ))
					end

				end
			end
		end

		if smooth_vars[3]==1
			println("temperature smoothing function")
			println("------------------------------")
			#for each ion species, find position of max and min gradients, then smooth out density
			ok = [0, 0]
			for j = 1:nspec
				maxgrad = 0.
				mingrad = 1.e20
				for k = 1:nz-1
					grad = hydro.T[k,j] / hydro.T[k+1,j]
					if grad > maxgrad
						maxgrad = grad
						kgrad[1] = k
						TA[1] = hydro.T[k,j]
						TB[1] = hydro.T[k+1,j]
						if (k+tsmooth < nz) & (k-tsmooth > 0)
							ok[1] = 1
						end
					end
					if grad < mingrad
						mingrad = grad
						kgrad[2] = k
						TA[2] = hydro.T[k,j]
						TB[2] = hydro.T[k+1,j]
						if (k+tsmooth < nz) & (k-tsmooth > 0)
							ok[2] = 1
						end
					end
				end

				for k = 1:nedges
					if smooth_type[k,j]=="max"
						mm[k] = 1
					elseif smooth_type[k,j]=="min"
						mm[k] = 2
					end
				end

				for jj = 1:nedges
					m = mm[jj]
					if ok[m]!=1
                        continue #because it is a bad point...
                    end
					for k = kgrad[m]-sgn1[m]*tsmooth:kgrad[m] + sgn2[m]*tsmooth
							#now smooth out interface
							y1 = log10(TA[m])
							y2 = log10(TB[m])
							x1 = log10(r[kgrad[m] - sgn1[m]*tsmooth]) #note r is in micron
							x2 = log10(r[kgrad[m] + sgn2[m]*tsmooth]) #note r is in micron
							b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
							a  = ( y1 - b ) / x1

							hydro.T[k,j] = 10^(a*log10(r[k]) + b)
					end

				end
			end
		end

		#now update pressure after all this smoothing...
		for j = 1:nspec
			hydro.p[:,j] = hydro.rho[:,j] / mi[j] .* hydro.T[:,j]
		end
	end
    return hydro
end


#-----------------------------------
function init_fluxes(hydro)

	println()
	println(" WARNING - forcing quasi-neutrality and zero-current conditions")
	println()

	hydro.rho[:,nspec+1] = 0.
	hydro.u[:,nspec+1] = 0.

	#electrons
    for k = 1:nz
    	for j = 1:nspec
			hydro.rho[k,nspec+1] = hydro.rho[k,nspec+1] + me * Zi[j] * hydro.rho[k,j] / mi[j] #quasi-neutrality
			hydro.u[k,nspec+1] = hydro.u[k,nspec+1] + Zi[j] * hydro.rho[k,j] * hydro.u[k,j] / mi[j]   #zero-current condition
    	end
    	if !restart #take species 1 as reference
			hydro.p[k,nspec+1] = hydro.rho[k,nspec+1] / me * hydro.T[k,1] #electron pressure
			hydro.T[k,nspec+1] = hydro.T[k,1] #electron temperature
		else #we have read the electron temperature from file
			hydro.p[k,nspec+1] = hydro.rho[k,nspec+1] / me * hydro.T[k,nspec+1] #electron pressure
		end
    end

#   correct for electron velocity
    hydro.u[:,nspec+1] = hydro.u[:,nspec+1] ./ (hydro.rho[:,nspec+1] / me)
#	apply BC at piston side
    for j = 1:nspec+1
        hydro.u[1, j] = 0.
    end

#	finally, calculate U1D and F1D for all species
	for j = 1:nspec
		hydro.U1D[:, 3*(j-1)+1] = hydro.rho[:,j]
		hydro.U1D[:, 3*(j-1)+2] = hydro.rho[:,j] .* hydro.u[:,j]
		hydro.U1D[:, 3*(j-1)+3] = hydro.p[:,j] / (g-1) + 0.5 * hydro.rho[:,j] .* hydro.u[:,j].^2
		hydro.F1D[:, 3*(j-1)+1] = hydro.rho[:,j] .* hydro.u[:,j]
		hydro.F1D[:, 3*(j-1)+2] = hydro.rho[:,j] .* hydro.u[:,j].^2 + hydro.p[:,j]
		hydro.F1D[:, 3*(j-1)+3] = hydro.u[:,j] .* ( g / (g-1) * hydro.p[:,j] + 0.5 * hydro.rho[:,j] .* hydro.u[:,j].^2 )
	end
	#electrons
	hydro.U1D[:,neqi+1] = hydro.p[:,nspec+1] ./ (g-1) + 0.5 * hydro.rho[:,nspec+1] .* hydro.u[:,nspec+1].^2
	hydro.F1D[:,neqi+1] = hydro.u[:,nspec+1] .* ( g / (g-1) * hydro.p[:,nspec+1] + 0.5 * hydro.rho[:,nspec+1] .* hydro.u[:,nspec+1].^2 )

    return hydro
end
