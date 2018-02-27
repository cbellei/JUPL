using Constants

drx = 1.e-4

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
function calculate_collisions!(hydro, erf_table)

	#--- calculate collision coefficients ---
	#----------------------------------------
	hydro, xiab, k_DT, ke = collision_coefficients(hydro, erf_table, drx)

	if friction_switch
		hydro, Qextra = friction(hydro, xiab)
	end

	return hydro, k_DT, ke, Qextra
end

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
function collision_coefficients(hydro, erf_table, drx)

	pi = 2. * asin(1.)

	mu_ion = Array{Float64}(nz, nspec)

    xiab = zeros(Float64, nz, nspec, nspec)
    L_ab = zeros(Float64, nz, nspec, nspec)
    k_DT = zeros(Float64, nz, nspec+1, nspec+1)
    nu_DT = zeros(Float64, nz, nspec+1, nspec+1)

    mi_eV = Array{Float64}(nspec)
    muab = Array{Float64}(nspec, nspec)
    muab_eV = Array{Float64}(nspec, nspec)
    L_ie = Array{Float64}(nz, nspec)
    ni_cc = Array{Float64}(nz, nspec)
    ne_cc = Array{Float64}(nz)
    Te_eV = Array{Float64}(nz)
    taue = Array{Float64}(nz)
    taui = Array{Float64}(nz)
    Mab = Array{Float64}(nz)
    Tab = Array{Float64}(nz)
    Dab = Array{Float64}(nz)
    Qab = Array{Float64}(nz)
    meanL_ie = Array{Float64}(nz)
    T_eV = Array{Float64}(nz, nspec+1)

    erf = Array{Float64}(nz)
    i_erf = Array{Int64}(nz)

	for i = 1:nspec
		mi_eV[i] = 3.e8^2 * mi[i] / qe_C
		ni_cc[:,i] = 1.e-3 * hydro.rho[:,i] / mi_g[i]  #cm-3
	end
	ne_cc = 1.e-3 * hydro.rho[:,nspec+1] / me_g #cm-3

	for i = 1:nspec+1
	    T_eV[:,i] = max.(hydro.T[:,i] / qe_C, 0.) #avoid negative temperature
	end
	Te_eV = T_eV[:,nspec+1]

	for i = 1:nspec
		for j = 1:nspec
		    du = abs.(hydro.u[:,i] - hydro.u[:,j])
			muab[i,j] = 1.e-3 / ( mi_g[i] + mi_g[j] ) * mi_g[j] * mi_g[i]  #kg
			muab_eV[i,j] =  3.e8^2 / 1.6e-19 * muab[i,j]
			#The following L_ab is the same as the one used in LSP
			L_ab[:,i,j] = max.(Lab_min,  23. - log.(Zi[i]*Zi[j] *
			  sqrt.(ni_cc[:,i] * Zi[i]^2 ./ T_eV[:,i] +
                    ni_cc[:,j] * Zi[j]^2 ./ T_eV[:,j] ) ./
			  (muab_eV[i,j] * (T_eV[:,i]/mi_eV[i] +
			      T_eV[:,j]/mi_eV[j] + du.^2./3./3.e8^2 )))
				)
		end
	end

	#L_ab = dmax1(1.0, 23.0 - log( Z1*Z2*(m1_g + m2_g) / (m1_g*T2_eV + m2_g*T1_eV) * sqrt( n1_cc*Z1^2/T1_eV + n2_cc*Z2^2/T2_eV )  ));

	for i = 1:nspec #logLambda for ion-electrons
		L_ie[:,i] = max.(
            Lab_min,
            log.(3/2 * (Te_eV*1.6e-19*1.e7).^1.5 ./ sqrt.(pi * ne_cc) / (Zi[i]*qe^3)) #Krall&Trivelpiece
        )
	end

	for i = 1:nspec
		for j = 1:nspec
		    du = abs.(hydro.u[:,i] - hydro.u[:,j])
			Tab = ( mi_g[i] * hydro.T[:,j] + mi_g[j] * hydro.T[:,i] ) / ( mi_g[i] + mi_g[j] )
			Mab = sqrt.( 0.5 * muab[i,j] ./ Tab  ) .* du
			println(nspec, i, j)
			i_erf = min.(max.(1, floor.(Int, Mab / drx) + 1), 100000)
			for k = 1:nz
				erf[k] = erf_table[i_erf[k]] +
				(erf_table[i_erf[k]+1] - erf_table[i_erf[k]]) / drx *
                (Mab[k] - (i_erf[k]-1) * drx)
			end
			erf = max.(0., min.(erf, 1.0)) #so that's between 0. and 1.
			Qab = 1. / (32 * pi) * (Zi[i]*Zi[j]*qe_C^2/eps0 ./ Tab).^2 .* L_ab[:,i,j] *
                    3./2./Mab.^3 .*( 0.5 * sqrt(pi) * erf - Mab .* exp.(-Mab.^2) )

			Dab = 3. * pi / 16 * sqrt.( 2 * Tab / pi ./ muab[i,j] ) ./ Qab

			if check_frequency
				xiab[:,i,j] = mi[i] * 1.e6 * ni_cc[:,i] .*
								min.(1.e6 * ni_cc[:,j] .* Tab ./ Dab / mi[i], 1./(smax*dtm))  #attenzione mi(i) si elide####
			else
				xiab[:,i,j] = (1.e6 * ni_cc[:,i]) .* (1.e6 * ni_cc[:,j]) .* Tab ./ Dab
			end

			xiab[:,i,j] = max.(0., xiab[:,i,j]) #avoid negative values too
		end
	end

	for i = 1:nz
		meanL_ie[i] = sum(L_ie[i,:]) / nspec
	end
	taue = 3.44e5 * Te_eV.^1.5 ./ meanL_ie #NRL page 37	- excluding density because it cancels out anyways
	ke = 3.2 * (1.6e-12 * Te_eV) .* taue / me_g  #all cgs
	#multiply by 100 to get SI units ([ke] = [L^-1*T^-1])
	ke = 1.e2 * ke

	#calculate ion viscosity
	for i = 1:nspec
		#NRL page 38 - excluding density because it cancels out anyways
		taui =  2.09e7 * T_eV[:,i].^1.5 ./ L_ab[:,i,i] * sqrt(Ai[i])
		mu_ion[:,i] = 0.96 * ( 1.6e-12 * T_eV[:,i] ) .* taui #all cgs
		mu_ion[:,i] = 1.e-1 * mu_ion[:,i] # [mu_ion] = [M*L^-1*T^-1]
	end
	#Redefine ion viscosity to include the factor 4./3 in front of derivative
	mu_ion = 4./3 * mu_ion

	#coefficients for heat transfer (NRL page 34)
	for i = 1:nspec
		for j = 1:nspec
			if i != j
				nu_DT[:,i,j] = 1.8e-19 * sqrt(mi_g[i]) * sqrt(mi_g[j]) * Zi[i]^2 * Zi[j]^2 *
                    ni_cc[:,j] .* L_ab[:,i,j] ./ (mi_g[i] * T_eV[:,j] + mi_g[j] * T_eV[:,i]).^1.5
			end
		end
	end
	#heat transfer with electrons
	for i = 1:nspec
		nu_DT[:,i,nspec+1] = 1.8e-19 * sqrt(mi_g[i]) * sqrt(me_g) * Zi[i]^2 *
			ne_cc .* L_ie[:,i] ./ (mi_g[i] * Te_eV + me_g * T_eV[:,i]).^1.5
		nu_DT[:,nspec+1,i] = 1.8e-19 * sqrt(me_g) * sqrt(mi_g[i]) * Zi[i]^2 *
			ni_cc[:,i] .* L_ie[:,i] ./ (me_g * T_eV[:,i] + mi_g[i] * Te_eV).^1.5
	end


	if check_frequency
		for i = 1:nspec
			for j = 1:nspec
				if i != j
					nu_DT[:,i,j] = min.(nu_DT[:,i,j], 1.e-1/dtm)
				end
			end
		end
		#heat transfer with electrons
		for i = 1:nspec
			nu_DT[:,i,nspec+1] = min.(nu_DT[:,i,nspec+1], 1.e-1/dtm)
			nu_DT[:,nspec+1,i] = min.(nu_DT[:,nspec+1,i], 1.e-1/dtm)
		end
	end

	for i = 1:nspec
		for j = 1:nspec
			k_DT[:,i,j] = 1.e6 * 3./2 * ni_cc[:,i] .* nu_DT[:,i,j]   #SI units
		end
	end
	#with electrons now
	for i = 1:nspec
		k_DT[:,i,nspec+1] = 1.e6 * 3./2 * ni_cc[:,i] .* nu_DT[:,i,nspec+1]   #SI units
		k_DT[:,nspec+1,i] = 1.e6 * 3./2 * ne_cc .* nu_DT[:,nspec+1,i]   #SI units
	end

	return hydro, xiab, k_DT, ke
end


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
function friction(hydro, xiab)

	R1D = zeros(Float64, nz, neqi+1)
	Qextra = zeros(Float64, nz)

	for i = 1:nspec
		for j = 1:nspec
			if i != j
				R1D[:,3*(i-1)+2] = R1D[:,3*(i-1)+2] -
                        xiab[:,i,j] .* (hydro.u[:,i] - hydro.u[:,j])
    		end
    	end
    	R1D[:,3*(i-1)+3] = hydro.u[:,i] .* R1D[:,3*(i-1)+2]
    end

 	for i = 1:nspec
		Qextra = Qextra - R1D[:,3*(i-1)+3]
	end

	return hydro, Qextra
end
