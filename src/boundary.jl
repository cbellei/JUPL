function apply_BC!(hydro)

	for kk = 1:2
		if kk==1
			nn = 1
			nn1 = 2
		else
			nn = nz
			nn1 = nz-1
		end
		if boundary[kk]=="reflection"
			if geom=="slab"
				#boundary condition for reflection at wall side
				for i = 1:nspec
					hydro.U1D[nn,3*(i-1)+1] = hydro.U1D[nn1,3*(i-1)+1] - lm * hydro.F1D[nn1,3*(i-1)+1] #continuity
					hydro.U1D[nn,3*(i-1)+2] = 0. #momentum
					hydro.U1D[nn,3*(i-1)+3] = hydro.U1D[nn1,3*(i-1)+3] - lm * hydro.F1D[nn1,3*(i-1)+3] #energy
				end
				#now electrons
				hydro.U1D[nn,neqi+1] =  hydro.U1D[nn1,neqi+1] - lm * hydro.F1D[nn1,neqi+1]
			elseif geom=="spherical" #spherical geometry
				#boundary condition for reflection at wall side
				for i = 1:nspec
					hydro.U1D[nn,3*(i-1)+1] = hydro.U1D[nn1,3*(i-1)+1] - lm * hydro.F1D[nn1,3*(i-1)+1] 	&
								- 2. * dtm_ / r[nz] * phi * hydro.F1D[nz-1,3*(i-1)+1]
					hydro.U1D[nn,3*(i-1)+2] = 0. #momentum
					hydro.U1D[nn,3*(i-1)+3] = hydro.U1D[nn1,3*(i-1)+3] - lm * hydro.F1D[nn1,3*(i-1)+3]  &
								- 2. * dtm_ / r[nz] * phi * hydro.F1D[nn1,3*(i-1)+3]
				end
				#now electrons
				hydro.U1D[nz,neqi+1] =  hydro.U1D[nz-1,neqi+1]
			end
		elseif boundary[kk]=="open"
			#boundary condition for reflection at wall side
			for i = 1:nspec
				hydro.U1D[nn,3*(i-1)+1] = hydro.rho[nn,i]
				hydro.U1D[nn,3*(i-1)+2] = hydro.rho[nn,i] * hydro.u[nn,i] #momentum
				hydro.U1D[nn,3*(i-1)+3] = hydro.p[nn,i] / (g-1) + 0.5 * hydro.rho[nn,i] * hydro.u[nn,i]^2
			end
			#now electrons
			hydro.U1D[nn,neqi+1] =  hydro.U1D[nn,neqi+1]
		else
			println("wrong boundary condition at boundary #", kk)
		end

		#update variables
		for i = 1:nspec
			hydro.rho[nn,i] = hydro.U1D[nn,3*(i-1)+1]
			hydro.u[nn,i] = hydro.U1D[nn,3*(i-1)+2] / hydro.rho[nn,i]
			hydro.T[nn,i] = (g-1) * mi[i] * ( hydro.U1D[nn,3*(i-1)+3]/hydro.rho[nn,i] - 0.5*hydro.u[nn,i]^2 )
			hydro.p[nn,i] = hydro.rho[nn,i] * hydro.T[nn,i] / mi[i]
		end
	end

	#BC for predictor and corrector
	hydro.U1D_p[1,:] = hydro.U1D[1,:]
	hydro.U1D_c[1,:] = hydro.U1D[1,:]
	hydro.U1D_p[nz,:] = hydro.U1D[nz,:]
	hydro.U1D_c[nz,:] = hydro.U1D[nz,:]

	return hydro
end
