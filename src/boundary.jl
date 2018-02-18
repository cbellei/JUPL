function apply_BC!()

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
					U1D[nn,3*(i-1)+1] = U1D[nn1,3*(i-1)+1] - lm * F1D[nn1,3*(i-1)+1] #continuity
					U1D[nn,3*(i-1)+2] = 0. #momentum
					U1D[nn,3*(i-1)+3] = U1D[nn1,3*(i-1)+3] - lm * F1D[nn1,3*(i-1)+3] #energy
				end
				#now electrons
				U1D[nn,neqi+1] =  U1D[nn1,neqi+1] - lm * F1D[nn1,neqi+1]
			elseif geom=="spherical" #spherical geometry
				#boundary condition for reflection at wall side
				for i = 1:nspec
					U1D[nn,3*(i-1)+1] = U1D[nn1,3*(i-1)+1] - lm * F1D[nn1,3*(i-1)+1] 	&
								- 2. * dtm_ / r[nz] * phi * F1D[nz-1,3*(i-1)+1]
					U1D[nn,3*(i-1)+2] = 0. #momentum
					U1D[nn,3*(i-1)+3] = U1D[nn1,3*(i-1)+3] - lm * F1D[nn1,3*(i-1)+3]  &
								- 2. * dtm_ / r[nz] * phi * F1D[nn1,3*(i-1)+3]
				end
				#now electrons
				U1D[nz,neqi+1] =  U1D[nz-1,neqi+1]
			end
		elseif boundary[kk]=="open"
			#boundary condition for reflection at wall side
			for i = 1:nspec
				U1D[nn,3*(i-1)+1] = rho[nn,i]
				U1D[nn,3*(i-1)+2] = rho[nn,i] * u[nn,i] #momentum
				U1D[nn,3*(i-1)+3] = p[nn,i] / (g-1) + 0.5 * rho[nn,i] * u[nn,i]^2
			end
			#now electrons
			U1D[nn,neqi+1] =  U1D[nn,neqi+1]
		else
			println("wrong boundary condition at boundary #", kk)
		end

		#update variables
		for i = 1:nspec
			rho[nn,i] = U1D[nn,3*(i-1)+1]
			u[nn,i] = U1D[nn,3*(i-1)+2] / rho[nn,i]
			T[nn,i] = (g-1) * mi[i] * ( U1D[nn,3*(i-1)+3]/rho[nn,i] - 0.5*u[nn,i]^2 )
			p[nn,i] = rho[nn,i] * T[nn,i] / mi[i]
		end
	end

	#BC for predictor and corrector
	U1D_p[1,:] = U1D[1,:]
	U1D_c[1,:] = U1D[1,:]
	U1D_p[nz,:] = U1D[nz,:]
	U1D_c[nz,:] = U1D[nz,:]
end
