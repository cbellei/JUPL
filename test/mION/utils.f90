subroutine allocate_arrays()

	use constants	
	implicit none

	allocate( r(nz), xi(nz), p(nz,nspec+1), u(nz,nspec+1), rho(nz,nspec+1), T(nz,nspec+1), q_FS(nz), q_diff(nz), xi12(nz), &
		nu0(nz), ke(nz), du(nz), vth_e(nz), Te(nz), Efield(nz), eps_visc(nz,neqi+1) )
	allocate( U1D(nz,neqi+1), F1D(nz,neqi+1), R1D(nz,neqi+1), U1D_p(nz,neqi+1), U1D_c(nz,neqi+1), G1D(nz,neqi+1) )
	allocate( C1D(nz,neqi), mu_ion(nz,nspec) )
	allocate( ne(nz), ni(nz,nspec) )	
	allocate( k_DT(nz,nspec+1,nspec+1), Q_DT(nz,nspec+1), Qextra(nz), xiab(nz,nspec,nspec) )
	allocate( L_ab(nz,nspec,nspec)  )

	allocate( mi(nspec), mi_g(nspec), qi(nspec), Zi(nspec), Ai(nspec) )	
	allocate( nleft(nspec), nmiddle(nspec), nright(nspec), r_shell(nspec), dr_shell(nspec) )
	allocate( r0(nspec,nregions+1), den0(nspec,nregions), temp0(nspec,nregions), vel0(nspec,nregions) )
	allocate( T0(nz,nspec), N0(nz,nspec), V0(nz,nspec) )
	allocate( LDIAG(nz-3), DIAG(nz-2), UDIAG(nz-3), BB(nz-2) )
	allocate( AA(nz-2,nz-2), IPIV(nz-2) )
	allocate( LDIAG2(nz-3), DIAG2(nz-2), UDIAG2(nz-3), BB2(nz-2) )
	allocate( AA2(nz-2,nz-2), IPIV2(nz-2) )
	
	allocate( smooth_type(nedges,nspec) )
		
end subroutine allocate_arrays


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine predictor(nq,jj)

	use constants
	implicit none
    integer, intent(in) :: nq, jj
	integer :: i, j, nn
	real*8 :: lm_
	real*8, dimension(nz) :: phi_visc
	
	lm_ = dtm_ / dr

	!take care of ions
	do i = 1, neqi
		U1D_p(2:nq-1,i) = U1D(2:nq-1,i) - lm_ * ( F1D(3:nq,i) - F1D(2:nq-1,i) ) &
				+ dtm_ * ( phi * G1D(3:nq,i) + (1-phi) * G1D(2:nq-1,i) )
		if (geom=="spherical") then !add correction for spherical geometry
			U1D_p(2:nq-1,i) = U1D_p(2:nq-1,i) - 2. * dtm_ / r(2:nq-1)  & 
						* ( phi * ( F1D(3:nq,i) - C1D(3:nq,i) ) &
							+ (1-phi) * ( F1D(2:nq-1,i) - C1D(2:nq-1,i) )  )
		endif
	enddo

	!take care of electrons
	if ( mod(jj,2)==1 ) then 
		U1D_p(2:nq-1,neqi+1) = U1D(2:nq-1,neqi+1) - lm_ * ( F1D(3:nq,neqi+1) - F1D(2:nq-1,neqi+1) ) &
				+ dtm_ * ( phi * G1D(3:nq,neqi+1) + (1-phi) * G1D(2:nq-1,neqi+1) )
		if (geom=="spherical") then !add correction for spherical geometry
			U1D_p(2:nq-1,neqi+1) = U1D_p(2:nq-1,neqi+1) - 2. * dtm_ / r(2:nq-1)  & 
						* ( phi * F1D(3:nq,neqi+1) + (1-phi) * F1D(2:nq-1,neqi+1)   )	
		endif
	else
		U1D_p(2:nq-1,neqi+1) = U1D(2:nq-1,neqi+1) - lm_ * ( F1D(2:nq-1,neqi+1) - F1D(1:nq-2,neqi+1) ) &
				+ dtm_ * ( (1-phi) * G1D(2:nq-1,neqi+1) + phi * G1D(1:nq-2,neqi+1) )
		if (geom=="spherical") then !add correction for spherical geometry
			U1D_p(2:nq-1,neqi+1) = U1D_p(2:nq-1,neqi+1) - 2. * dtm_ / r(2:nq-1)  & 
						* ( (1-phi) * F1D(2:nq-1,neqi+1) + phi * F1D(1:nq-2,neqi+1)   )	
		endif
	endif
			
end subroutine predictor


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine corrector(nq,jj)

	use constants
	implicit none
    integer, intent(in) :: nq, jj
	integer :: i,j,nn
	real*8 :: lm_
	real*8, dimension(nz) :: phi_visc
	
	lm_ = dtm_ / dr
	


	do i = 1, neqi
		U1D_c(2:nq-1,i) = U1D_p(2:nq-1,i) - lm_ *  ( F1D(2:nq-1,i) - F1D(1:nq-2,i) ) &
				+ dtm_ * ( (1-phi) * G1D(2:nq-1,i) + phi * G1D(1:nq-2,i) ) 
		if (geom=="spherical") then
			U1D_c(2:nq-1,i) = U1D_c(2:nq-1,i) - 2. * dtm_ / r(2:nq-1)  &
						*   ( (1-phi) * ( F1D(2:nq-1,i) - C1D(2:nq-1,i) ) &
								+ phi * ( F1D(1:nq-2,i) - C1D(1:nq-2,i) )  )	
		endif
	enddo 

	!take care of electrons
	if ( mod(jj,2)==1 ) then		
		U1D_c(2:nq-1,neqi+1) = U1D_p(2:nq-1,neqi+1) - lm_ * ( F1D(2:nq-1,neqi+1) - F1D(1:nq-2,neqi+1) )  &
				+ dtm_ * ( (1-phi) * G1D(2:nq-1,neqi+1) + phi * G1D(1:nq-2,neqi+1) )
		if (geom=="spherical") then !add correction for spherical geometry
			U1D_c(2:nq-1,neqi+1) = U1D_c(2:nq-1,neqi+1) - 2. * dtm_ / r(2:nq-1)  & 
						* ( (1-phi) * F1D(2:nq-1,neqi+1) + phi * F1D(1:nq-2,neqi+1)   )	
		endif
	else
		U1D_c(2:nq-1,neqi+1) = U1D_p(2:nq-1,neqi+1) - lm_ *  ( F1D(3:nq,neqi+1) - F1D(2:nq-1,neqi+1) )  &
				+ dtm_ * ( phi * G1D(3:nq,neqi+1) + (1-phi) * G1D(2:nq-1,neqi+1) )
		if (geom=="spherical") then !add correction for spherical geometry
			U1D_c(2:nq-1,neqi+1) = U1D_c(2:nq-1,neqi+1) - 2. * dtm_ / r(2:nq-1)  & 
						* ( phi * F1D(3:nq,neqi+1) + (1-phi) * F1D(2:nq-1,neqi+1)   )	
		endif
	endif	
	
end subroutine corrector


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine calculate_viscous_term(nq,jj)

	use constants
	
	implicit none
    integer, intent(in) :: nq, jj
    real*8 :: rr, m1, m2, m3, mu_ion_ip12, mu_ion_im12, r_ip12, r_im12, muold_ip12, muold_im12
    real*8 :: muold_ion_ip12, muold_ion_im12
    real*8 ::  u_ip12, u_im12, uoldip12, uoldim12
    real*8, dimension(nz) :: uold, Tnew, tion, TeV, Tdiff, udiff, phi_visc, UU
    real*8, dimension(nz) :: taui, U1Dnew, T_before, muold
    integer :: i, j, k, m, count, INFO
    real*8 :: CC, res, m1old, m2old, m3old, lambda
    logical :: is_diag_dominant, write_to_file=.true., first_iteration, abort
    logical :: first_try
    integer :: nn, nnn, minE, w
    real :: toll1 = 1.e-4, toll2 = 1.e-10
    
	!find position of r = xmax_visc
!	do i = 1, nz
!		if (r(i) <= xmax_visc) then !we found where r~=xmax_visc
!			nn = i
!			exit
!		endif
!	enddo
	nn = maxloc(rho(:,1),1)		
    nnn = nq-1-nn

	!Define rr
	rr = dtm / dr**2 / nsubcycle

	!Solve Crank-Nicolson system	
	do j = viscous_species, viscous_species!nspec		

		uold = u(:,j)
		muold = mu_ion(:,j)
		res = 1.e5
		count = 1		

		!MOMENTUM EQUATION: Solution of Crank-Nicolson system
		!----------------------------------------------------
		BB = 0.
		AA = 0.
		UU = 0.
		phi_visc = 0.
		is_diag_dominant = .true.
		do i = nn, nq - 2 !Because AA has size (nz-2)x(nz-2)
			mu_ion_ip12 = 0.5 * (  mu_ion(i+1,j) + mu_ion(i+2,j) )
			mu_ion_im12 = 0.5 * (  mu_ion(i,j) + mu_ion(i+1,j) )
			muold_ion_ip12 = 0.5 * (  muold(i+1) + muold(i+2) )
			muold_ion_im12 = 0.5 *  (  muold(i) + muold(i+1) )

			if (geom=="spherical") then
				r_ip12 = 0.5 * ( r(i+2) + r(i+1) )
				r_im12 = 0.5 * ( r(i+1) + r(i) )
				m1 = mu_ion_ip12*r_ip12**2
				m1old = muold_ion_ip12*r_ip12**2
				m2 = mu_ion_ip12*r_ip12**2 + mu_ion_im12*r_im12**2
				m2old = muold_ion_ip12*r_ip12**2 + muold_ion_im12*r_im12**2
				m3 = mu_ion_im12*r_im12**2
				m3old = muold_ion_im12*r_im12**2			
			else !slab
				m1 = mu_ion_ip12	!==mu_ion@i+1/2	
				m2 = ( mu_ion_ip12 + mu_ion_im12 )
				m3 = mu_ion_im12	  !==mu_ion@i-1/2			
			endif	 

			DIAG(i) = 2.*r(i+1)**2*rho(i+1,j) + rr * m2   &
				+ r(i+1)*dtm * ( 2.*mu_ion(i+1,j)/r(i+1) + 0.5*(mu_ion(i+1,j)-mu_ion(i-1,j))/dr  )
			AA(i,i) = DIAG(i)
			if (i < nq-2) then
				UDIAG(i) = - rr * m1
				AA(i,i+1) = UDIAG(i)
			endif
			if (i > nn) then
				LDIAG(i-1) = - rr * m3
				AA(i,i-1) = LDIAG(i-1)
			endif
			if (i == nn) then !we are at R = L -> add B.C. with u/=0
				BB(i) = (2.*r(i+1)**2*rho(i+1,j) - rr*m2old)*uold(i+1)    &
						 + rr*m1old*uold(i+2) + 2.*rr*m3old*uold(i)
			else if (i == nq-2) then !we are at bottom of array
				BB(i) = rr*m3old*uold(i) + (2.*r(i+1)**2*rho(i+1,j) - rr*m2old)*uold(i+1)				 
			else
				BB(i) = rr*m3old*uold(i) + (2.*r(i+1)**2*rho(i+1,j) - rr*m2old )*uold(i+1)   &
						+ rr*m1old*uold(i+2)  &
						- r(i+1)*dtm*uold(i+1) * ( 2.*muold(i+1)/r(i+1) + 0.5*(muold(i+1) - muold(i-1))/dr )
			endif

			if(is_diag_dominant) then
				if ( (i>nn) .and. (i<nq-2) .and. &
						 ( abs(DIAG(i)) < abs(UDIAG(i))+abs(LDIAG(i-1)) ) ) then
					is_diag_dominant = .false.
				elseif ( (i==nn) .and. ( abs(DIAG(i)) < abs( UDIAG(i))) ) then
					is_diag_dominant = .false.
				elseif ( (i==nq-2) .and. ( abs(DIAG(i)) < abs( LDIAG(i-1)) ) ) then
					is_diag_dominant = .false.
				endif
			endif
			
		enddo !do i = 1, nq - 2	
			

			if(is_diag_dominant) then !matrix is diagonally dominant -> use Thomas' algorithm
				!Thomas' algorithm
				call DGTSV(nnn,1,LDIAG(nn:nz-3),DIAG(nn:nz-2),UDIAG(nn:nz-3),BB(nn:nz-2),nnn,INFO)	
			else
				!LU decomposition
				call DGESV(nnn,1,AA(nn:nz-2,nn:nz-2),nnn,IPIV(nn:nz-2),BB(nn:nz-2),nnn,INFO)
			endif	


		!IF COMMENTED, DO NOT SOLVE MOMENTUM EQUATION
		u(nn+1:nq-1,j) = BB(nn:nz-2)

		U1D(nn+1:nq-1,3*(j-1)+2) = rho(nn+1:nq-1,j) * u(nn+1:nq-1,j)
	enddo !j=1,nspec

	!Add artificial viscosity
	U1D(2:nq-1,1:neqi+1) = U1D(2:nq-1,1:neqi+1)   &
					+ eps_visc(1:nq-2,1:neqi+1) / nsubcycle * ( U1D(3:nq,1:neqi+1) - 2*U1D(2:nq-1,1:neqi+1) + U1D(1:nq-2,1:neqi+1)  )


end subroutine calculate_viscous_term



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine source_terms(nq)

	use constants
	
    implicit none    
    integer, intent(in) :: nq
    integer :: i, j, m, nn
    real*8, dimension(nz) :: q_visc, phi_visc

    Q_DT = 0.
    G1D = 0.
    q_visc = 0.
    
	if (heattransfer_switch .and. electron_switch) then
		do i = 1, nspec+1
			do j = 1, nspec+1
				if (j .ne. i) then
					Q_DT(:,i) = Q_DT(:,i) + k_DT(:,i,j) * ( T(:,j) -  T(:,i))
				endif
			enddo
		enddo
		Q_DT(:,nspec+1) = Q_DT(:,nspec+1) + Qextra !electrons
	elseif( heattransfer_switch .and. (.not.electron_switch) ) then  !only ions
		do i = 1, nspec
			do j = 1, nspec
				if (j .ne. i) then
					Q_DT(:,i) = Q_DT(:,i) + k_DT(:,i,j) * ( T(:,j) -  T(:,i))
				endif
			enddo
		enddo
	endif
	
	if (Efield_switch .and. electron_switch ) then
		Efield(2:nq-1) = - 0.5 * ( p(3:nq,nspec+1) - p(1:nq-2,nspec+1) ) / dr  &
							 / (qe_C * ne(2:nq-1) ) 
	else
		Efield = 0.
	endif
	

	do i = 1, nspec
		G1D(:,3*(i-1)+2) = Zi(i) * qe_C * ni(:,i) * Efield !source term
		G1D(:,3*(i-1)+2) = G1D(:,3*(i-1)+2) + R1D(:,3*(i-1)+2)
		G1D(:,3*(i-1)+3) = Q_DT(:,i) + Zi(i) * qe_C * ni(:,i) * Efield * u(:,i)
		G1D(:,3*(i-1)+3) = G1D(:,3*(i-1)+3) + R1D(:,3*(i-1)+3)
	enddo

	if (ion_viscosity) then

		!find position of r = xmax_visc
!		do i = 1, nz
!			if (r(i) <= xmax_visc) then !we found where r~=xmax_visc
!				nn = i
!				exit
!			endif
!		enddo

		nn = maxloc(rho(:,1),1)

		do j = viscous_species, viscous_species!nspec								
			G1D(nn+1:nq-2,3*(j-1)+3) = G1D(nn+1:nq-2,3*(j-1)+3) 	&
					+  0.25*mu_ion(nn+1:nq-2,j)*(u(nn+2:nq-1,j)-u(nn:nq-3,j))**2/dr**2 &
					+  mu_ion(nn+1:nq-2,j)*(u(nn+1:nq-2,j)/r(nn+1:nq-2))**2		&
					-  0.5*2.*mu_ion(nn+1:nq-2,j)*u(nn+1:nq-2,j)/r(nn+2:nq-1)*(u(nn+2:nq-1,j)-u(nn:nq-3,j))/dr					
		enddo
	endif



	q_diff(2:nq-1) = 1. / dr * ke(2:nq-1) * ( T(3:nq,nspec+1) - T(2:nq-1,nspec+1) ) 

	vth_e = sqrt( T(:,nspec+1) / me  )
	q_FS(2:nq-1) = rho(3:nq,nspec+1) / me * T(3:nq,nspec+1) * vth_e(3:nq) !ensures that it is forward-differencing later
	q_diff = sign(  min( flimit * q_FS, abs(q_diff) ), q_diff)

	G1D(:,neqi+1) = Q_DT(:,nspec+1) - qe_C * ne * Efield * u(:,nspec+1)
	
	if (electron_heat_conduction) then
		G1D(2:nq-1,neqi+1) = G1D(2:nq-1,neqi+1) + 1. / dr * ( q_diff(2:nq-1) - q_diff(1:nq-2) )
		if (geom=="spherical") then !add correction due to spherical geometry
			G1D(:,neqi+1) = G1D(:,neqi+1) + 2. * q_diff / r(:)		
		endif
	endif

	!Lindl, ch. 3 eq. (26) - rho in g/cm3, T in keV
	if (bremsstrahlung) then
		G1D(:,neqi+1) = G1D(:,neqi+1) - 3.0e16 * (1.e3*me*ne)**2 	&
			* (1.e-3*T(:,nspec+1)/qe_C)
	endif
	
		
end subroutine source_terms

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine update_variables(UU)

	use constants
	
	implicit none
	real*8, intent(in), dimension(nz,neqi+1) :: UU
	integer :: i,j   
	    
	!need to initialize electron density and speed to zero
	rho(:,nspec+1) = 0. 
	ne = 0.
	u(:,nspec+1) = 0.
	

	do i = 1, nspec
		rho(:,i) = UU(:,3*(i-1)+1)
		u(:,i) = UU(:,3*(i-1)+2) / rho(:,i)
		T(:,i) = (g-1) * mi(i) * ( UU(:,3*(i-1)+3)/rho(:,i) - 0.5*u(:,i)**2 )
		p(:,i) = rho(:,i) * T(:,i) / mi(i)

		F1D(:,3*(i-1)+1) = rho(:,i) * u(:,i)
		F1D(:,3*(i-1)+2) = rho(:,i) * u(:,i)**2 + p(:,i)
		F1D(:,3*(i-1)+3) = rho(:,i) * u(:,i) * UU(:,3*(i-1)+3) / rho(:,i) + p(:,i) * u(:,i)

		if (geom=="spherical") then !define correcting matrix
			C1D(:,3*(i-1)+1) = 0.
			C1D(:,3*(i-1)+2) = p(:,i)
			C1D(:,3*(i-1)+3) = 0.	
		endif	

	    !electrons: quasi-neutrality + zero-current equation
	    rho(:,nspec+1) = rho(:,nspec+1) + me * Zi(i) * rho(:,i) / mi(i)  !quasi-neutrality
	    u(:,nspec+1) = u(:,nspec+1) + Zi(i) * rho(:,i) * u(:,i) / mi(i)  !zero-current condition  

		ni(:,i) = rho(:,i) / mi(i)
		ne(:) = ne(:) + Zi(i) * ni(:,i)
	enddo
	u(:,nspec+1) = u(:,nspec+1) / ne !correct with ne
	
	!electron temperature and pressure											
	T(:,nspec+1) = (g-1) * me * ( UU(:,neqi+1)/rho(:,nspec+1) - 0.5*u(:,nspec+1)**2 )				
	p(:,nspec+1) = rho(:,nspec+1) * T(:,nspec+1) / me	
	F1D(:,neqi+1) = u(:,nspec+1) * UU(:,neqi+1)  + p(:,nspec+1) * u(:,nspec+1)
	
	!limit max temperature by Tmax and avoid negative temperatures
	do j = 1, nspec+1
		T(:,j) = min( T(:,j), Tmax  )
	    T(:,j) = max( T(:,j), 0.  ) !avoid negative temperature
	enddo
	
end subroutine update_variables

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine 	calculate_collisions(drx,erf_table)

	use constants
	
	implicit none	
	real*8, intent(in) :: drx
	real*8, intent(in), dimension(100001) :: erf_table
		
	!--- calculate collision coefficients ---
	!----------------------------------------							
	call collision_coefficients( erf_table, drx )

	if (friction_switch) then
		call friction()
	endif
	!----------------------------------------

end subroutine calculate_collisions


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine collision_coefficients(erf_table, drx)

	use constants

	implicit none
	real*8, intent(in) :: drx
	real*8, intent(in), dimension(100001) :: erf_table
	real*8, dimension(nspec,nspec) :: muab, muab_eV
	real*8, dimension(nspec) :: mi_eV
	real*8, dimension(nz,nspec) :: L_ie
	real*8, dimension(nz) :: meanL_ie, ne_cc, Te_eV, taue, taui, Mab, Tab, Dab, Qab
	real*8, dimension(nz,nspec+1,nspec+1) :: nu_DT
	real*8, dimension(nz,nspec) :: ni_cc
	real*8, dimension(nz,nspec+1) :: T_eV
	real*8 :: erf(nz), pi 
	integer :: i, j, k, i_erf(nz)
		
	pi = 2. * asin(1.)
	
	nu_DT = 0.
	k_DT = 0.
	xiab = 0.
			
	do i = 1, nspec
		mi_eV(i) = 3.e8**2 * mi(i) / qe_C
		ni_cc(:,i) = 1.e-3 * rho(:,i) / mi_g(i)  	!cm-3
	enddo
	ne_cc = 1.e-3 * rho(:,nspec+1) / me_g		!cm-3
			
	do i = 1, nspec + 1
	    T_eV(:,i) = max(  T(:,i) / qe_C, 0.  ) !avoid negative temperature
	enddo
	Te_eV = T_eV(:,nspec+1)
		
	
	do i = 1, nspec
		do j = 1, nspec
		    du = abs(  u(:,i) - u(:,j)  )
			muab(i,j) = 1.e-3 / ( mi_g(i) + mi_g(j) ) * mi_g(j) * mi_g(i)  !kg
			muab_eV(i,j) =  3.e8**2 / 1.6e-19 * muab(i,j)
			!The following L_ab is the same as the one used in LSP	
			L_ab(:,i,j) = max(  Lab_min,  23. - log( Zi(i)*Zi(j)   &
			  *   sqrt( ni_cc(:,i)*Zi(i)**2/T_eV(:,i) + ni_cc(:,j)*Zi(j)**2/T_eV(:,j) ) &
			  / ( muab_eV(i,j) * ( T_eV(:,i)/mi_eV(i) &
			      + T_eV(:,j)/mi_eV(j) + du**2./3./3.e8**2 )  )  )    &
				)
		enddo
	enddo

	!L_ab = dmax1(1.0, 23.0 - log( Z1*Z2*(m1_g + m2_g) / (m1_g*T2_eV + m2_g*T1_eV) * sqrt( n1_cc*Z1**2/T1_eV + n2_cc*Z2**2/T2_eV )  )); 

	do i = 1, nspec !logLambda for ion-electrons	
		L_ie(:,i) = max( Lab_min, log(1.5 * (Te_eV*1.6e-19*1.e7)**1.5 / sqrt(pi * ne_cc) / (Zi(i)*qe**3) ) )   !Krall&Trivelpiece
	enddo
				
	do i = 1, nspec
		do j = 1, nspec
		    du = abs(  u(:,i) - u(:,j)  )
			Tab   = ( mi_g(i) * T(:,j) + mi_g(j) * T(:,i) ) / ( mi_g(i) + mi_g(j) )
			Mab   = sqrt( 0.5 * muab(i,j) / Tab  ) * (du + 1.e-20)
			i_erf = min0( max0(1, floor(Mab / drx)  + 1), 100000) 	
			do k = 1, nz
				erf(k) = erf_table(i_erf(k)) + &
				(  erf_table(i_erf(k)+1) - erf_table(i_erf(k))  ) / drx * ( Mab(k) - (i_erf(k)-1) * drx  )
			enddo
			erf = max(0.,min(erf, 1.0)) !so that's between 0. and 1.	
			Qab = 1. / (32 * pi) * (Zi(i)*Zi(j)*qe_C**2/eps0/Tab)**2 * L_ab(:,i,j) * 3./2./Mab**3 			&
					* ( 0.5 * sqrt(pi) * erf - Mab * exp(- Mab**2) )
			
			Dab = 3. * pi / 16 * sqrt( 2 * Tab / pi / muab(i,j) ) / Qab
		
			if (check_frequency) then
				xiab(:,i,j) = mi(i) * 1.e6 * ni_cc(:,i) * 	&   
								min( 1.e6 * ni_cc(:,j) * Tab / Dab / mi(i), 1./(smax*dtm) )  !attenzione mi(i) si elide!!!!
			else
				xiab(:,i,j) = (1.e6 * ni_cc(:,i)) * (1.e6 * ni_cc(:,j)) * Tab / Dab
			endif
			
			xiab(:,i,j) = max(0.,xiab(:,i,j)) !avoid negative values too
		enddo
	enddo			
	
	do i = 1, nz
		meanL_ie(i) = sum(L_ie(i,:)) / nspec
	enddo
	taue = 3.44e5 * Te_eV**1.5 / meanL_ie !NRL page 37	- excluding density because it cancels out anyways
	ke = 3.2 * (1.6e-12 * Te_eV) * taue / me_g  !all cgs	
	!multiply by 100 to get SI units ([ke] = [L^-1*T^-1])
	ke = 1.e2 * ke
		
	!calculate ion viscosity
	do i = 1, nspec
		!NRL page 38 - excluding density because it cancels out anyways
		taui =  2.09e7 * T_eV(:,i)**1.5 / L_ab(:,i,i) * sqrt(Ai(i)) 
		mu_ion(:,i) = 0.96 * ( 1.6e-12 * T_eV(:,i) )  * taui	  !all cgs
		mu_ion(:,i) = 1.e-1 * mu_ion(:,i) ! [mu_ion] = [M*L^-1*T^-1]
	enddo
	!Redefine ion viscosity to include the factor 4./3 in front of derivative
	mu_ion = 4./3 * mu_ion
		
	!coefficients for heat transfer (NRL page 34)
	do i = 1,nspec
		do j = 1,nspec
			if (i .ne. j) then
				nu_DT(:,i,j) = 1.8e-19 * sqrt(mi_g(i)) * sqrt(mi_g(j)) * Zi(i)**2 * Zi(j)**2 &
					* ni_cc(:,j) * L_ab(:,i,j) / (mi_g(i) * T_eV(:,j) + mi_g(j) * T_eV(:,i))**1.5
			endif
		enddo
	enddo
	!heat transfer with electrons
	do i = 1,nspec
		nu_DT(:,i,nspec+1) = 1.8e-19 * sqrt(mi_g(i)) * sqrt(me_g) * Zi(i)**2 &
			* ne_cc * L_ie(:,i) / (mi_g(i) * Te_eV + me_g * T_eV(:,i))**1.5					
		nu_DT(:,nspec+1,i) = 1.8e-19 * sqrt(me_g) * sqrt(mi_g(i)) * Zi(i)**2 &
			* ni_cc(:,i) * L_ie(:,i) / (me_g * T_eV(:,i) + mi_g(i) * Te_eV)**1.5				
	enddo
	
	
	if (check_frequency) then
		do i = 1,nspec
			do j = 1,nspec
				if (i .ne. j) then
					nu_DT(:,i,j) = min(nu_DT(:,i,j),1.e-1/dtm)
				endif
			enddo
		enddo
		!heat transfer with electrons
		do i = 1,nspec
			nu_DT(:,i,nspec+1) = min(nu_DT(:,i,nspec+1),1.e-1/dtm)			
			nu_DT(:,nspec+1,i) = min(nu_DT(:,nspec+1,i),1.e-1/dtm)			
		enddo		
	endif
		
	do i = 1,nspec
		do j = 1,nspec
			k_DT(:,i,j) = 1.e6 * 3./2 * ni_cc(:,i) * nu_DT(:,i,j)   !SI units
		enddo
	enddo
	!with electrons now
	do i = 1,nspec
		k_DT(:,i,nspec+1) = 1.e6 * 3./2 * ni_cc(:,i) * nu_DT(:,i,nspec+1)   !SI units
		k_DT(:,nspec+1,i) = 1.e6 * 3./2 * ne_cc * nu_DT(:,nspec+1,i)   !SI units
	enddo

	open(unit=31, file = 'ne_cc.csv', action = 'write')
	do i = 1, nz
		write(31,'(E15.8E3)') ne_cc(i)
	enddo
	close(31)

	open(unit=31, file = 'ni_cc.csv', action = 'write')
	do i = 1, nz
		write(31,'(E15.8E3, A2, E15.8E3)') ni_cc(i,1), ",", ni_cc(i,2)
	enddo
	close(31)

    open(unit=31, file = 'Te_eV.csv', action = 'write')
    do i = 1, nz
        write(31,'(E15.8E3)') Te_eV(i)
    enddo
    close(31)

    open(unit=31, file = 'T_eV.csv', action = 'write')
    do i = 1, nz
        write(31,'(E15.8E3, A2, E15.8E3, A2, E15.8E3)') T_eV(i,1), ",", T_eV(i,2), ",", T_eV(i,3)
    enddo
    close(31)

	open(unit=31, file = 'L_ab.csv', action = 'write')
	do i = 1, nz
		write(31,'(3(E15.8E3,A2) E15.8E3) ') L_ab(i,1,1), ",", L_ab(i,1,2), ",", L_ab(i,2,1), ",", L_ab(i,2,2)
	enddo
	close(31)

	open(unit=31, file = 'L_ie.csv', action = 'write')
	do i = 1, nz
		write(31,'(E15.8E3, A2, E15.8E3)') L_ie(i,1), ",", L_ie(i,2)
	enddo
	close(31)

	open(unit=31, file = 'xiab.csv', action = 'write')
	do i = 1, nz
		write(31,'(3(E15.8E3,A2) E15.8E3) ') xiab(i,1,1), ",", xiab(i,1,2), ",", xiab(i,2,1), ",", xiab(i,2,2)
	enddo
	close(31)

	open(unit=31, file = 'taue.csv', action = 'write')
	do i = 1, nz
		write(31,'(E15.8E3)') taue(i)
	enddo
	close(31)

	open(unit=31, file = 'ke.csv', action = 'write')
	do i = 1, nz
		write(31,'(E15.8E3)') ke(i)
	enddo
	close(31)

	open(unit=31, file = 'nu_DT.csv', action = 'write')
	do i = 1, nz
		write(31,'(8(E15.8E3,A2) E15.8E3) ') nu_DT(i,1,1), ",", nu_DT(i,1,2), ",", nu_DT(i,1,3), ",",  &
				nu_DT(i,2,1), ",", nu_DT(i,2,2), ",", nu_DT(i,2,3), ",", &
				nu_DT(i,3,1), ",", nu_DT(i,3,2), ",", nu_DT(i,3,3)
	enddo
	close(31)

	open(unit=31, file = 'k_DT.csv', action = 'write')
	do i = 1, nz
		write(31,'(8(E15.8E3,A2) E15.8E3) ') k_DT(i,1,1), ",", k_DT(i,1,2), ",", k_DT(i,1,3), ",",  &
				k_DT(i,2,1), ",", k_DT(i,2,2), ",", k_DT(i,2,3), ",", &
				k_DT(i,3,1), ",", k_DT(i,3,2), ",", k_DT(i,3,3)
	enddo
	close(31)

end subroutine collision_coefficients


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine friction()
	use constants

	implicit none
	integer :: i, j
	
	R1D = 0.
	Qextra = 0.

	do i = 1, nspec
		do j = 1, nspec
			if (i .ne. j) then
				R1D(:,3*(i-1)+2) = R1D(:,3*(i-1)+2) - xiab(:,i,j) * (u(:,i) - u(:,j))
    		endif
    	enddo
    	R1D(:,3*(i-1)+3) = u(:,i) * R1D(:,3*(i-1)+2)
    enddo			
	
 	do i = 1, nspec	
		Qextra = Qextra - R1D(:,3*(i-1)+3)
	enddo

end subroutine friction

