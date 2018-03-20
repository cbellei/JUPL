!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!    SUBROUTINES FOR BOUNDARY CONDITIONS    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_BC()
	use constants

	implicit none
	integer :: i,kk,nn,nn1
	
	
	do kk = 1,2
		if(kk==1)then
			nn = 1
			nn1 = 2
		else
			nn = nz
			nn1 = nz-1
		endif
		if (boundary(kk)=='reflection') then
			if (geom=="slab") then
				!boundary condition for reflection at wall side
				do i = 1, nspec
					U1D(nn,3*(i-1)+1) = U1D(nn1,3*(i-1)+1) - lm * F1D(nn1,3*(i-1)+1) !continuity
					U1D(nn,3*(i-1)+2) = 0. !momentum
					U1D(nn,3*(i-1)+3) = U1D(nn1,3*(i-1)+3) - lm * F1D(nn1,3*(i-1)+3) !energy
				enddo
				!now electrons
				U1D(nn,neqi+1) =  U1D(nn1,neqi+1) - lm * F1D(nn1,neqi+1)
			else if(geom=="spherical") then!spherical geometry
				!boundary condition for reflection at wall side
				do i = 1, nspec
					U1D(nn,3*(i-1)+1) = U1D(nn1,3*(i-1)+1) - lm * F1D(nn1,3*(i-1)+1) 	&
								- 2. * dtm_ / r(nz) * phi * F1D(nz-1,3*(i-1)+1)
					U1D(nn,3*(i-1)+2) = 0. !momentum
					U1D(nn,3*(i-1)+3) = U1D(nn1,3*(i-1)+3) - lm * F1D(nn1,3*(i-1)+3)  &
								- 2. * dtm_ / r(nz) * phi * F1D(nn1,3*(i-1)+3) 
				enddo
				!now electrons
				U1D(nz,neqi+1) =  U1D(nz-1,neqi+1) 
			endif		
		else if (boundary(kk)=='open') then
			!boundary condition for reflection at wall side
			do i = 1, nspec
				U1D(nn,3*(i-1)+1) = rho(nn,i) 
				U1D(nn,3*(i-1)+2) = rho(nn,i) * u(nn,i) !momentum
				U1D(nn,3*(i-1)+3) = p(nn,i) / (g-1) + 0.5 * rho(nn,i) * u(nn,i)**2
			enddo
			!now electrons
			U1D(nn,neqi+1) =  U1D(nn,neqi+1) 
		else
			write(*,*) "wrong boundary condition at boundary #", kk
		endif
		
		!update variables
		do i = 1, nspec
			rho(nn,i) = U1D(nn,3*(i-1)+1)
			u(nn,i) = U1D(nn,3*(i-1)+2) / rho(nn,i)
			T(nn,i) = (g-1) * mi(i) * ( U1D(nn,3*(i-1)+3)/rho(nn,i) - 0.5*u(nn,i)**2 )
			p(nn,i) = rho(nn,i) * T(nn,i) / mi(i)
		enddo
		
	enddo

	!BC for predictor and corrector
	U1D_p(1,:) = U1D(1,:)
	U1D_c(1,:) = U1D(1,:)
	U1D_p(nz,:) = U1D(nz,:)
	U1D_c(nz,:) = U1D(nz,:)
	
end subroutine apply_BC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBROUTINE ART_VISCOSITY    !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine artviscosity()

	use constants
	implicit none
	integer :: m,j,k
	real*8 :: a, b, y1, y2, x1, x2

	if (geom=="slab") then
		do k = 1, nz
		   eps_visc(k,1:neqi+1) = eps_visc_max
		enddo
	else !spherical
		y1 = log10(eps_visc_max) 
		y2 = log10(eps_visc_max / eps_compress)
		x1 = log10(r(nz))
		x2 = log10(r(1))
		b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
		a  = ( y1 - b ) / x1 
		do k = 1,nz
			eps_visc(k,1:neqi+1) = r(k)**a * 10**b
		enddo
		write(*,*) "eps_visc(1,1) = ", eps_visc(1,1), "eps_visc(nz,1) = ", eps_visc(nz,1)
	endif


end subroutine artviscosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBROUTINE  INITSPACEGRID    !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initspacegrid()
	use constants

	implicit none
	integer :: m,j,k
	real*8 :: a, b, y1, y2, x1, x2
	
	if(geom=="spherical") then	
		r(1) = L + rmin
		r0(:,nregions+1) = r0(:,nregions+1) + rmin !correct also r0 array
	else
		r(1) = 0.
	endif
	
	do k = 2, nz
		r(k) = r(k-1) + dr !note that dr is negative
	enddo
	
	if(geom=="spherical") then	
		r(nz) = rmin
	endif
		
end subroutine initspacegrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBROUTINE  INITVARIABLES    !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initvariables()

	use constants

	implicit none
	integer :: m,j,k

	do j = 1, nspec
		do m = 1, nregions
			do k = 1, nz
				if ( r(k) >= r0(j,m) .and. r(k) <= r0(j,m+1)) then
					T0(k,j) = temp0(j,m)
					N0(k,j) = den0(j,m)
					V0(k,j) = vel0(j,m)
				endif
			enddo
		enddo
	enddo
	
	rho = 0.
	u   = 0.
	p   = 0.
	T   = 0.
	
	!ions				
    do k = 1, nz
    	do j = 1, nspec
			rho(k,j) = N0(k,j) * mi(j)
			u(k,j) = V0(k,j)
			p(k,j) = rho(k,j) / mi(j) * T0(k,j)
			T(k,j) = T0(k,j)		
		enddo
    enddo

	u(1,:) = 0.0

end subroutine initvariables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!    SUBROUTINE  DO_SMOOTHING    !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine do_smoothing()

	use constants

	implicit none
	integer :: m, j, k, jj, kgrad(2)=0, m1, m2
	real*8 :: a, b, y1, y2, x1, x2
	real*8 :: grad, mingrad, maxgrad
	real*8 :: nA(2)=0., nB(2)=0., uA(2)=0., uB(2)=0., TA(2)=0., TB(2)=0.
	integer :: sgn1(2) = (/0,1/), sgn2(2) = (/1,0/), ok(2) = 0
	real*8, parameter :: epsilon = 1.e-10
	real*8 :: vel(nz)
	integer, dimension(nedges) :: mm
		
            
	if (smoothing) then
			
		!now smooth out interface between nz0 and nz0+nsmooth
		!----------------------------------------------------
		if (smooth_vars(1)==1) then
			write(*,*) "density smoothing function"
			write(*,*) "--------------------------"
			!for each ion species, find position of max and min gradients, then smooth out density	
			do j = 1, nspec
				maxgrad = 0.
				mingrad = 1.e20
				do k = 1, nz - 1
					grad = rho(k,j) / rho(k+1,j)
					if (grad > maxgrad) then
						maxgrad = grad
						kgrad(1) = k
						nA(1) = rho(k,j) / mi(j)
						nB(1) = rho(k+1,j) / mi(j)
					endif
					if (grad < mingrad) then
						mingrad = grad
						kgrad(2) = k
						nA(2) = rho(k,j) / mi(j)
						nB(2) = rho(k+1,j) / mi(j)
					endif			
				enddo
				
				do k = 1,nedges
					if(smooth_type(k,j)=="max") then
						mm(k)=1
					elseif(smooth_type(k,j)=="min") then
						mm(k)=2
					endif
				enddo
								
				do jj = 1,nedges
					m = mm(jj)
					do k = kgrad(m)-sgn1(m)*nsmooth, kgrad(m) + sgn2(m)*nsmooth
							!now smooth out interface
							y1 = log10( nA(m)  ) 
							y2 = log10( nB(m)  )
							x1 = log10( r(kgrad(m) - sgn1(m)*nsmooth) ) !note r is in micron
							x2 = log10( r(kgrad(m) + sgn2(m)*nsmooth) ) !note r is in micron
							b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
							a  = ( y1 - b ) / x1 

							rho(k,j) = (  10**( a*log10(r(k)) + b )   ) * mi(j)
					enddo	
			
				enddo
			enddo
		endif
		
		if (smooth_vars(2)==1) then 
			write(*,*) "velocity smoothing function"
			write(*,*) "--------------------------"
			!for each ion species, find position of max and min gradients, then smooth out density
			ok = 0	
			do j = 1, nspec
				maxgrad = 0.
				mingrad = 1.e20
				vel = max(epsilon,abs(u(:,j)))
				do k = 1, nz - 1
					grad = vel(k) / vel(k+1)
					if (grad > maxgrad) then
						maxgrad = grad
						kgrad(1) = k
						uA(1) = vel(k) 
						uB(1) = vel(k+1)
						if ( k+tsmooth < nz .and. k-tsmooth>0 ) then
							ok(1) = 1
						endif					
					endif
					if (grad < mingrad) then
						mingrad = grad
						kgrad(2) = k
						uA(2) = vel(k)   
						uB(2) = vel(k+1)
						if ( k+tsmooth < nz .and. k-tsmooth>0 ) then
							ok(2) = 1
						endif					
					endif			
				enddo
				
				do k = 1,nedges
					if(smooth_type(k,j)=="max") then
						mm(k)=1
					elseif(smooth_type(k,j)=="min") then
						mm(k)=2
					endif
				enddo
								
				do jj = 1,nedges
					m = mm(jj)
					if(ok(m).ne.1) cycle !because it is a bad point...
					do k = kgrad(m)-sgn1(m)*vsmooth, kgrad(m) + sgn2(m)*vsmooth
							!now smooth out interface
							y1 = log10( uA(m) ) 
							y2 = log10( uB(m) )
							x1 = log10( r(kgrad(m) - sgn1(m)*vsmooth) ) !note r is in micron
							x2 = log10( r(kgrad(m) + sgn2(m)*vsmooth) ) !note r is in micron
							b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
							a  = ( y1 - b ) / x1 

							u(k,j) = - (  10**( a*log10(r(k)) + b )   ) 
					enddo	
			
				enddo
			enddo
		endif
		
		if(smooth_vars(3)==1) then		
		
			write(*,*) "temperature smoothing function"
			write(*,*) "------------------------------"
			!for each ion species, find position of max and min gradients, then smooth out density
			ok = 0
			do j = 1, nspec
				maxgrad = 0.
				mingrad = 1.e20
				do k = 1, nz - 1
					grad = T(k,j) / T(k+1,j)
					if (grad > maxgrad) then
						maxgrad = grad
						kgrad(1) = k
						TA(1) = T(k,j) 
						TB(1) = T(k+1,j)
						if ( k+tsmooth < nz .and. k-tsmooth>0 ) then
							ok(1) = 1
						endif					
					endif
					if (grad < mingrad) then
						mingrad = grad
						kgrad(2) = k
						TA(2) = T(k,j)   
						TB(2) = T(k+1,j)
						if ( k+tsmooth < nz .and. k-tsmooth>0 ) then
							ok(2) = 1
						endif					
					endif			
				enddo
			
				do k = 1,nedges
					if(smooth_type(k,j)=="max") then
						mm(k)=1
					elseif(smooth_type(k,j)=="min") then
						mm(k)=2
					endif
				enddo
								
				do jj = 1,nedges
					m = mm(jj)
					if(ok(m).ne.1) cycle !because it is a bad point...
					do k = kgrad(m)-sgn1(m)*tsmooth, kgrad(m) + sgn2(m)*tsmooth
							!now smooth out interface
							y1 = log10( TA(m) ) 
							y2 = log10( TB(m) )
							x1 = log10( r(kgrad(m) - sgn1(m)*tsmooth) ) !note r is in micron
							x2 = log10( r(kgrad(m) + sgn2(m)*tsmooth) ) !note r is in micron
							b  = ( y1*x2 - y2*x1 ) / ( x2 - x1 )
							a  = ( y1 - b ) / x1 

							T(k,j) = (  10**( a*log10(r(k)) + b )   ) 
					enddo	
			
				enddo			
			enddo
		endif
		
		!now update pressure after all this smoothing...
		do j = 1, nspec
			p(:,j) = rho(:,j)/mi(j) * T(:,j)
		enddo				
		
	endif		
end subroutine do_smoothing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------- SUBROUTINE READ_FILES -----------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_files()

	use constants
		
	implicit none
	
	real*8, dimension(nz) :: var
	real*8 :: Tmean
	integer :: i,k,c1,c2,c3,ierr
	character(len=20) :: filename, x1
	character(len=8) :: fmt ! format descriptor

	fmt = '(I1)' ! an integer of width 5 with zeros at the left	
	c1 = 0
	c2 = 0
	c3 = 0
	
	write(*,*) "--------------------------------"
	do i = 1,(nspec+1)*3+2
		if (i==1) then
			filename = 'r.dat'
			write(*,*) "reading ", filename
			write(*,*)
			open(unit=unt, file = filename, action = 'read', iostat = ierr)
				if(ierr/=0) then
					write(*,*) 'problems in opening file r.dat'
					stop
				else
					do k = 1, nz
						read(unt,'(E20.10E3)') r(k)
					enddo
				endif
			close(unt)
		else if (i==(nspec+1)*3+2) then
			filename = 'efield.dat'
			write(*,*) "reading ", filename
			open(unit=unt, file = filename, action = 'read', iostat = ierr)
				if(ierr/=0) then
					write(*,*) '  warning - efield.dat does not exist'
				else
					do k = 1, nz
						read(unt,'(E20.10E3)') Efield(k)
					enddo
				endif
			close(unt)
		else if (i>1 .and. i<=(nspec+1)+1) then
			c1 = c1 +1
			write (x1,fmt) c1
			filename = trim('vel'//trim(x1)//'.dat')
			write(*,*) "reading ", filename
			open(unit=unt, file = filename, action = 'read', iostat = ierr)
				if(ierr/=0) then
					write(*,*) 'problems in opening file', filename
					stop
				else
					do k = 1, nz
						read(unt,'(E20.10E3)') u(k,c1)  !m/s
					enddo
				endif
			close(unt)
			if (i==nspec+2) write(*,*)
		else if (i>(nspec+1)+1 .and. i<=2*(nspec+1)+1) then
			c2 = c2 +1
			write (x1,fmt) c2
			filename = trim('rho'//trim(x1)//'.dat') 
			write(*,*) "reading ", filename
			open(unit=unt, file = filename, action = 'read', iostat = ierr)
				if(ierr/=0) then
					write(*,*) 'problems in opening file', filename
					stop
				else
					do k = 1, nz
						read(unt,'(E20.10E3)') rho(k,c2)  !kg/cm3
					enddo
				endif
			close(unt)
			if (i==2*(nspec+1)+1) write(*,*)
		else if (i>2*(nspec+1)+1 .and. i<=3*(nspec+1)+1) then
			c3 = c3 +1
			write (x1,fmt) c3
			filename = trim('temp'//trim(x1)//'.dat')			
			write(*,*) "reading ", filename
			open(unit=unt, file = filename, action = 'read', iostat = ierr)
				if(ierr/=0) then
					write(*,*) 'problems in opening file', filename
					stop
				else
					do k = 1, nz
						read(unt,'(E20.10E3)') T(k,c3)  !Joule
					enddo
				endif
			close(unt)
			if (i==3*(nspec+1)+1) write(*,*)
		endif
	enddo
	write(*,*) "--------------------------------"
	
	!define position of boundary
	L = r(1)
	
	!make temperature equal at border
	Tmean = sum(T(1,:))/(nspec+1)
	T(1,:) = Tmean

end subroutine read_files


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------- SUBROUTINE DEFINE FLUXES ----------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine define_fluxes()

	use constants
	
	implicit none
	integer :: k,j

	U1D = 0.
	F1D = 0.
	R1D = 0.
	
	write(*,*)
	write(*,*) " WARNING - forcing quasi-neutrality and zero-current conditions"
	write(*,*)
	
	rho(:,nspec+1) = 0.
	u(:,nspec+1) = 0.
		
	!electrons
    do k = 1, nz
    	do j = 1, nspec
			rho(k,nspec+1) = rho(k,nspec+1) + me * Zi(j) * rho(k,j) / mi(j) !quasi-neutrality
			u(k,nspec+1) = u(k,nspec+1) + Zi(j) * rho(k,j) * u(k,j) / mi(j)   !zero-current condition  
    	enddo
    	if(.not.(restart)) then !take species 1 as reference
			p(k,nspec+1) = rho(k,nspec+1) / me * T(k,1) !electron pressure
			T(k,nspec+1) = T(k,1) !electron temperature	
		else !we have read the electron temperature from file
			p(k,nspec+1) = rho(k,nspec+1) / me * T(k,nspec+1) !electron pressure			
		endif	
    enddo
        
!   correct for electron velocity
    u(:,nspec+1) = u(:,nspec+1) / ( rho(:,nspec+1) / me )
!	apply BC at piston side
!    do j = 1, nspec+1
!        u(1,j) = 0.!V0(1,j)!0.
!    enddo
	

!	finally, calculate U1D and F1D for all species
	do j = 1, nspec
		U1D(:,3*(j-1)+1) = rho(:,j)
		U1D(:,3*(j-1)+2) = rho(:,j) * u(:,j) 
		U1D(:,3*(j-1)+3) = p(:,j) / (g-1) + 0.5 * rho(:,j) * u(:,j)**2
		F1D(:,3*(j-1)+1) = rho(:,j) * u(:,j)
		F1D(:,3*(j-1)+2) = rho(:,j) * u(:,j)**2 + p(:,j) 
		F1D(:,3*(j-1)+3) = u(:,j) * ( g / (g-1) * p(:,j) + 0.5 * rho(:,j) * u(:,j)**2 )
	enddo
	!electrons
	U1D(:,neqi+1) = p(:,nspec+1) / (g-1) + 0.5 * rho(:,nspec+1) * u(:,nspec+1)**2
	F1D(:,neqi+1) = u(:,nspec+1) * ( g / (g-1) * p(:,nspec+1) + 0.5 * rho(:,nspec+1) * u(:,nspec+1)**2 )	

end subroutine define_fluxes




