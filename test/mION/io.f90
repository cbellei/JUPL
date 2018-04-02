!---------------------------------------------
subroutine newton_raphson(j,rr,nn,nnn,u_old,nq)

	use constants	

	implicit none
	
	integer, intent(in) :: j, nn, nnn,nq
    real*8, intent(inout) :: rr, u_old(nz)

    integer :: i, k, m
    real*8 :: mu_ion_ip12, mu_ion_im12, r_ip12, r_im12, res
    real*8 :: muold_ip12, muold_im12, m1,m2,m3,m1old,m2old,m3old,muold_ion_im12,muold_ion_ip12
    real*8 ::  u_ip12, u_im12, u_oldip12, u_oldim12
    integer :: count, INFO, ierr
    real*8 :: CC, res1, res2
    logical :: is_diag_dominant
    real*8, dimension(nz) :: Tnew, tion, TeV, taui, U1Dnew, T_before, u_before, muold, pnew
    real*8, dimension(nz) :: UU, phi_visc
	
	open(unit=31, file = 'rho.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') rho(i,j)
		enddo
	close(31)
	open(unit=31, file = 'u.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') u(i,j)
		enddo
	close(31)
	open(unit=31, file = 'u_old.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') u_old(i)
		enddo
	close(31)
	open(unit=31, file = 'U1D.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') U1D(i,j)
		enddo
	close(33)
	open(unit=31, file = 'L_ab.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') L_ab(i,j,j)
		enddo
	close(31)
	open(unit=31, file = 'T.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') T(i,j)
		enddo
	close(31)
	open(unit=31, file = 'mu_ion.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') mu_ion(i,j)
		enddo
	close(31)
	open(unit=31, file = 'U1D.csv', action = 'read', iostat = ierr)
		do i = 1,nz
			read(31,'(E30.15E3)') U1D(i,3*(j-1)+3)
		enddo
	close(33)

	!Construct diagonals (including upper and lower) and vector of knowns B
	
	CC = 5.016e22
	Tnew = T(:,j)
	p(:,j) = rho(:,j) * T(:,j) / mi(j)
	pnew = p(:,j)
	muold = mu_ion(:,j)
	
	write(*,*) "rr = ", rr
	
	
	count = 1
	do k = 1, 30
			!MOMENTUM EQUATION: Solution of Crank-Nicolson system
			!----------------------------------------------------
			BB = 0.
			AA = 0.
			UU = 0.
			phi_visc = 0.
			is_diag_dominant = .true.
			do i = nn, nz - 2 !Because AA has size (nz-2)x(nz-2)
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

				DIAG(i) = 2.*r(i)**2*rho(i+1,j) + rr * m2    &
					+ r(i)*dtm * ( 2.*mu_ion(i,j)/r(i) + (mu_ion(i+1,j)-mu_ion(i,j))/dr  )
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
					BB(i) = (2.*r(i)**2*rho(i+1,j) - rr*m2old)*u_old(i+1)    &
							 + rr*m1old*u_old(i+2) + 2.*rr*m3old*u_old(i)
				else if (i == nq-2) then !we are at bottom of array
					BB(i) = rr*m3old*u_old(i) + (2.*r(i)**2*rho(i+1,j) - rr*m2old)*u_old(i+1)				 
				else
					BB(i) = rr*m3old*u_old(i) + (2.*r(i)**2*rho(i+1,j) - rr*m2old )*u_old(i+1)   &
							+ rr*m1old*u_old(i+2)  &
							- r(i)*dtm * ( 2.*muold(i)/r(i) + (muold(i+1) - muold(i))/dr )
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


				phi_visc(i) = 3./2 * 0.5 * (mu_ion(i,j)+muold(i)) * (   &
					0.5/dr**2*( (u(i+1,j)-u(i,j))**2. + (u_old(i+1)-u_old(i))**2. ) &		
					+ ( (u(i,j)/r(i))**2. + (u_old(i)/r(i))**2. )  &
					- 0.5/3 * ( 1./dr/r(i)**2.*(r(i+1)**2*u(i+1,j) - r(i)**2*u(i,j)))**2  &
					- 0.5/3 * ( 1./dr/r(i)**2.*(r(i+1)**2*u_old(i+1) - r(i)**2*u_old(i)))**2 )
												
			enddo !do i = 1, nq - 2	

			if (k==1) then
					write(*,*) mu_ion(nn:nn+5,4)
					write(*,*) "****--------- writing to file, tm = ", tm
					call write_variable('LDIAGmio.csv',LDIAG,nz-3)
					call write_variable('UDIAGmio.csv',UDIAG,nz-3)
					call write_variable('DIAGmio.csv',DIAG,nz-2)
					call write_variable('BBmio.csv',BB,nz-2)
					call write_variable('phimio.csv',phi_visc,nz)
					call write_variable('T.csv',T(:,j),nz)
					call write_variable('U1D.csv',U1D(:,3*(j-1)+3),nz)
					call write_variable('rho.csv',rho(:,j),nz)
					call write_variable('u.csv',u(:,j),nz)
					call write_variable('u_old.csv',u_old,nz)
					call write_variable('L_ab.csv',L_ab(:,j,j),nz)
					call write_variable('radius.csv',r,nz)
					call write_variable('mu_ion.csv',mu_ion(:,j),nz)
					call write_variable('E_field.csv',Efield,nz)
					!stop
					i = nn
					write(*,*) r(i),rho(i+1,j), rr, mu_ion(i,j), dtm	
					!stop				
			endif

				if(is_diag_dominant) then !matrix is diagonally dominant -> use Thomas' algorithm
					!Thomas' algorithm
					call DGTSV(nnn,1,LDIAG(nn:nz-3),DIAG(nn:nz-2),UDIAG(nn:nz-3),BB(nn:nz-2),nnn,INFO)	
				else
					!LU decomposition
					call DGESV(nnn,1,AA(nn:nz-2,nn:nz-2),nnn,IPIV(nn:nz-2),BB(nn:nz-2),nnn,INFO)
				endif	

			u_before = u(:,j)
			!IF COMMENTED, DO NOT SOLVE MOMENTUM EQUATION
			u(nn+1:nq-1,j) = BB(nn:nz-2)
		
			T_before = Tnew !first,record old value of Tnew

			UU(nn+1:nq-1) = U1D(nn+1:nq-1,3*(j-1)+3) + dtm * phi_visc(nn+1:nq-1)		
			Tnew(nn+1:nq-1) = (g-1) * mi(j) * ( UU(nn+1:nq-1)/rho(nn+1:nq-1,j) - 0.5*u(nn+1:nq-1,j)**2 )
			Tnew = max( Tnew, Tmin ) !avoid negative temperatures
		
			TeV = Tnew / qe_C
			pnew = rho(:,j) * Tnew / mi(j)
			!NRL page 38 - excluding density because it cancels out anyways
			taui =  2.09e7 * TeV**1.5 / L_ab(:,j,j) * sqrt(Ai(j)) 
			mu_ion(:,j) = 0.96 * ( 1.6e-12 * TeV )  * taui	  !all cgs
			mu_ion(:,j) = 1.e-1 * mu_ion(:,j) ! [mu_ion] = [M*L^-1*T^-1]
			!Redefine ion viscosity to include the factor 4./3 in front of derivative
			mu_ion(:,j) = 4./3 * mu_ion(:,j)

			if (k==1) then
				write(*,*) "aqui"
				write(*,*) "U1D = ", U1D(nn+1,3*(j-1)+3), "UU=", UU(nn+1)
				write(*,*) "phi_visc =", phi_visc(nn+1), "dtm =", dtm
				write(*,*) "Tnew = ", Tnew(nn+1), "mi =", mi(j), "u = ", u(nn+1,j)
				call write_variable('BBsol.csv',BB,nz-2)
				call write_variable('U1Dsol.csv',U1D(:,3*(j-1)+3),nz)
				call write_variable('UU.csv',UU,nz)
				call write_variable('musol.csv',mu_ion(:,j),nz)
				call write_variable('TeV.csv',TeV,nz)
				call write_variable('taui.csv',taui,nz)
				call write_variable('Tnew.csv',Tnew,nz)
				call write_variable('T_before.csv',T_before,nz)
				!stop
			endif

			
			if (count > 1) then
				res1 = maxval(abs(u(nn+1:nq-1,j)-u_before(nn+1:nq-1))/abs(u(nn+1:nq-1,j)))
				res2 = maxval(abs(Tnew(nn+1:nq-1)-T_before(nn+1:nq-1))/Tnew(nn+1:nq-1))
				write(*,*) "res1 = ", res1,"res2 = ", res2
				stop
			endif
			
			count = count + 1
			
			if (k==2) then
				stop
			endif
			
			
			if (count>10000) then
				write(*,*) "--------------------------------------------------------------"
				write(*,*) "Newton-Raphson not converging - Aborting viscosity calculation"
				write(*,*) "species  = ", j
				write(*,*) "count = ", count
				write(*,*) "tm = ", tm
				write(*,*) "res = ", res
				write(*,*) "--------------------------------------------------------------"	
				stop				
			endif
						
		enddo !while loop
		
end subroutine newton_raphson


!---------------------------------------------
subroutine write_variable(filename,var,len)
	use constants
	
	implicit none
	integer :: k, ierr
	integer, intent(in) :: len
	real*8, intent(in) :: var(len)
	character(len=*), intent(in) :: filename
	
	write(*,*) "filename = ", filename
	
	open(unit = 9999, file = trim(filename), status = 'replace', iostat=ierr)
		do k = 1, len				
			write(9999,1005) var(k)
		enddo				
	close(9999)
	
	1005 format(E30.15E3)

end subroutine write_variable


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine write_data()
	use constants
		
	implicit none
	real*8, dimension(nz) :: var
	integer :: i,k,c1,c2,c3
		
	c1 = 0
	c2 = 0
	c3 = 0
	do i = 1, (nspec+1)*3+2
		if (i==1) then
			var = r
		else if (i==(nspec+1)*3+2) then
			var = efield
		else if (i>1 .and. i<=(nspec+1)+1) then
			c1 = c1 + 1
			var = u(:,c1)		
		else if (i>(nspec+1)+1 .and. i<=2*(nspec+1)+1) then
			c2 = c2 + 1
			if (c2<=nspec) then
				var = 1.e-6 * rho(:,c2)	/ mi(c2)	
			else !electrons
				var = 1.e-6 * rho(:,c2)	/ me			
			endif
		else if (i>2*(nspec+1)+1 .and. i<=3*(nspec+1)+1) then
			c3 = c3 + 1
			var = T(:,c3)/qe_C
		endif
		do k = 1, nz				
			write(unt+i,1005) var(k)
		enddo				
	enddo

	1005 format(E34.22E3)

end subroutine write_data


subroutine open_files()
	use constants

	implicit none
	integer :: i, c1, c2, c3, ierr
	character(len=20) :: filename, x1
	character(len=8) :: fmt ! format descriptor

	fmt = '(I1)' ! an integer of width 5 with zeros at the left

	c1 = 0
	c2 = 0
	c3 = 0

	do i = 1, (nspec+1)*3+2
		if (i==1) then
			filename = 'r.csv'
		else if (i==(nspec+1)*3+2) then
			filename = 'efield.csv'
		else if (i>1 .and. i<=(nspec+1)+1) then
			c1 = c1 +1
			write (x1,fmt) c1
			filename = trim('vel'//trim(x1)//'.csv')
		else if (i>(nspec+1)+1 .and. i<=2*(nspec+1)+1) then
			c2 = c2 +1
			write (x1,fmt) c2
			filename = trim('den'//trim(x1)//'.csv')
		else if (i>2*(nspec+1)+1 .and. i<=3*(nspec+1)+1) then
			c3 = c3 +1
			write (x1,fmt) c3
			filename = trim('temp'//trim(x1)//'.csv')
		endif
		open(unit = unt+i, file = filename, status = 'replace', iostat=ierr)
		if(ierr/=0) then
			write(*,*) 'problems in opening the file', filename
			stop
		endif
	enddo

end subroutine open_files

subroutine close_files()
	use constants

	implicit none
	integer :: i

	do i = 1, (nspec+1)*3+2
		close(unt + i)
	enddo

end subroutine close_files

!--------------------------------------------

subroutine write_all_data(U1Dv, U1D_pv, U1D_cv, F1Dv, G1Dv, C1Dv)

    use constants

    integer :: i
	real*8, intent(in), dimension(nz,neqi+1) :: U1Dv, U1D_pv, U1D_cv, F1Dv, G1Dv
    real*8, intent(in), dimension(nz,neqi) :: C1Dv


	open(unit=31, file = 'U1D.csv', action = 'write')
	do i = 1, nz
		write(31,'(6(E34.22E3,A2) E34.22E3)') &
                U1Dv(i,1), ",", U1Dv(i,2), ",", U1Dv(i,3), ",", U1Dv(i,4), ",", U1Dv(i,5), &
                ",", U1Dv(i,6), ",", U1Dv(i,7)
	enddo
	close(31)

	open(unit=31, file = 'U1D_p.csv', action = 'write')
	do i = 1, nz
		write(31,'(6(E34.22E3,A2) E34.22E3)') &
				U1D_pv(i,1), ",", U1D_pv(i,2), ",", U1D_pv(i,3), ",", U1D_pv(i,4), ",", U1D_pv(i,5), &
				",", U1D_pv(i,6), ",", U1D_pv(i,7)
	enddo
	close(31)

	open(unit=31, file = 'U1D_c.csv', action = 'write')
	do i = 1, nz
		write(31,'(6(E34.22E3,A2) E34.22E3)') &
				U1D_cv(i,1), ",", U1D_cv(i,2), ",", U1D_cv(i,3), ",", U1D_cv(i,4), ",", U1D_cv(i,5), &
				",", U1D_cv(i,6), ",", U1D_cv(i,7)
	enddo
	close(31)

	open(unit=31, file = 'F1D.csv', action = 'write')
	do i = 1, nz
		write(31,'(6(E34.22E3,A2) E34.22E3)') &
				F1Dv(i,1), ",", F1Dv(i,2), ",", F1Dv(i,3), ",", F1Dv(i,4), ",", F1Dv(i,5), &
				",", F1Dv(i,6), ",", F1Dv(i,7)
	enddo
	close(31)

	open(unit=31, file = 'G1D.csv', action = 'write')
	do i = 1, nz
		write(31,'(6(E34.22E3,A2) E34.22E3)') &
				G1Dv(i,1), ",", G1Dv(i,2), ",", G1Dv(i,3), ",", G1Dv(i,4), ",", G1Dv(i,5), &
				",", G1Dv(i,6), ",", G1Dv(i,7)
	enddo
	close(31)

	open(unit=31, file = 'C1D.csv', action = 'write')
	do i = 1, nz
		write(31,'(6(E34.22E3,A2))') &
				C1Dv(i,1), ",", C1Dv(i,2), ",", C1Dv(i,3), ",", C1Dv(i,4), ",", C1Dv(i,5), &
				",", C1Dv(i,6)
	enddo
	close(31)

	open(unit=31, file = 'T.csv', action = 'write')
	do i = 1, nz
		write(31,'(2(E34.22E3,A2) E34.22E3)') &
				T(i,1), ",", T(i,2), ",", T(i,3)
	enddo
	close(31)


end subroutine write_all_data

