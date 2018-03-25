!July 20th, 2015 - adding a third ion species
!Sep 12th, 2014 - shock driven by a piston in a neutralized plasma of e- and i+
!Sep 4th, 2014 - making MacCormack work
!August 17th, 2014 - trying with MacCormack
!include table for psi_integral - August 8th, 2014
!adding two species - August 1st, 2014
!adding diffusion equation - July 27, 2014
!piston01 - July 21, 2014
!Sod02 - July 20, 2014
!Riemann problem in Eulerian coordinates
!July 19, 2014 - Claudio Bellei

!-------------------------------------------------------------
!  ALGORITHM OF SOLUTION IS MacCormack (Laney page 361)  !
!-------------------------------------------------------------

program shock

	use constants

	implicit none

	real*8 :: dxx = 1.e-4, dummy
	real*8, dimension(100001) :: erf_table
	integer :: ierr, i, j, k
	integer :: nq !for quiet start
	integer*4 :: today(3), now(3)
	logical :: time1 = .true., time2 = .true.

	call idate(today)   ! today(1)=day, (2)=month, (3)=year
	call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
	write ( *, 999 )  today(2), today(1), today(3), now

	call initialize()

	!check geometry
	if (geom/="slab" .and. geom/="spherical") then
		write(*,*)
		write(*,*) "Invalid geometry"
		write(*,*)
		stop
	endif

	!-------------------------------
	!-------------------------------		
	write(*,*) 'dtm = ', dtm

	call initspacegrid()

	call artviscosity()

	if(.not.(restart)) then
		call initvariables()
	endif

	call do_smoothing()

	call define_fluxes()

	!initialize predictor and corrector
	U1D_p = U1D
	U1D_c = U1D
	Efield = 0.
	G1D = 0.
	Q_DT = 0.
	q_diff = 0.

	open(unit=unt, file = './erf_integral.dat', action = 'read', iostat = ierr)
	if(ierr/=0) then
		write(*,*) 'problems in opening file erf_integral.dat'
		stop
	else
		do i = 1, 100001
			read(unt,'(2E20.10E3)') dummy, erf_table(i)
		enddo
	endif
	close(unt)

	open(unit=unt, file = 'parameters.csv', status = 'replace', iostat = ierr)
	if(ierr/=0) then
		write(*,*) 'problems in opening the file'
		stop

	else
		write(unt,'(E15.6E3,A2,I5,A2,E15.6E3,A2,I13,A2,E15.6E3,A2,E15.6E3)') &
				dtm, ', ', nz, ', ', dt_print, ', ', maxind, ', ', tm_quiet, ', ', L
	endif
	close(unt)

	call open_files()
	call write_data()

	write(*,*) "C'est parti!"
	do j = 1, maxind

		if (tm>maxTime) then
			write(*,*) "---Simulation ended---"
			stop
		endif
		if (mod(tm,dt_print) < dtm) then
			ifprint = .true.
			write(*,*) "--writing files,   t = ", tm
		endif


		if ( mod(100*j, maxind)==0 ) then
			write(*,'(I3,A26)') int( float(j*100) / maxind ), "  % of the simulation done"
			write(*,'(A15,E10.5E2,I10)') "      t, j = ", tm, j
		endif

		tm = tm + dtm

		if (time1 .and. tm < tm_quiet) then
			if(geom=="slab") then
				nq = nquiet
			else if (geom=="spherical") then
				nq = nz-nquiet
			endif
			write(*,*) "nz = ", nz, "nquiet = ", nquiet, "tm_quiet=", tm_quiet
			time1 = .false.
		elseif (time2 .and. tm >= tm_quiet) then
			write(*,*) "end of quiet time"
			nq = nz
			dtm = dtm / dt_multi
			lm = dtm / dr
			time2 = .false.
		endif

		do i = 1, nz
			do k = 1, nspec
				if (  u(i,k) .ne. u(i,k)  ) then
					write(*,*) 'here is the flaw, i,j=', i,j, "tm = ", tm
					write(*,*) "species = ", k
					write(*,*) "r = ", r(i)
					write(*,*) u(i-3:i+3,1)
					stop
				endif
			enddo
		enddo

		call apply_BC()

		!----- predictor -----------
		!---------------------------
		call predictor(nq,j)
		call update_variables(U1D_p)
		call calculate_collisions(dxx, erf_table)
		call source_terms(nq)

		!----- corrector -----------
		!---------------------------
		call corrector(nq,j)

		!----- FINAL STEP -------
		U1D(2:nq-1,1:neqi+1) = 0.5 * ( U1D(2:nq-1,1:neqi+1) + U1D_c(2:nq-1,1:neqi+1) ) 	&
				+ eps_visc(1:nq-2,1:neqi+1) * ( U1D(3:nq,1:neqi+1) - 2*U1D(2:nq-1,1:neqi+1) + U1D(1:nq-2,1:neqi+1)  )

		call update_variables(U1D)
        call calculate_collisions( dxx,erf_table)
        call source_terms(nq)
        !--------------------------------
        call write_all_data(U1D, U1D_p, U1D_c, F1D, G1D, C1D)
        stop
        !--------------------------------


		if (ion_viscosity) then !second step in splitting operator algorithm

			call calculate_viscous_term(nq,j)

			call update_variables(U1D)
			call calculate_collisions(dxx,erf_table)
			call source_terms(nq)

			call apply_BC()
			call predictor(nq,j)
			call update_variables(U1D_p)
			call calculate_collisions( dxx, erf_table )
			call source_terms(nq)

			call corrector(nq,j)

			U1D(2:nq-1,1:neqi+1) = 0.5 * ( U1D(2:nq-1,1:neqi+1) + U1D_c(2:nq-1,1:neqi+1) ) 	&
					+ eps_visc(1:nq-2,1:neqi+1) * ( U1D(3:nq,1:neqi+1) - 2*U1D(2:nq-1,1:neqi+1) + U1D(1:nq-2,1:neqi+1)  )

			call update_variables(U1D)
			call calculate_collisions( dxx,erf_table)
			call source_terms(nq)

		endif

		if (ifprint) then
			call write_data()
			ifprint = .false.
		endif

	enddo

	write(*,*) 'Final time = ', tm

	call close_files()

	999  format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', i2.2, ':', i2.2, ':', i2.2 )

end program  shock