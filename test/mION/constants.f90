module constants

	real*8, parameter :: mp_g = 1.*1836.*9.11e-28; !proton mass in grams
	real*8, parameter :: me = 9.11e-31; !kg
    real*8, parameter :: me_g = 9.11e-28
	real*8, parameter :: qe_C = 1.6e-19;
	real*8, parameter :: qe = 4.8032e-10; !statcoulomb
	real*8, parameter :: g = 5./3. !gamma
	real*8, parameter :: eps0 = 8.854187e-12
	integer, parameter :: unt = 50 !unit for opening files
	
	integer :: nspec, neqi, nregions, nsubcycle
	real*8 :: v_piston, v_shell, n_left, T_left, n2n1, T2T1, n_shocked, T_shocked, T_ablation, Tmin
	real*8 :: n_right, T_right
	real*8 :: L, Lab_min
	integer :: maxind, nz, dnz, dt_multi, nquiet, nsmooth, vsmooth, tsmooth, nz0, nz00
	integer :: viscous_species
	real*8 :: dr, eps_visc_max, eps_compress, eps_n, phi, tm_quiet
	logical :: heattransfer_switch, Efield_switch, electron_switch, friction_switch
	logical :: smoothing, bremsstrahlung, ion_viscosity, restart, electron_heat_conduction
	logical :: check_frequency, ifprint
	real*8 :: flimit, ni_frac, nmin, rmin, dt_print
	real*8 :: lm, dtm, dtm_, CFL, cs, coeff_CFL, Tmax, smax, maxTime
	character*9 :: geom
	
	real*8, dimension(:), allocatable :: mi, mi_g, qi, Ai, nleft, nmiddle, nright
	real*8, dimension(:), allocatable :: Zi, r_shell, dr_shell
    real*8, dimension(:,:), allocatable :: p, u, rho, T, ni, r0, den0, temp0, vel0
    real*8, dimension(:,:), allocatable :: U1D, F1D, R1D, U1D_p, U1D_c, G1D, Q_DT, eps_visc, C1D
    real*8, dimension(:,:), allocatable :: T0, N0, V0, mu_ion, AA, AA2
    real*8, dimension(:), allocatable   :: LDIAG, DIAG, UDIAG, BB, IPIV
    real*8, dimension(:), allocatable   :: LDIAG2, DIAG2, UDIAG2, BB2, IPIV2
    real*8, dimension(:), allocatable :: r, xi, xi12, nu0, du, ke, vth_e, q_FS, q_diff, Te, Efield
    real*8, dimension(:), allocatable :: ne, Qextra 	
    real*8, dimension(:,:,:), allocatable :: k_DT, xiab, L_ab	
    real*8 :: tm = 0.0, xmax_visc
    integer, dimension(3) :: smooth_vars
	character*9, dimension(:,:), allocatable :: smooth_type
	!character(len=*), dimension(*) :: boundary 
	character(len=12), dimension(2) :: boundary
    integer :: nedges
    
    save
	
	contains
	
	subroutine initialize()
		implicit none
		integer :: i
		real :: dr_smooth

		geom = "slab" 
		boundary = [character(len=12) :: "reflection", "open"]
		flimit = 0.01 !flux limit


		heattransfer_switch = .true. !if .false., no thermal equilibration
		efield_switch = .true.
		electron_switch = .true.   ! if .false., electrons are uncoupled to ions (no E-field, no thermalization with electrons)	
		electron_heat_conduction = .true.
		friction_switch = .true.
		bremsstrahlung = .true.
		check_frequency = .true. !check if put a max limit on frequency for friction
		smax = 5. ! max time step for friction (only used when check_frequency==true)
		ion_viscosity = .false.
		restart = .false.
		viscous_species = 4
		nsubcycle = 1
		smoothing = .true.
 		dr_smooth = 2.5e-4 !smoothing range

		coeff_CFL = 1.e-4
		dt_print = 5.e-10

		nspec = 2
		nz = 2000 !number of zones
		neqi = 3*nspec !total # of equations for the ions
		dnz = 40
		Tmax = 500.e3 * qe_C !in Joules
		Lab_min = 2.

		nregions = 2
		nedges = 1
		
		L = 2.5e-2
		maxTime = 2.e-7
		
		call allocate_arrays()

		smooth_type(1:nedges,1) = (/"max"/) !smooth type for species 1
		smooth_type(1:nedges,2) = (/"min"/) !smooth type for species 2
		smooth_vars(1:3) = (/1,0,1/) !density,velocity,temperature
		
		Tmin = 0.5 * qe_C
				
!  M = 10.0714285714
!  <A> = 2.5
!  Z   = 1.0
!  T1  = 1.0
!  v_piston = 84466.2798362
!  n2n1 = 3.8850945332
!  T2T1 = 32.5711744606

  		!DT ions
		Ai(1) = 2.5
		Zi(1) = 1.0 		
  		!3He ions
		Ai(2) = 3.0
		Zi(2) = 2.0  
		n_left = 2*5.38e22 * 1.e6!m-3 this is the density of DT-ice
		T_left = 1. * qe_C!J
	
n2n1 = 3.8850945332
T2T1 = 1.
		n_shocked = n_left * n2n1
T_shocked = 1.0 * qe_C * T2T1 
		
		n_right = 7.18e19 * 1.e6!m-3 this is the density of the gas
		T_right = T_left !J
		
		eps_n = 1.e-2  !should not be too much lower than 1.e-4	
		nmin = n_right * eps_n

		!-------------------------------
		!-- parameters for collisions --
		!-------------------------------
		do i = 1, nspec
			mi_g(i) = Ai(i) * mp_g; !mass of <DT> in grams
			mi(i) = 1.e-3 * mi_g(i)
			qi(i) = Zi(i) * qe; !charge of species i in statcoulomb
		enddo
		
		
		!define initial density and temperature
		!--------------------------------------
!		!DT
		r0(1,:)  = (/ double precision :: 0.0, 1.2e-2, L /)
		den0(1,:)   = (/ double precision :: 1.e26,1.e20/)
		temp0(1,:)  = (/ double precision :: 100.*qe_C,1.*qe_C/)
		vel0(1,:)     = (/ double precision :: -0.3e6, -0.3e6 /)
		
		!3He
		r0(2,:)  = (/ double precision :: 0.0, 1.2e-2, L /)
		den0(2,:)   = (/ double precision :: 1.e20,1.e25/)
		temp0(2,:)  = (/ double precision :: 100.*qe_C,1.*qe_C/)
		vel0(2,:)     = (/ double precision :: -0.3e6, -0.3e6 /)
						
		xmax_visc = 180.e-6 !for r > xmax_visc, ion_visc = ion_visc * 1.e-3

		maxind  = 21200000
		dt_multi = 1
		
		dr = L / nz !kg/m^2 !note that dr is negative
		if(geom=="spherical") then
			dr = -dr !i.e., dr is negative
		endif
		

eps_visc_max = 5.5e-3!artificial viscosity
eps_compress = 1.0
		if (geom=="slab") then
			rmin = 0.
		else !spherical
			rmin = 2.5e-6
		endif

		phi = 0.5		
		tm_quiet = 9.0e-8 
		nquiet = floor(  1.18e-2 / abs(dr) )

		cs = max( sqrt( g * T_left / minval(mi) ), sqrt( g * T_right / minval(mi) ) )

		CFL = abs(dr) / cs

		dtm = coeff_CFL * CFL !may need to be varied to avoid NaN when tm > tm_coll
		
		lm = dtm / dr

		if (ion_viscosity) then
			dtm_ = 0.5 * dtm
			eps_visc_max = 1./3 * eps_visc_max!artificial viscosity
		else
			dtm_ = dtm
		endif		

		nsmooth = abs(int(dr_smooth/dr)) !only used if smoothing==.true.
		vsmooth = nsmooth  !smoothing for velocity
		tsmooth = nsmooth  !smoothing for temperature
		
		write(*,*) "smoothing over", nsmooth, "grid points"

	end subroutine initialize

end module