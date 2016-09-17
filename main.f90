! Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2014-2016

! Model parameters (pre-compiler constants) are located in shellmodels.f90 which defines the shell model to be integrated. 
! shellmodels.f90 is generated by setup.pl using the template shellmodels_plain.f90.
#include "shellmodels.f90"

module helshell
	
	implicit none
	public                  
	integer, private :: ii 
	
        ! ---------- MODEL SETUP ----------
	
	! Weights
	real,parameter,dimension(NUM_TRIAD_GEOMS)   :: Gpq = G_CONSTRUCTOR  ! Fortran is case-insensitive, we therefore rename G -> Gpq
	real,parameter,dimension(NUM_TRIAD_GEOMS,4) :: g   = g_CONSTRUCTOR
	real,parameter,dimension(NUM_TRIAD_GEOMS,4) :: Gg  = [Gpq*g(:,1), Gpq*g(:,2), Gpq*g(:,3), Gpq*g(:,4)]
	real,parameter,dimension(NUM_TRIAD_GEOMS,4) :: eps = eps_CONSTRUCTOR
	real,parameter,dimension(NUM_TRIAD_GEOMS,4) :: xi  = xi_CONSTRUCTOR
	
	! Truncation, interaction geometries, forcing and dissp
	integer, parameter      :: npad         = q_MAX            ! Number of shells as padding (ie. boundary shells having u(n) = 0 so that trid-interaction with these have no contr.)
	integer, parameter      :: nsh 	        = NSH + 2*npad	   ! Total number of shells
	real,    parameter      :: up0          = 1.0D0            ! Initial velocity shell 1 pos-amplitude
	real,    parameter      :: un0          = up0*MODEL_PARITY       
	integer, parameter,dimension(NUM_TRIAD_GEOMS) :: plist = p_LIST ! Active p,q pairs (triad geometries)
	integer, parameter,dimension(NUM_TRIAD_GEOMS) :: qlist = q_LIST

        ! k-vectors and modal interaction coeficients (labelled "d")
	real,parameter,dimension(nsh)                   :: k  = [(KZERO*LAMBDA**(ii-npad), ii=1, nsh)] ! Wave vectors, padding adjusted so that first non padding element has k = lambda^1
	real,parameter,dimension(nsh,NUM_TRIAD_GEOMS,4) :: d1 = d1_CONSTRUCTOR ! Interaction term 1 of 3
	real,parameter,dimension(nsh,NUM_TRIAD_GEOMS,4) :: d2 = d2_CONSTRUCTOR ! Interaction term 2 of 3
	real,parameter,dimension(nsh,NUM_TRIAD_GEOMS,4) :: d3 = d3_CONSTRUCTOR ! Interaction term 3 of 3
	
        ! ---------- NUMERICS ----------
	complex, dimension(nsh,4) :: dup = cmplx(0,0), dun = cmplx(0,0) ! d(up)dt, d(un)dt increments
        complex, dimension(nsh,4) ::  up = cmplx(0,0),  un = cmplx(0,0) ! Instanteneous positive and negative shell shell velocities
	integer, parameter        :: nt 	 = NT        ! Outer time-step loop, i.e. number of velocity profiles/states saved.
	integer, parameter        :: nti 	 = NTI       ! Inner time-step loop
	integer                   :: nt_prev         = 0         ! If resumed, this is the total number of steps taken so far (nt).
	integer                   :: aggr_steps_prev = 0         ! If resumed, this is the total number of agrregated steps taken so far (AGGREGATION_STEPS).  
	real,    parameter        :: dt          = DT        ! Time step
	real,    parameter        :: dt2=dt/2, dt3=dt/3, dt6=dt/6 ! RK4 time steps
	real,    parameter, dimension(nsh) :: nuk2x    = NUK2X
	integer, parameter        :: ismultithreaded = MULTITHREADED ! OpenMP flag

        ! ---------- AGGREGATES ----------
!        integer, parameter                     :: num_anomscaling = 5                   ! 
!        real, dimension(nsh,num_anomscaling)   :: E_anomscaling = 0                     ! 
        real, dimension(nsh)                   :: up_abs = 0, un_abs = 0                ! Aggregated Up_n*conj(Up_n), Un_n*conj(Un_n) fields.
        real, dimension(nsh,NUM_TRIAD_GEOMS,5) :: corr_p = 0, corr_n = 0                ! Aggregated correlators: Delta_n^{+,s',s''} and Delta_n^{-,s',s''}
        real, dimension(nsh,NUM_TRIAD_GEOMS,5) :: Ecorr = 0                             ! Instanteneous energy tripple-correlator:  Delta_{n,p,q}^{+,s',s''} + MODEL_PARITY* Delta_{n,p,q}^{-,s',s''}
        real, dimension(nsh,NUM_TRIAD_GEOMS,5) :: Eflux_instant = 0                     ! Instanteneous energi fluxes (Pi^E) for each submodel and triad geometry (p,q pair)
        real, dimension(nsh,NUM_TRIAD_GEOMS,5) :: Eflux = 0                             ! Aggregated energi fluxes (Pi^E) for each submodel and triad geometry (p,q pair)
        real, dimension(nsh,5)                 :: Eflux_instant__submodel = 0           ! Temp var used for structfunc calculations
        integer, parameter                         :: num_structfunc = 2                ! Number of calculated structure functions
        real, parameter, dimension(num_structfunc) :: structfuncs_expo = [2/3.0,4/3.0]  ! ... structfunc exponents
        real, dimension(nsh,5,num_structfunc)      :: structfuncs                       ! Aggregated structure functions: S_l(k_n) = (k_n^{-1}*abs(Eflux_n))^{1/l}. The "+1"-entries are submodel- and triad geometry summed contributions.
        
        ! ---------- SHORTCUTS ----------
	integer, parameter                        :: nshtop=nsh-npad, nshbot=npad+1       ! Top (highest) and bottom (lowest) active shells.
        integer, parameter                        :: nshfrc     = FSH + nshbot-1          ! ACTUAL Forcing shell (counting boundary shells too). FSH = 1 gives the bottom active shell.
	integer, parameter                        :: nsh_active = nsh - 2*npad	          ! Number of active shells.
	integer, parameter, dimension(nsh_active) :: rng_0      = [(ii,ii=nshbot,nshtop)] ! Full range of active shells.
	! -------------------------------
	real,parameter,dimension(nsh)               :: kinv    = k**(-1)
	real,parameter,dimension(nsh)               :: k2      = k*2
	real,parameter,dimension(NUM_TRIAD_GEOMS,4) :: GgepsSp = [Gpq*g(:,1)*eps(:,1)*Sp_1, Gpq*g(:,2)*eps(:,2)*Sp_2, Gpq*g(:,3)*eps(:,3)*Sp_3, Gpq*g(:,4)*eps(:,4)*Sp_4]

        ! ---------- SUBROUTINES ----------
	contains
	subroutine initialise() ! Initialise velocity profiles
                up(rng_0,1) = up0*cmplx(k(rng_0))**(-1.0D0/3)
        	un(rng_0,1) = un0*cmplx(k(rng_0))**(-1.0D0/3)	        
#if DISABLE_N_SHELL
        	un(rng_0,1) = 0
#endif

! Some debugging stuff...
!                up(nshbot:(nshfrc+1),1) = 1.0D1*up0*cmplx(k(nshbot:(nshfrc+1)))**(-1.0D0/3)
!	        up((nshfrc+2):nshtop,1) = 1.0D-20
!	        un(rng_0,1)             = 0
	end subroutine
end module

program helicalshellmodel

	use ieee_arithmetic
	use helshell
	implicit none
	integer                                 :: ii,jj, nn,mm,pp,qq, geom,submodel,sfunc, start ! Loop indices
	integer :: m0, m1, m2
	real                                    :: time_start, time_finish  ! CPU-time used for integration.
        complex, dimension(:,:,:), allocatable  :: saved
        integer                                 :: nt10 = nt/100;     ! for stdout pct progress bar
        integer                                 :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM ! OpenMP stuff

	call print_setup()
        allocate(saved(nsh,nt+1,2))
        
#if RESUME == 0
        call initialise()
	saved(:,1,1) = up(:,1) ! First saved entry is initial profile
	saved(:,1,2) = un(:,1)
        start = 2 ! Don't overwrite initial profile
#else
	call resume(start) ! The value of "start" is now the resumed length of the velocity history "saved"
	up(:,1) = saved(:,start,1) ! Last saved state
	un(:,1) = saved(:,start,2) 
	start = start+1 ! Set starting point of outer loop according to saved state
        write(*,"(A19,I9)"), ' *** Resuming step ', start-1
        print * 
#endif

        call cpu_time(time_start) ! Save starting time

        !$OMP PARALLEL SHARED(un,up,dun,dup,saved), PRIVATE(ii,jj,nn)
        if (ismultithreaded .eq. 1) then
                TID = OMP_GET_THREAD_NUM()
                print *, 'Thread ', TID
                if (TID .EQ. 0) then
                        NTHREADS = OMP_GET_NUM_THREADS()
                        print *, '#Threads = ', NTHREADS
                end if
        end if
        
        print *, "Starting integration..."
	do ii = start, nt+1 
		do jj = 1,nti
		
		        ! Writing out loops hints the optimizer to vectorize. This choice of loop structures below is HIGHLY intentional, do not change
		        
		        ! ---------- STEP 1 ---------------
		        !$OMP DO SCHEDULE(STATIC)
		        do nn = nshbot,nshtop 
                                dup(nn,1) = DUDT_P(1,nn)                        
#if (DTYPE == 3) && !(MULTITHREADED)
                        end do
		        do nn = nshbot,nshtop 
#endif
                                dun(nn,1) = DUDT_N(1,nn)
                        end do
                        !$OMP END DO
                        !$OMP SINGLE
!                       up(:,2) = up(:,1) + dt2*dup(:,1)  
!                       un(:,2) = un(:,1) + dt2*dun(:,1)
		        do nn = nshbot,nshtop
                		up(nn,2) = nuk2x(nn)*( up(nn,1) + dt2*dup(nn,1) )
                                un(nn,2) = nuk2x(nn)*( un(nn,1) + dt2*dun(nn,1) )
                        end do
                        !$OMP END SINGLE

		        ! ---------- STEP 2 ---------------
		        !$OMP DO SCHEDULE(STATIC)
		        do nn = nshbot,nshtop
        			dup(nn,2) = DUDT_P(2,nn)
#if (DTYPE == 3) && !(MULTITHREADED)
                        end do
		        do nn = nshbot,nshtop 
#endif
                                dun(nn,2) = DUDT_N(2,nn)
                        end do
                        !$OMP END DO
                        !$OMP SINGLE
!         		up(:,3) = up(:,1) + dt2*dup(:,2)
!                       un(:,3) = un(:,1) + dt2*dun(:,2)
		        do nn = nshbot,nshtop
                 		up(nn,3) = nuk2x(nn)*up(nn,1) + dt2*dup(nn,2)
                                un(nn,3) = nuk2x(nn)*un(nn,1) + dt2*dun(nn,2)
                        end do
                        !$OMP END SINGLE

		        ! ---------- STEP 3 ---------------
		        !$OMP DO SCHEDULE(STATIC)
		        do nn = nshbot,nshtop
        		        dup(nn,3) = DUDT_P(3,nn)
#if (DTYPE == 3) && !(MULTITHREADED)
                        end do
		        do nn = nshbot,nshtop 
#endif
                                dun(nn,3) = DUDT_N(3,nn)
                        end do
                        !$OMP END DO
                        !$OMP SINGLE
! 		        up(:,4) = up(:,1) + dt*dup(:,3)
!                       un(:,4) = un(:,1) + dt*dun(:,3)
		        do nn = nshbot,nshtop
        		        up(nn,4) = nuk2x(nn)*( nuk2x(nn)*up(nn,1) + dt*dup(nn,3) )
                                un(nn,4) = nuk2x(nn)*( nuk2x(nn)*un(nn,1) + dt*dun(nn,3) )
                        end do
                        !$OMP END SINGLE

		        ! ---------- STEP 4 ---------------
		        !$OMP DO SCHEDULE(STATIC)
		        do nn = nshbot,nshtop
        			dup(nn,4) = DUDT_P(4,nn)
#if (DTYPE == 3) && !(MULTITHREADED)
                        end do
		        do nn = nshbot,nshtop 
#endif
                                dun(nn,4) = DUDT_N(4,nn)
                        end do
                        !$OMP END DO
                        !$OMP SINGLE
!        		up(:,1) = up(:,1) + dt6*dup(:,1) + dt3*dup(:,2) + dt3*dup(:,3) + dt6*dup(:,4) ! RK4 solution
!                       un(:,1) = un(:,1) + dt6*dun(:,1) + dt3*dun(:,2) + dt3*dun(:,3) + dt6*dun(:,4) ! RK4 solution
		        do nn = nshbot,nshtop
                		up(nn,1) = nuk2x(nn)*( nuk2x(nn)*( up(nn,1) + dt6*dup(nn,1) ) + dt3*dup(nn,2) + dt3*dup(nn,3) ) + dt6*dup(nn,4) ! RK4 solution
                                un(nn,1) = nuk2x(nn)*( nuk2x(nn)*( un(nn,1) + dt6*dun(nn,1) ) + dt3*dun(nn,2) + dt3*dun(nn,3) ) + dt6*dun(nn,4) ! RK4 solution
                        end do
                        
                        ! ---------- ADD FORCING ---------------
                        
                        do nn = nshfrc,(nshfrc+DFSH-1)
	        	        up(nn,1)   = up(nn,1)   + (FORCING_P(nn)/DFSH)
                                un(nn,1)   = un(nn,1)   + (FORCING_N(nn)/DFSH)
                        end do

                        ! ---------- AGGREGATED STATS ---------------
#if DO_AGGREGATIONS
                        ! Spectra components
		        do nn = nshbot,nshtop
        		        up_abs(nn) = up_abs(nn) + up(nn,1)*CONJG(up(nn,1))
        		        un_abs(nn) = un_abs(nn) + un(nn,1)*CONJG(un(nn,1))
#if !(DISABLE_E_ANOMSCALING)
!                                do mm = 2,num_anomscaling
!                                !E_anomscaling(nn,mm) =  
!                                end do
#endif
                        end do

                        ! Correlators Delta_n^{+,s',s''}, Delta_n^{-,s',s''} and energy-like correlators Delta_{n,p,q}^{+,s',s''} + MODEL_PARITY* Delta_{n,p,q}^{-,s',s''}
                        Ecorr(:,:,:) = 0
		        do nn = nshbot,nshtop
		        
	                        PLACEHOLDER__SUM_CORR_P_1
	                        PLACEHOLDER__SUM_CORR_P_2
	                        PLACEHOLDER__SUM_CORR_P_3
	                        PLACEHOLDER__SUM_CORR_P_4

                                PLACEHOLDER__SUM_CORR_N_1
                                PLACEHOLDER__SUM_CORR_N_2
                                PLACEHOLDER__SUM_CORR_N_3
                                PLACEHOLDER__SUM_CORR_N_4
#if !(DISABLE_FLUX_CALC)
                                PLACEHOLDER__ECORR_1
                                PLACEHOLDER__ECORR_2
                                PLACEHOLDER__ECORR_3
                                PLACEHOLDER__ECORR_4
#endif
                        end do
#if !(DISABLE_FLUX_CALC)
                        ! Energy and helicity fluxes and structure functions
                        Eflux_instant__submodel(:,:) = 0
                        Eflux_instant(:,:,:) = 0;
                        do geom = 1,NUM_TRIAD_GEOMS
                        
                                pp = plist(geom)
                                qq = qlist(geom)

        		        do nn = nshbot,nshtop
		        
!		                        m0 = nn+1+((nn-qq)<0)*(qq-nn)
		                        m0 = nn+1+(((nn-npad)-qq)<0)*(qq-(nn-npad))
		                        m1 = min(nn+qq, nshtop)
		                        m2 = min(nn+qq-pp, nshtop)
                                        
                                        do mm = m0,m1
#if (SUBMODEL_CLASS==5)
                                                Eflux_instant(nn,geom,1) = Eflux_instant(nn,geom,1) + Gg(geom,1) * k2(mm-qq) * Ecorr(mm-1,geom,1) 
                                                Eflux_instant(nn,geom,2) = Eflux_instant(nn,geom,2) + Gg(geom,2) * k2(mm-qq) * Ecorr(mm-1,geom,2) 
                                                Eflux_instant(nn,geom,3) = Eflux_instant(nn,geom,3) + Gg(geom,3) * k2(mm-qq) * Ecorr(mm-1,geom,3) 
                                                Eflux_instant(nn,geom,4) = Eflux_instant(nn,geom,4) + Gg(geom,4) * k2(mm-qq) * Ecorr(mm-1,geom,4) 
#else
                                                Eflux_instant(nn,geom,SUBMODEL_CLASS) = Eflux_instant(nn,geom,SUBMODEL_CLASS) + Gg(geom,SUBMODEL_CLASS) * k2(mm-qq) * Ecorr(mm-1,geom,SUBMODEL_CLASS) 
#endif
                                        end do
                                        
                                        do mm = m0,m2
#if (SUBMODEL_CLASS==5)
                                                Eflux_instant(nn,geom,1) = Eflux_instant(nn,geom,1) - GgepsSp(geom,1) * k2(mm-qq) * Ecorr(mm-1,geom,1) 
                                                Eflux_instant(nn,geom,2) = Eflux_instant(nn,geom,2) - GgepsSp(geom,2) * k2(mm-qq) * Ecorr(mm-1,geom,2) 
                                                Eflux_instant(nn,geom,3) = Eflux_instant(nn,geom,3) - GgepsSp(geom,3) * k2(mm-qq) * Ecorr(mm-1,geom,3) 
                                                Eflux_instant(nn,geom,4) = Eflux_instant(nn,geom,4) - GgepsSp(geom,4) * k2(mm-qq) * Ecorr(mm-1,geom,4) 
#else
                                                Eflux_instant(nn,geom,SUBMODEL_CLASS) = Eflux_instant(nn,geom,SUBMODEL_CLASS) - GgepsSp(geom,SUBMODEL_CLASS) * k2(mm-qq) * Ecorr(mm-1,geom,SUBMODEL_CLASS) 
#endif
                                        end do

                                        Eflux(nn,geom,:) = Eflux(nn,geom,:) + Eflux_instant(nn,geom,:)

                                end do
#if !(DISABLE_STRUCTURE_FUNCS)
                                Eflux_instant__submodel(:,:) = Eflux_instant__submodel(:,:) + Eflux_instant(:,geom,:)
#endif

                        end do
                        
#if !(DISABLE_STRUCTURE_FUNCS)
                        Eflux_instant__submodel(:,5)   = abs(sum(Eflux_instant__submodel(:,1:4),2))
                        Eflux_instant__submodel(:,1:4) = abs(Eflux_instant__submodel(:,1:4))
                        
                        do sfunc = 1,num_structfunc       
#if (SUBMODEL_CLASS==5)
                                structfuncs(:,1,sfunc) = structfuncs(:,1,sfunc) + ( kinv*Eflux_instant__submodel(:,1) )**(structfuncs_expo(sfunc))
                                structfuncs(:,2,sfunc) = structfuncs(:,2,sfunc) + ( kinv*Eflux_instant__submodel(:,2) )**(structfuncs_expo(sfunc))
                                structfuncs(:,3,sfunc) = structfuncs(:,3,sfunc) + ( kinv*Eflux_instant__submodel(:,3) )**(structfuncs_expo(sfunc))
                                structfuncs(:,4,sfunc) = structfuncs(:,4,sfunc) + ( kinv*Eflux_instant__submodel(:,4) )**(structfuncs_expo(sfunc))
                                structfuncs(:,5,sfunc) = structfuncs(:,5,sfunc) + ( kinv*Eflux_instant__submodel(:,5) )**(structfuncs_expo(sfunc))
#else
                                structfuncs(:,SUBMODEL_CLASS,sfunc) = structfuncs(:,SUBMODEL_CLASS,sfunc) + ( kinv*Eflux_instant__submodel(:,SUBMODEL_CLASS) )**(structfuncs_expo(sfunc))
#endif
                        end do
#endif                              
#endif
#endif
                        !$OMP END SINGLE
		end do

                !$OMP SINGLE
		saved(:,ii,1) = up(:,1) ! Save state
		saved(:,ii,2) = un(:,1)
		if (MOD(ii-1,nt10) == 0) then
		        write(*,"(I3,A)"), INT(FLOAT(ii-1)/nt * 100.), '%'
		if (ieee_is_nan( REAL(up(nshbot,1)) )) then
				write(*,"(A45)"), '*** NaNs encountered, halting integration.'
				saved(:,nt+1,1) = saved(:,ii-1,1)
				saved(:,nt+1,2) = saved(:,ii-1,2)
				exit
	        end if
		end if
		!$OMP END SINGLE
		
	end do 

        !$OMP SINGLE
        call cpu_time(time_finish) ! Save starting time
	print *
        write(*,"(A34,F12.3)"), 'Integration CPU-time in seconds: ',time_finish-time_start
	print *
        call dump()
	!$OMP END SINGLE
	!$OMP END PARALLEL		

	contains 
#include "dump.f90"
	subroutine print_setup()
		print *, '-------------- MODEL ------------------------'
	        write(*,"(A19,I2)"),                     'MODEL: ',             MODEL
	        write(*,"(A19,I3,A2,I2,A1)"),            'SHELLS (PADDING): ',  NSH, ' (',npad,')'
	        write(*,"(A19,F4.2)"),                   'LAMBDA: ',            LAMBDA
	        write(*,"(A19,I2)"),                     'INTERACTION SCOPE: ', q_MAX
	        print *, '-------------- FORCING/DISSIPATION ----------'
	        write(*,"(A19,A5,I2,A12,I3)"),           'FORCING: ',           'Type ',FTYPE,' at shell n=',nshfrc
                write(*,"(A19,A5,I1)"),                  'DISSIPATION: ',       'Type ',DTYPE
	        write(*,"(A19,E7.1)"),                   'NU (SMALL SCALE): ',  VISC_SMALLSCALE
	        write(*,"(A19,E7.1)"),                   'NU (LARGE SCALE): ',  VISC_LARGESCALE
	        print *, '-------------- NUMERICS ---------------------'
	        write(*,"(A19,E7.1)"),                   'DT: ',                dt
	        write(*,"(A19,I5,A10,I5,A16)"),          'TSTEPS: ',            nt/1000, 'E+3 (with ',nti/1000,'E+3 inner steps)'
	        print *, '---------------------------------------------'
                if (DISABLE_VEL_SAMPLING .eq. 1)    write(*,"(A71)"),'*** Warning: Saving only end members of the velocity profiles sampled.'
                if (.not.DO_AGGREGATIONS .eq. 1)    write(*,"(A40)"),'*** Warning: Omitting aggregated stats.'
                if (DISABLE_STRUCTURE_FUNCS .eq. 1) write(*,"(A43)"),'*** Warning: Omitting structure functions.'
                if (DISABLE_FLUX_CALC .eq. 1)       write(*,"(A41)"),'*** Warning: Omitting flux calculations.'
                if (DISABLE_N_SHELL .eq. 1)         write(*,"(A40)"),'*** Warning: Negative shell disabled.'
                print *        
	end subroutine

end program helicalshellmodel
