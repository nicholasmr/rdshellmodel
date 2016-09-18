! Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2014-2016

! This file is a template. It must be parsed by the setup.pl pre-compiler to replace all "PREPLACE__*" symbols, which generates the proper shellmodels.f90 file used by main.f90. The shellmodels.f90 file then contains the specific model setup wanted. 
! NOTE: capital "lc" is a marker for the pre-compiler to substitute for a newline character. See the Makefile.

!-----------------------------------
! FLAGS AND DEBUGGING
!-----------------------------------
#define DISABLE_N_SHELL         0       /* Disable u^-_n (negative) shells*/
#define DISABLE_VEL_SAMPLING    1       /* If disabled, only first and last simulated velocities (u_n^{+,-}) are saved */

#define DISABLE_FLUX_CALC       1       /* If doing aggreations, skip Pi_^E_n (Eflux) */
#define DISABLE_STRUCTURE_FUNCS 1       /* If doing aggreations, skip structure functions (structfuncs) */
#define DISABLE_E_ANOMSCALING   1       /* NOT YET IMPLEMENTED */

!-----------------------------------
! FORCING TYPE
!-----------------------------------
! If ftype is positive, both helical shells (+,-) are forced, if negative, only positive shell(s) are forced (+).
! ftype=0 => no forcing, ftype=1 or -1 => f_n = const, ftype=2 or -2 => f_n = const/conj(u_n)

#if FTYPE == 0
#define FORCING_0   0
#else
#define FORCING_0   FMAGNITUDE*cmplx(1, 1)
#endif

#define FORCING__CONST_U_INCREMENT           dt*FORCING_0
#define FORCING__CONST_E_INCREMENT__P(N)     dt*FORCING_0/CONJG(up(N,1))
#define FORCING__CONST_E_INCREMENT__N(N)     dt*FORCING_0/CONJG(un(N,1))

#if FTYPE == 0
#define FORCING_P(N)  0
#define FORCING_N(N)  0
#elif FTYPE == 1
#define FORCING_P(N)  FORCING__CONST_U_INCREMENT
#define FORCING_N(N)  FORCING__CONST_U_INCREMENT
#elif FTYPE == -1
#define FORCING_P(N)  FORCING__CONST_U_INCREMENT
#define FORCING_N(N)  0
#elif FTYPE == 2
#define FORCING_P(N)  FORCING__CONST_E_INCREMENT__P(N)
#define FORCING_N(N)  FORCING__CONST_E_INCREMENT__N(N)
#elif FTYPE == -2
#define FORCING_P(N)  FORCING__CONST_E_INCREMENT__P(N)
#define FORCING_N(N)  0
#endif

#if DISABLE_N_SHELL
#define FORCING_N 0
#endif

!-----------------------------------
! DISSIPATION TYPE
!-----------------------------------
! dtype=0 => no dissipation, dtype=1 => small-scale dissipation, dtype=2 => large-scale dissipation, dtype=3 => large- and small-scale dissipation

! Uses the "direct integration of dissipation trick" (nuk2x) 
#define NUK2X__SMALLSCALE   -dt2*VISC_SMALLSCALE*k(ii)**2
#define NUK2X__LARGESCALE   -dt2*VISC_LARGESCALE*k(ii)**(LSKDEP)

#if DTYPE == 0
#define NUK2X                   [(1, ii=1, nsh)]
#elif DTYPE == 1
#define NUK2X                   [(exp(NUK2X__SMALLSCALE), ii=1, nsh)] 
#elif DTYPE == 2
#define NUK2X                   [(exp(NUK2X__LARGESCALE), ii=1, nsh)] 
#elif DTYPE == 3
#define NUK2X                   [(exp(NUK2X__SMALLSCALE NUK2X__LARGESCALE), ii=1, nsh)] 
#endif

!-----------------------------------
! MODEL TYPES
!-----------------------------------
! Model 1-4    Helical Sabra models    (unpublished master thesis: N. M. Rathmann (2014), De Pietro et al., 2015) 
! Model 11-14  New helical models      (Rathmann and Ditlevsen, 2016)
! Model 15     New mixed helical model (Rathmann and Ditlevsen, 2016)
!-----------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! CUSTOM PRE-COMPILER STUFF WHICH IS REPLACED BY setup.pl
!----------------------------------------------------------------------
!----------------------------------------------------------------------
#define MODEL                   REPLACE__MODEL
#define NSH                     REPLACE__NSH
#define LAMBDA                  REPLACE__LAMBDA
#define q_MAX                   REPLACE__q_MAX
#define NUM_TRIAD_GEOMS         REPLACE__NUM_TRIAD_GEOMS
#define SUBMODEL_CLASS          REPLACE__SUBMODEL_CLASS
#define q_LIST                  REPLACE__q_LIST  
#define p_LIST                  REPLACE__p_LIST
! Interaction coefficients (array constructors)
#define G_CONSTRUCTOR           REPLACE__G_CONSTRUCTOR
#define g_CONSTRUCTOR           REPLACE__g_CONSTRUCTOR
#define eps_CONSTRUCTOR         REPLACE__eps_CONSTRUCTOR
#define xi_CONSTRUCTOR          REPLACE__xi_CONSTRUCTOR
#define d1_CONSTRUCTOR          REPLACE__d1_CONSTRUCTOR
#define d2_CONSTRUCTOR          REPLACE__d2_CONSTRUCTOR
#define d3_CONSTRUCTOR          REPLACE__d3_CONSTRUCTOR

! PARITY HERE REFERS TO THE SIGN SYMMETRY WITHIN THE MODEL EQUATIONS AND TRIPLE-CORRELATORS
#if MODEL < 10 
#define MODEL_PARITY  +1
#else
#define MODEL_PARITY  -1
#endif

!-----------------------------------
! Non-linear sigma-interaction models (sigma_1, ..., sigma_6) (Ditlevsen, 2010)
!-----------------------------------
! -1, +1, -1, -1, +1, -1  (Model 1)
! -1, -1, -1, +1, -1, +1  (Model 2)
! +1, -1, +1, -1, -1, -1  (Model 3)
! +1, +1, +1, +1, +1, +1  (Model 4)
!-----------------------------------

! Model 1-4 {s',s''} pairs
#define Sp_1 -1
#define Sq_1 +1
#define Sp_2 -1
#define Sq_2 -1
#define Sp_3 +1
#define Sq_3 -1
#define Sp_4 +1
#define Sq_4 +1

#define IDX_1(N,p,q) N+p
#define IDX_2(N,p,q) N+q
#define IDX_3(N,p,q) N-p
#define IDX_4(N,p,q) N+q-p
#define IDX_5(N,p,q) N-q
#define IDX_6(N,p,q) N+p-q

!-----------------------------------

#define DUDT_P_1(STEP,N,triad,p,q)  d1(N,triad,1)*CONJG(un(IDX_1(N,p,q),STEP))*up(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,1)*CONJG(un(IDX_3(N,p,q),STEP))*un(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,1)*      up(IDX_5(N,p,q),STEP) *un(IDX_6(N,p,q),STEP) 

#define DUDT_N_1(STEP,N,triad,p,q)  d1(N,triad,1)*CONJG(up(IDX_1(N,p,q),STEP))*un(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,1)*CONJG(up(IDX_3(N,p,q),STEP))*up(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,1)*      un(IDX_5(N,p,q),STEP) *up(IDX_6(N,p,q),STEP) 
				  
!-----------------------------------
#define DUDT_P_2(STEP,N,triad,p,q)  d1(N,triad,2)*CONJG(un(IDX_1(N,p,q),STEP))*un(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,2)*CONJG(un(IDX_3(N,p,q),STEP))*up(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,2)*      un(IDX_5(N,p,q),STEP) *up(IDX_6(N,p,q),STEP) 

#define DUDT_N_2(STEP,N,triad,p,q)  d1(N,triad,2)*CONJG(up(IDX_1(N,p,q),STEP))*up(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,2)*CONJG(up(IDX_3(N,p,q),STEP))*un(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,2)*      up(IDX_5(N,p,q),STEP) *un(IDX_6(N,p,q),STEP) 
				  
!-----------------------------------
#define DUDT_P_3(STEP,N,triad,p,q)  d1(N,triad,3)*CONJG(up(IDX_1(N,p,q),STEP))*un(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,3)*CONJG(up(IDX_3(N,p,q),STEP))*un(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,3)*      un(IDX_5(N,p,q),STEP) *un(IDX_6(N,p,q),STEP)        

#define DUDT_N_3(STEP,N,triad,p,q)  d1(N,triad,3)*CONJG(un(IDX_1(N,p,q),STEP))*up(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,3)*CONJG(un(IDX_3(N,p,q),STEP))*up(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,3)*      up(IDX_5(N,p,q),STEP) *up(IDX_6(N,p,q),STEP)        
				  
!-----------------------------------
#define DUDT_P_4(STEP,N,triad,p,q)  d1(N,triad,4)*CONJG(up(IDX_1(N,p,q),STEP))*up(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,4)*CONJG(up(IDX_3(N,p,q),STEP))*up(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,4)*      up(IDX_5(N,p,q),STEP) *up(IDX_6(N,p,q),STEP)        

#define DUDT_N_4(STEP,N,triad,p,q)  d1(N,triad,4)*CONJG(un(IDX_1(N,p,q),STEP))*un(IDX_2(N,p,q),STEP) + LC\
        			    d2(N,triad,4)*CONJG(un(IDX_3(N,p,q),STEP))*un(IDX_4(N,p,q),STEP) + LC\
				    d3(N,triad,4)*      un(IDX_5(N,p,q),STEP) *un(IDX_6(N,p,q),STEP)        

!-----------------------------------

#define DUDT_P(STEP,N)  REPLACE__DUDT_P 
#if !(DISABLE_N_SHELL)
#define DUDT_N(STEP,N)  REPLACE__DUDT_N 
#else
#define DUDT_N(STEP,N)  0
#endif

!-----------------------------------
! AGGREGATED TRIPLE CORRELATIONS
!-----------------------------------

#define CORR_P_1(N,p,q)        REAL(CONJG(up(IDX_5(N,p,q)+1,1))*CONJG(un(IDX_6(N,p,q)+1,1))*up(N+1,1))
#define CORR_N_1(N,p,q)        REAL(CONJG(un(IDX_5(N,p,q)+1,1))*CONJG(up(IDX_6(N,p,q)+1,1))*un(N+1,1))

#define CORR_P_2(N,p,q)        REAL(CONJG(up(IDX_5(N,p,q)+1,1))*CONJG(un(IDX_6(N,p,q)+1,1))*un(N+1,1))
#define CORR_N_2(N,p,q)        REAL(CONJG(un(IDX_5(N,p,q)+1,1))*CONJG(up(IDX_6(N,p,q)+1,1))*up(N+1,1))

#define CORR_P_3(N,p,q)        REAL(CONJG(up(IDX_5(N,p,q)+1,1))*CONJG(up(IDX_6(N,p,q)+1,1))*un(N+1,1))
#define CORR_N_3(N,p,q)        REAL(CONJG(un(IDX_5(N,p,q)+1,1))*CONJG(un(IDX_6(N,p,q)+1,1))*up(N+1,1))

#define CORR_P_4(N,p,q)        REAL(CONJG(up(IDX_5(N,p,q)+1,1))*CONJG(up(IDX_6(N,p,q)+1,1))*up(N+1,1))
#define CORR_N_4(N,p,q)        REAL(CONJG(un(IDX_5(N,p,q)+1,1))*CONJG(un(IDX_6(N,p,q)+1,1))*un(N+1,1))

! These placeholders are just 
#define PLACEHOLDER__SUM_CORR_P_1  REPLACE__SUM_CORR_P_1
#define PLACEHOLDER__SUM_CORR_P_2  REPLACE__SUM_CORR_P_2
#define PLACEHOLDER__SUM_CORR_P_3  REPLACE__SUM_CORR_P_3
#define PLACEHOLDER__SUM_CORR_P_4  REPLACE__SUM_CORR_P_4
#define PLACEHOLDER__SUM_CORR_N_1  REPLACE__SUM_CORR_N_1
#define PLACEHOLDER__SUM_CORR_N_2  REPLACE__SUM_CORR_N_2
#define PLACEHOLDER__SUM_CORR_N_3  REPLACE__SUM_CORR_N_3
#define PLACEHOLDER__SUM_CORR_N_4  REPLACE__SUM_CORR_N_4
#define PLACEHOLDER__ECORR_1       REPLACE__ECORR_1
#define PLACEHOLDER__ECORR_2       REPLACE__ECORR_2
#define PLACEHOLDER__ECORR_3       REPLACE__ECORR_3
#define PLACEHOLDER__ECORR_4       REPLACE__ECORR_4


