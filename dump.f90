! Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2014-2016

! See for netCDF routine examples see: http://www.unidata.ucar.edu/software/netcdf/examples/programs/

subroutine dump()

        use netcdf
        implicit none

        character (len = *), parameter :: FILE_NAME = "./dump.nc"
        integer :: ID_nc, ID_up_re, ID_up_im, ID_un_re, ID_un_im, ID_corr_p, ID_corr_n, ID_up_abs, ID_un_abs, ID_Eflux, ID_structfuncs
        integer :: ID_dim_time, ID_dim_shell, ID_dim_model, ID_dim_p, ID_dim_q, ID_dim_triad, ID_dim_structfunc, ID_dims(2), ID_dims3(3), ID_dims4(4)
        
#if DISABLE_VEL_SAMPLING == 0
#define SAVESLICE    :
#define SAVESLICELEN nt+1
#else
#define SAVESLICE    [1,nt+1]
#define SAVESLICELEN 2
#endif

        real, dimension(nsh,SAVESLICELEN) :: up_re, up_im, un_re, un_im
        
        up_re(:,:) = REAL( saved(:,SAVESLICE,1))
        up_im(:,:) = AIMAG(saved(:,SAVESLICE,1))
        un_re(:,:) = REAL( saved(:,SAVESLICE,2))
        un_im(:,:) = AIMAG(saved(:,SAVESLICE,2))
        
print *, "*** NETCDF: creating file"
        call check( nf90_create(FILE_NAME, NF90_CLOBBER, ID_nc) )
        
print *, "*** NETCDF: saving model setup"
        ! Model
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "model",   MODEL) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "lambda",  LAMBDA) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "parity",  MODEL_PARITY) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "visc",    VISC_SMALLSCALE) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "viscinv", VISC_LARGESCALE) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "kdepinv", LSKDEP) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "k",       k) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "k0",      KZERO) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "nsh",     nsh) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "nshtop",  nshtop) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "nshbot",  nshbot) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "q_max",   q_MAX) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "p_list",  p_LIST) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "q_list",  q_LIST) )
        
        ! Forcing/dissipation
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "fsh",     FSH) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "ftype",   FTYPE) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "dtype",   DTYPE) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "up_0",    up0) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "un_0",    un0) )

        ! Weights (dims flattened)      
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "eps",   RESHAPE(eps, [SIZE(eps)]) ) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "xi",    RESHAPE(xi,  [SIZE(xi)]) ) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "g",     RESHAPE(g,   [SIZE(g)]) ) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "G",     RESHAPE(Gpq, [SIZE(Gpq)]) ) )
        
        ! Numerics
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "nt",      nt) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "nti",     nti) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "dt",      dt) )
        
        ! Other
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "DO_AGGREGATIONS",      DO_AGGREGATIONS) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "RESUME_AGGREGATIONS",  RESUME_AGGREGATIONS) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "AGGREGATION_STEPS",    DO_AGGREGATIONS*(nt-nt_prev) + RESUME_AGGREGATIONS*aggr_steps_prev) )
        call check( nf90_put_att(ID_nc, NF90_GLOBAL, "DISABLE_N_SHELL",      DISABLE_N_SHELL) )
        
print *, "*** NETCDF: defining dimensions"
!        call check( nf90_def_dim(ID_nc, "time",   nt+1, ID_dim_time) )
        call check( nf90_def_dim(ID_nc, "time",  NF90_UNLIMITED, ID_dim_time) )
        call check( nf90_def_dim(ID_nc, "shell", nsh,            ID_dim_shell) )
        call check( nf90_def_dim(ID_nc, "model", 5,              ID_dim_model) ) ! 5'th entry is for multi-submodel aggregates in coupled configurations.
        call check( nf90_def_dim(ID_nc, "p",     q_MAX,          ID_dim_p) )
        call check( nf90_def_dim(ID_nc, "q",     q_MAX,          ID_dim_q) )
        call check( nf90_def_dim(ID_nc, "triad", NUM_TRIAD_GEOMS,ID_dim_triad) )
        call check( nf90_def_dim(ID_nc, "structfuncs", num_structfunc, ID_dim_structfunc) )
        
print *, "*** NETCDF: defining variables"
        ID_dims =  (/ ID_dim_shell, ID_dim_time /)
        call check( nf90_def_var(ID_nc, "up_re", NF90_DOUBLE, ID_dims, ID_up_re) )
        call check( nf90_def_var(ID_nc, "up_im", NF90_DOUBLE, ID_dims, ID_up_im) )
        call check( nf90_def_var(ID_nc, "un_re", NF90_DOUBLE, ID_dims, ID_un_re) )
        call check( nf90_def_var(ID_nc, "un_im", NF90_DOUBLE, ID_dims, ID_un_im) )
        
        ID_dims = (/ ID_dim_shell, ID_dim_model/)
        call check( nf90_def_var(ID_nc, "up_abs_sum", NF90_DOUBLE, ID_dim_shell, ID_up_abs) )
        call check( nf90_def_var(ID_nc, "un_abs_sum", NF90_DOUBLE, ID_dim_shell, ID_un_abs) )

        ID_dims3 = (/ ID_dim_shell, ID_dim_triad, ID_dim_model/)
        call check( nf90_def_var(ID_nc, "corr_p_sum", NF90_DOUBLE, ID_dims3, ID_corr_p) )
        call check( nf90_def_var(ID_nc, "corr_n_sum", NF90_DOUBLE, ID_dims3, ID_corr_n) )
        call check( nf90_def_var(ID_nc, "Eflux_sum", NF90_DOUBLE, ID_dims3, ID_Eflux) )

        ID_dims3 = (/ ID_dim_shell, ID_dim_model, ID_dim_structfunc/)
        call check( nf90_def_var(ID_nc, "structfuncs_sum", NF90_DOUBLE, ID_dims3, ID_structfuncs) )
        
        call check( nf90_enddef(ID_nc) )
        
print *, "*** NETCDF: saving velocity variables"
        call check( nf90_put_var(ID_nc, ID_up_re, up_re) )
        call check( nf90_put_var(ID_nc, ID_up_im, up_im) )
        call check( nf90_put_var(ID_nc, ID_un_re, un_re) )
        call check( nf90_put_var(ID_nc, ID_un_im, un_im) )
        call check( nf90_put_var(ID_nc, ID_corr_p, corr_p) )
        call check( nf90_put_var(ID_nc, ID_corr_n, corr_n) )
        call check( nf90_put_var(ID_nc, ID_up_abs, up_abs) )
        call check( nf90_put_var(ID_nc, ID_un_abs, un_abs) )
        call check( nf90_put_var(ID_nc, ID_Eflux, Eflux) )
        call check( nf90_put_var(ID_nc, ID_structfuncs, structfuncs) )
        
        call check( nf90_close(ID_nc) )
        print *, "*** SUCCEEDED dumping to netCDF"

end subroutine

subroutine resume(len_time)

        use netcdf
        implicit none

        integer, intent (out) :: len_time
        character (len = *), parameter :: FILE_NAME = "./resume.nc"
        integer :: ID_nc, ID_up_re, ID_up_im, ID_un_re, ID_un_im, ID_corr_p, ID_corr_n, ID_up_abs, ID_un_abs, ID_Eflux, ID_structfuncs
        integer :: ID_dim_time
        real, dimension(:,:), allocatable :: up_re, up_im, un_re, un_im
       
        print *, "*** NETCDF: Loading file to resume from"
        call check( nf90_open(FILE_NAME, NF90_NOWRITE, ID_nc) )

        ! Get record dimension length
        call check( nf90_inq_dimid(ID_nc, "time", ID_dim_time) )
        call check( nf90_inquire_dimension(ID_nc, ID_dim_time, len = len_time) )

        ! Get var IDs
        call check( nf90_inq_varid(ID_nc, "up_re", ID_up_re) )
        call check( nf90_inq_varid(ID_nc, "up_im", ID_up_im) )
        call check( nf90_inq_varid(ID_nc, "un_re", ID_un_re) )
        call check( nf90_inq_varid(ID_nc, "un_im", ID_un_im) )
        call check( nf90_inq_varid(ID_nc, "corr_p_sum", ID_corr_p) )
        call check( nf90_inq_varid(ID_nc, "corr_n_sum", ID_corr_n) )
        call check( nf90_inq_varid(ID_nc, "up_abs_sum", ID_up_abs) )
        call check( nf90_inq_varid(ID_nc, "un_abs_sum", ID_un_abs) )
        call check( nf90_inq_varid(ID_nc, "Eflux_sum", ID_Eflux) )
        call check( nf90_inq_varid(ID_nc, "structfuncs_sum", ID_structfuncs) )

        ! Read state history
        allocate(up_re(nsh,len_time))
        allocate(up_im(nsh,len_time))
        allocate(un_re(nsh,len_time))
        allocate(un_im(nsh,len_time))
        call check( nf90_get_var(ID_nc, ID_up_re, up_re) )
        call check( nf90_get_var(ID_nc, ID_up_im, up_im) )
        call check( nf90_get_var(ID_nc, ID_un_re, un_re) )
        call check( nf90_get_var(ID_nc, ID_un_im, un_im) )
#if RESUME_AGGREGATIONS
        print *, "*** Resuming aggregations"
        call check( nf90_get_var(ID_nc, ID_corr_p, corr_p) )
        call check( nf90_get_var(ID_nc, ID_corr_n, corr_n) )
        call check( nf90_get_var(ID_nc, ID_up_abs, up_abs) )
        call check( nf90_get_var(ID_nc, ID_un_abs, un_abs) )
        call check( nf90_get_var(ID_nc, ID_Eflux, Eflux) )
        call check( nf90_get_var(ID_nc, ID_structfuncs, structfuncs) )
#else
        print *, "*** NOT resuming from saved aggregations"
#endif
        ! Attributes
        call check( nf90_get_att(ID_nc, NF90_GLOBAL, "nt", nt_prev) ) ! Assume when resuming that nti is preserved between runs.
        call check( nf90_get_att(ID_nc, NF90_GLOBAL, "AGGREGATION_STEPS", aggr_steps_prev) )

        ! Load into model structures
        if (len_time .eq. 2) then ! Resuming state without fully save velocity profile? 
                len_time = nt_prev+1
                saved(:,[1, len_time],1) = cmplx(up_re,up_im)
                saved(:,[1, len_time],2) = cmplx(un_re,un_im)                  
        else
                saved(:,1:len_time,1) = cmplx(up_re,up_im)
                saved(:,1:len_time,2) = cmplx(un_re,un_im)  
        end if
        
end subroutine

subroutine check(status)
    
        use netcdf
        implicit none
    
        integer, intent (in) :: status
    
        if (status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
        end if
end subroutine  

