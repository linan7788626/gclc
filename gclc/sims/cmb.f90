module gclc_sims_cmb

    use healpix_types
    use paramfile_io

    implicit none

    type gclc_sc_params

        integer(I4B)    :: nside_ra, nside_dec
        real(DP)        :: scale_ra, scale_dec
        real(DP)        :: ra0, dec0

        logical         :: has_pol = .true.
        logical         :: use_nfw = .true.

        character(512)  :: input_cls_fits
        integer(I4B)    :: lmax
        
        logical         :: init = .false.
    end type

    type gclc_sc_data

        real(DP),       allocatable :: ra(:), dec(:)
        real(DP),       allocatable :: theta(:), phi(:)
        real(DP),       allocatable :: def_theta(:), def_phi(:)
        
        real(DP),       allocatable :: input_cls(:,:)
        real(DP),       allocatable :: map_TQU(:,:)
        complex(DPC),   allocatable :: alms_TEB(:,:,:)

        real(DP),       allocatable :: map_kappa(:,:), map_def(:,:)

        logical                     :: alloc = .false.
    end type

    type(gclc_sc_params)            :: scinfo
    type(gclc_sc_data),     target  :: scdata

    contains

    subroutine gclc_sc_init(ini_file)

        implicit none

        character(*),       intent(in)  :: ini_file
        type(paramfile_handle)          :: ini_handle

        logical         :: do_reso = .false.
        logical         :: do_scale = .false.

        ini_handle = parse_init(ini_file)

        scinfo%has_pol      = parse_lgt(ini_handle, 'has_pol', .true.)
        scinfo%use_nfw      = parse_lgt(ini_handle, 'use_nfw', .true.)

        scinfo%lmax         = parse_long(ini_handle, 'lmax')

        scinfo%nside_ra     = parse_long(ini_handle, 'Nside_RA')
        scinfo%nside_dec    = parse_long(ini_handle, 'Nside_DEC')

        scinfo%ra0          = parse_double(ini_handle, 'RA0')
        scinfo%dec0         = parse_double(ini_handle, 'DEC0')

        if (scinfo%reso_arcmin .gt. 1.00d-6) do_reso = .true.

        scinfo%scale_ra     = parse_double(ini_handle, 'scale_RA', 0.00_dp)
        scinfo%scale_dec    = parse_double(ini_handle, 'scale_dec', 0.00_dp)

        scinfo%input_cls_fits   = parse_string(ini_handle, 'input_Cls')

        scinfo%init = .true.

        return
    end subroutine


    subroutine gclc_sc_allocate()

        implicit none
        
        integer(I4B)        :: nside_ra, nside_dec
        integer(I4B)        :: npix, lmax
        integer(I4B)        :: np

        if (scdata%alloc) then
            return
        endif
        
        nside_ra = scinfo%nside_ra
        nside_dec = scinfo%nside_dec

        npix = nside_ra * nside_dec
        lmax = scinfo%lmax

        allocate(scdata%ra(0:nside_ra-1))
        allocate(scdata%dec(0:nside_dec-1))

        allocate(scdata%theta(0:npix-1))
        allocate(scdata%phi(0:npix-1))

        allocate(scdata%def_theta(0:npix-1))
        allocate(scdata%def_phi(0:npix-1))

        if (scinfo%has_pol) then np = 3 else np = 1
        allocate(scdata%input_cls(0:lmax,1:np))

        allocate(scdata%map_TQU(0:npix-1,1:np))
        allocate(scdata%alms_TEB(1:np,lmax*(lmax+1)/2+lmax+1))

        allocate(scdata%map_kappa(0:npix-1,1:1))
        allocate(scdata%map_def(0:npix-1,1:1))

        scdata%alloc = .true.

    end subroutine


    subroutine gclc_sc_read()

        implicit none

        
    end subroutine

end module
