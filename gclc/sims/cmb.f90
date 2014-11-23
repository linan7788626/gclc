module gclc_sims_cmb

    use healpix_types
    use paramfile_io

    use gclc_maths

    implicit none

    type sc_params

        integer(I4B)    :: nside_ra, nside_dec
        real(DP)        :: scale_ra, scale_dec
        real(DP)        :: ra0, dec0

        logical         :: has_pol = .true.
        logical         :: use_nfw = .true.

        character(512)  :: input_cls_fits
        integer(I4B)    :: lmax
        
        logical         :: init = .false.
    end type

    type sc_data

        real(DP),       allocatable :: ra(:), dec(:)
        real(DP),       allocatable :: theta(:), phi(:)
        real(DP),       allocatable :: def_theta(:), def_phi(:)
        
        real(DP),       allocatable :: input_cls(:,:)
        real(DP),       allocatable :: map_TQU(:,:)
        complex(DPC),   allocatable :: alms_TEB(:,:,:)

        real(DP),       allocatable :: map_kappa(:,:), map_def(:,:)

        logical                     :: alloc = .false.
    end type

    type(sc_params)            :: scinfo
    type(sc_data),     target  :: scdata

    contains

    subroutine sc_init(ini_file)

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


    subroutine sc_allocate()

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
        lmax = scinfo%lmaxDEG2RAD*dec

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

        return
    end subroutine


    subroutine sc_setup_coord()
        
        integer(I4B)        :: npix, ipix
        integer(I4B)        :: ipix_ra, ipix_dec

        real(DP)            :: reso_ra, scale_ra, ra0
        real(DP)            :: reso_dec, scale_dec, dec0

        real(DP)            :: ra, dec
        real(DP)            :: theta, phi
    
        npix = scinfo%nside_ra * scinfo%nside_dec

        reso_ra = scinfo%scale_ra / scinfo%nside_ra
        reso_dec = scinfo%scale_dec / scinfo%nside_dec

        ra0 = scinfo%ra0
        dec0 = scinfo%dec0

        scale_ra = scinfo%scale_ra
        scale_dec = scinfo%scale_dec

        do ipix=0, nside_ra-1
            scinfo%ra(ipix) = scinfo%ra0 + 0.50_dp*reso_ra + ipix*reso_ra - 0.50_dp*scale_ra
        enddo

        do ipix=0, nside_dec-1
            scinfo%dec(ipix) = scinfo%dec0 + 0.50_dp*reso_dec + ipix*reso_dec - 0.50_dp*scale_dec
        enddo
        
        ipix = 0
        do ipix_ra=0, nside_ra-1
            ra = scinfo%ra(ipix_ra)
            do ipix_dec=0, nside_dec-1
                dec = scinfo%dec(ipix_dec)

                call radec2thetaphi(ra, dec, theta, phi)

                scdata%theta(ipix) = theta
                scdata%phi(ipix)   = phi
            enddo
        enddo

        return
    end subroutine


    subroutine sc_read()

        implicit none

    end subroutine

end module
