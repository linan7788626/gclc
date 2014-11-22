module gclc_sims_cmb

    use healpix_types
    use paramfile_io

    implicit none

    type gclc_sc_params

        integer(I4B)    :: nside_ra, nside_dec
        real(DP)        :: scale_ra, scale_dec
        real(DP)        :: ra0, dec0
        real(DP)        :: reso_arcmin

        logical         :: has_pol = .true.
        logical         :: use_nfw = .true.

        character(512)  :: input_cls_fits
        integer(I4B)    :: lmax

    end type

    type gclc_sc_data

        real(DP),       allocatable :: ra(:), dec(:)
        real(DP),       allocatable :: theta(:), phi(:)
        real(DP),       allocatable :: defle_theta(:), defle_phi(:)
        
        real(DP),       allocatable :: input_cls(:,:)
        real(DP),       allocatable :: map_TQU(:,:)
        complex(DPC),   allocatable :: alms_TEB(:,:,:)

        real(DP),       allocatable :: map_kappa(:), map_defle(:)
    end type

    type(gclc_sc_params)            :: scinfo
    type(gclc_sc_data),     target  :: scdata

    contains

    subroutine gclc_sc_init(ini_file)

        implicit none

        character(*),       intent(in)  :: ini_file
        type(paramfile_handle)          :: ini_handle

        ini_handle = parse_init(ini_file)

        scinfo%has_pol      = parse_lgt(ini_handle, 'has_pol', .true.)
        scinfo%use_nfw      = parse_lgt(ini_handle, 'use_nfw', .true.)

        scinfo%lmax         = parse_long(ini_handle, 'lmax')

        scinfo%nside_ra     = parse_long(ini_handle, 'Nside_RA')
        scinfo%nside_dec    = parse_long(ini_handle, 'Nside_DEC')

        scinfo%ra0          = parse_double(ini_handle, 'RA0')
        scinfo%dec0         = parse_double(ini_handle, 'DEC0')

        scinfo%reso_arcmin  = parse_double(ini_handle, 'reso_arcmin', 0.00_dp)
        scinfo%scale_ra     = parse_double(ini_handle, 'scale_RA', 0.00_dp)
        scinfo%scale_dec    = parse_double(ini_handle, 'scale_dec', 0.00_dp)

        if (scinfo%reso_arcmin .gt. 1.00d-6 .and. &
        &   scinfo%scale_ra*scinfo%scale_dec .gt. 1.00d-6) then
            write(*,*) 'Both reso_arcmin and scale_RA/DEC are defined.'
            write(*,*) 'Please choose one.'
            stop
        endif

        scinfo%input_cls_fits   = parse_string(ini_handle, 'input_Cls')

        return
    end subroutine

end module
