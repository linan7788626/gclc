module cmb_tools
    
    use healpix_types
    use rngmod,     only: planck_rng, rand_init

    implicit none

    complex(DPC),   allocatable     :: md_alms_TGC(:,:,:)

    contains

    subroutine create_map_unlensed(input_cls_file, seed, alms_TGC, map_TQU, theta, phi, lmax)

        implicit none

        character(*),   intent(in)  :: input_cls_file
        integer(I4B),   intent(in)  :: lmax
        integer(I4B),   intent(in)  :: seed(1:)
        real(DP),       intent(out) :: map_TQU(0:,1:)
        real(DP),       intent(in)  :: theta(0:), phi(0:)

        complex(DPC),   intent(out) :: alms_TGC(1:,0:)

    end subroutine


    subroutine synalms(input_cls_file, seed, lmax, alms_TGC)
        
        implicit none

        character(*),   intent(in)  :: input_cls_file
        integer(I4B),   intent(in)  :: seed(1:)
        integer(I4B),   intent(in)  :: lmax
        complex(DPC),   intent(out) :: alms_TGC(1:,0:)
    
        character(80)       :: header(1:60)
        type(planck_rng)    :: rng_handle
        integer(I4B)        :: np, p

        np = size(alms_TGC(:,0))

        if (np .eq. 3) then
            p = 1
        else if (np .eq. 1) then
            p = 0
        else
            write(*,*) "Check the dimension of ALMS_TGC"
            stop
        endif

        if (.not. allocated(md_alms_TGC)) then
            allocate(md_alms_TGC(1:np,0:lmax,0:lmax))
        else
            if (size(md_alms_TGC(:,0,0)) .ne. np) then
                write(*,*) 'WARNING: The dimension of allocaed MD_ALMS_TGC does NOT match ALMS_TGC'
                write(*,*) 'Now reallocate MD_ALMS_TGC'
                deallocate(md_alms_TGC)
                allocate(md_alms_TGC(1:np,0:lmax,0:lmax))
            endif
        endif

        header = ' '

        call rand_init(rng_handle, seed(1), seed(2), seed(3), seed(4))
        call create_alm(8192, lmax, lmax, p, trim(input_cls_file), rng_handle, 0.00_dp, md_alms_TGC, header)
        call convert_alms_f2cxx(md_alms_TGC, alms_TGC)

        return
    end subroutine


    subroutine convert_alms_f2cxx(alms_f, alms_cxx)

        implicit none

        complex(DPC),   intent(in)  :: alms_f(1:,0:,0:)
        complex(DPC),   intent(out) :: alms_cxx(1:,0:)

        integer(I4B)    :: il, im, idx
        integer(I4B)    :: lmax

        lmax = size(alms_f(1,:,0))

        if (lmax .ne. size(alms_f(1,0,:))) then
            write(*,*) "CONVERT_ALMS_F2CXX in cmb.F90"
            write(*,*) "Check the dimension of ALMS_F"
            stop
        endif
        
        do il=0, lmax
            do im=0, lmax
                idx = im*(2*lmax+1-im)/2+il
                alms_cxx(1:,idx) = alms_f(1:,il,im)
            enddo
        enddo

        return
    end subroutine

end module
