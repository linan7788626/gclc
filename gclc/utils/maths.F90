module gclc_maths

    use healpix_types

    implicit none

    contains

    subroutine radec2thetaphi(ra, dec, theta, phi)

        implicit none

        real(DP),   intent(in)  :: ra, dec
        real(DP),   intent(out) :: theta, phi

        theta = HALFPI - DEG2RAD*dec
        phi   = DEG2RAD*ra
        
        return
    end subroutine
end module
