!========================================================================
! This file based on EwE-F
!
!========================================================================

SUBROUTINE EULER (B, time, deltat, integrate, lat, lon)

use statevartypesecopath, only: RLEN

use statevartypesecosim, only: es_data, nvars, nstanzas, relax, &
                         StepsPerMonth, UpdateStanzas, BBAvg, LossAvg, &
                         EatenByAvg, EatenOfAvg, PredAvg

implicit none

! variables inherited from Ecosim-F model
real(RLEN), intent(inout)     :: B(nvars)
real(RLEN), intent(in)        :: time, deltat
integer, intent(in)           :: integrate(nvars)
integer, optional, intent(in) :: lat, lon

! in-subroutine variables
integer                 :: i
real(RLEN)              :: dydx(nvars)
real(RLEN)              :: loss(nvars), biomeq(nvars)
real(RLEN)              :: lossSt(nvars)

#ifdef _Ecospace_
 call derivs (time, B, dydx, biomeq, loss, integrate, lat, lon)
#else
 call derivs (time, B, dydx, biomeq, loss, integrate)
#endif

do i = 1, nvars
    EatenByAvg(i)  = EatenByAvg(i) + es_data(i)%qq
    EatenOfAvg(i)  = EatenOfAvg(i) + es_data(i)%M2
    PredAvg(i)     = PredAvg(i) + es_data(i)%pred
    lossSt(i)      = loss(i)
end do

do i = 1, nvars
    if (integrate(i) == 0) then
        B(i) = (1 - relax) * biomeq(i) + relax * B(i)
    else if (integrate(i) /= 0) then
        B(i) = B(i) + deltat * dydx(i)
    else
    end if
end do

if (StepsPerMonth > 1) then

#ifdef _Ecospace_
    call derivs (time + deltat, B, dydx, biomeq, loss, &
            integrate, lat, lon)
#else
    call derivs (time + deltat, B, dydx, biomeq, loss, &
            integrate)
#endif

    do i = 1, nvars
        loss(i) = lossSt(i)
    end do

end if

! averaging here for multistanza calculations
! and updating foraging times
do i = 1, nvars
    
    BBAvg(i)   = BBAvg(i) + B(i)
    LossAvg(i) = LossAvg(i) + loss(i)

    if (UpdateStanzas .eqv. .true.) then
        BBAvg(i)        = BBAvg(i) / StepsPerMonth
        LossAvg(i)      = LossAvg(i) / StepsPerMonth
        EatenByAvg(i)   = EatenByAvg(i) / StepsPerMonth
        EatenOfAvg(i)   = EatenOfAvg(i) / StepsPerMonth
        PredAvg(i)      = PredAvg(i) / StepsPerMonth
        lossSt(i)       = lossSt(i) / StepsPerMonth
    end if
    !WRITE(*,*) 'BB ',B(i)
end do

if (UpdateStanzas .eqv. .true.) then
    
    if (nstanzas > 0) then
        ! below was lossSt but corrected upon EwE6 rk4 routine check
        call multistanza (B)
    end if

    ! Update foraging times at the end of each time step
    call updateForagingTimes (integrate)

end if

END SUBROUTINE EULER







