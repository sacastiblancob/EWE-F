!========================================================================
! This file is part of EwE-F
! Copyright (C) 2011-2019 Middle East Technical University
! Institute of Marine Sciences (IMS-METU), Erdemli/Turkey and
! Istituto Nazionale di Oceanografia e di Geofisica Sperimentale (OGS),
! Trieste/Italy.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 as
! published by the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with  this program; if not,
! see <http://www.gnu.org/licenses/gpl-2.0.html>.
!========================================================================

subroutine rk4 (B, time, deltat, integrate, lat, lon)

#ifdef isWithBFM
  use global_mem
#else
  use statevartypesecopath, only: RLEN
#endif

!! Runge-Kutta 4th order ODE solver

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
real(RLEN)              :: dh, d6, th
real(RLEN)              :: yt(nvars),  dym(nvars), dydx(nvars), dyt(nvars)
real(RLEN)              :: loss(nvars), biomeq(nvars)
real(RLEN)              :: lossSt(nvars)

!dh = real(deltat / 2.0D0, 4)
!d6 = real(deltat / 6.0D0, 4)
dh = real(deltat / 2.0D0, RLEN)
d6 = real(deltat / 6.0D0, RLEN)
th = time + dh

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
        yt(i) = (1 - relax) * biomeq(i) + relax * B(i)
    else if (integrate(i) /= 0) then
        yt(i) = B(i) + dh * dydx(i)
    else
    end if
end do

#ifdef _Ecospace_
 call derivs (th, yt, dyt, biomeq, loss, integrate, lat, lon)
#else
 call derivs (th, yt, dyt, biomeq, loss, integrate)
#endif

do i = 1, nvars
    if (integrate(i) == 0) then
        yt(i) = (1 - relax) * biomeq(i) + relax * B(i)
    else if (integrate(i) /= 0) then
        yt(i) = B(i) + dh * dyt(i)
    else
    end if
end do

#ifdef _Ecospace_
 call derivs (th, yt, dym, biomeq, loss, integrate, lat, lon)
#else
 call derivs (th, yt, dym, biomeq, loss, integrate)
#endif

do i = 1, nvars
    if (integrate(i) == 0) then
        yt(i) = (1 - relax) * biomeq(i) + relax * B(i)
        if (yt(i) < 0) then
            print *, "Error in Runge Kutta!"
            call abort
        end if
    else if (integrate(i) /= 0) then
        yt(i) = B(i) + deltat * dym(i)
        dym(i) = dyt(i) + dym(i)
    else
    end if
end do

#ifdef _Ecospace_
 call derivs (time + deltat, yt, dyt, biomeq, loss, integrate, lat, lon)
#else
 call derivs (time + deltat, yt, dyt, biomeq, loss, integrate)
#endif

do i = 1, nvars
    loss(i) = lossSt(i) ! added later for compatibility with the EwE6
                        ! code but this did not make difference
    if (integrate(i) > 0) then
        if (integrate(i) == i) then
            B(i) = B(i) + d6 * (dydx(i) + dyt(i) + 2 * dym(i))
        end if
    else if (integrate(i) == 0) then
        B(i) = (1 - relax) * biomeq(i) + relax * B(i)
    else
        B(i) = B(i)
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

end do

if (UpdateStanzas .eqv. .true.) then

    if (nstanzas > 0) then
        ! below was lossSt but corrected upon EwE6 rk4 routine check
        call multistanza (B)
    end if

    ! Update foraging times at the end of each time step
    call updateForagingTimes (integrate)

end if

end subroutine rk4
