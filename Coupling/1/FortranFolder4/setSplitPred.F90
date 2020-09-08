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

subroutine setSplitPred (B)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

use statevartypesecopath, only: ms_data
use statevartypesecosim, only: es_data, es_ms_data, nvars, nstanzas

implicit none

! variables inherited
real(RLEN), intent(inout) :: B(nvars)

! in-subroutine variables
real(RLEN) :: Bt        ! B calc. by multiplying numbers by weight at age
real(RLEN) :: Pt        ! Q calculated by multiplying numbers by Q at age
real(RLEN) :: Nt        ! numbers at age
integer :: age       ! age (in months)
integer :: stanza    ! stanza number
integer :: substanza ! substanza number
integer :: age_last  ! maximum age of substanza group

do stanza = 1, nstanzas
!    write(*,*) stanza,nstanzas
!    write(*,*)
    do substanza = 1, ms_data(stanza)%substanzas

        Bt = 1.0e-30
        Pt = 1.0e-30
        Nt = 1.0e-30

        if  (substanza < ms_data(stanza)%substanzas) then
            age_last = ms_data(stanza)%age_start(substanza + 1) - 1
        else
            age_last = ms_data(stanza)%age_infinity
        end if

!        write(*,*) age_last
!        write(*,*)

        do age = ms_data(stanza)%age_start(substanza), age_last
            Bt = Bt + es_ms_data(stanza)%NageS(age) &
              * es_ms_data(stanza)%WageS(age)
            Pt = Pt + es_ms_data(stanza)%NageS(age) &
              * ms_data(stanza)%WWa(age)
            Nt = Nt + es_ms_data(stanza)%NageS(age)
        end do

        B(ms_data(stanza)%ep_groupno(substanza))            = Bt
        es_data(ms_data(stanza)%ep_groupno(substanza))%pred = Pt

!        write(*,*) Bt
!        write(*,*) ms_data(stanza)%ep_groupno(substanza)
!        write(*,*) Pt
!        write(*,*) ms_data(stanza)%ep_groupno(substanza)
!        write(*,*) Nt
!        write(*,*)


    end do
end do

end subroutine setSplitPred
