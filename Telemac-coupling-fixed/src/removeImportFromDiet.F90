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

              subroutine removeImportFromDiet &
  (nvars, ep_diet, ep_data, es_data)

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

!use statevartypesecopath, only: ep_data, ep_diet
use statevartypesecopath, only: ecopath_data
!use statevartypesecosim, only: es_data, nvars
use statevartypesecosim, only: ecosim_data

implicit none

! INTENT VARIABLES
INTEGER, INTENT(IN)               :: nvars
TYPE(ecopath_data), INTENT(IN)    :: ep_data(nvars)
REAL(RLEN), INTENT(INOUT)         :: ep_diet &
                             (nvars+1,count(ep_data(:)%org_type == 2))
TYPE(ecosim_data), INTENT(INOUT)  :: es_data(nvars)

! IN SUBROUTINE VARIABLES
integer :: iec, j
real(RLEN) :: fractionW0import

!!

do j = 1, nvars
    if (ep_data(j)%org_type == 2) then
        if (ep_diet(nvars + 1, j) > 0) then
            fractionW0import = (1 - ep_diet(nvars + 1, j) / 1)
        else
            fractionW0import = 1
        end if

        do iec = 1, nvars
            if (fractionW0import == 0) then
                ep_diet(iec, j) = 0
            else
                ep_diet(iec, j) = ep_diet(iec, j)
            end if
        end do
!        write(*,*) ep_diet(:,:)
        es_data(j)%QBoutside = ep_data(j)%QoB * (1 - fractionW0Import)
!        write(*,*) j, es_data(j)%QBoutside
    end if
end do

end subroutine removeImportFromDiet
