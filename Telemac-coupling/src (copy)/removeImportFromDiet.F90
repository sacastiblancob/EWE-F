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

subroutine removeImportFromDiet ()

!#ifdef isWithBFM
!  use global_mem
!#else
  use statevartypesecopath, only: RLEN
!#endif

use statevartypesecopath, only: ep_data, ep_diet
use statevartypesecosim, only: es_data, nvars

implicit none

integer :: i, j
real(RLEN) :: fractionW0import

do j = 1, nvars
    if (ep_data(j)%org_type == 2) then
        if (ep_diet(nvars + 1, j) > 0) then
            fractionW0import = (1 - ep_diet(nvars + 1, j) / 1)
        else
            fractionW0import = 1
        end if

        do i = 1, nvars
            if (fractionW0import == 0) then
                ep_diet(i, j) = 0
            else
                ep_diet(i, j) = ep_diet(i, j)
            end if
        end do
!        write(*,*) ep_diet(:,:)
        es_data(j)%QBoutside = ep_data(j)%QoB * (1 - fractionW0Import)
!        write(*,*) j, es_data(j)%QBoutside
    end if
end do

end subroutine removeImportFromDiet
