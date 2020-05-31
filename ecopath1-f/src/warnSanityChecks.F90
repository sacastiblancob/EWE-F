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

subroutine warnSanityChecks ()

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, nvars

implicit none

integer :: var

do var = 1, nvars
    if (ep_data(var)%EE > 1) then
        print *, "EE is bigger than unity for group:", var
        print *, "Consider to adjust input parameters accordingly!"
        print *, "ERROR in ECOPATH!"
        call abort
    end if
end do

do var = 1, nvars
    if (ep_data(var)%respiration < 0) then
        print *, "Respiration cannot be negative for group:", var
        print *, "Consider to adjust input parameters accordingly!"
        print *, "ERROR in ECOPATH!"
        call abort
    end if
end do

end subroutine warnSanityChecks
