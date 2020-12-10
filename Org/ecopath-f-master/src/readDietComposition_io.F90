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

subroutine readDietComposition_io (drows, dcols)

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: DietComp_fname, isASCIIinputFile, ep_diet

implicit none

integer, intent(inout) :: drows, dcols
integer                :: i
character ( len = 250 ) :: groupname, dummy

! open Ecopath diet matrix file
open(555, file = DietComp_fname, form = 'formatted')

! read number of rows and columns in the file
read(555, *) dummy, drows
read(555, *) dummy, dcols
read(555, *) 

! allocate size of Ecopath diet data in variable <ep_diet>
if (isASCIIinputFile) then
    allocate(ep_diet(drows, dcols))
end if 

! read Ecopath diet data from file
if (isASCIIinputFile) then
    do i = 1, drows
        read(555, *, end = 500) groupname, ep_diet(i, :)
    end do
end if

500 continue
 close(555)

end subroutine readDietComposition_io
