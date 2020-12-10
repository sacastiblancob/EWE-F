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

subroutine readDetritusFate_io (detrows, detcols, ndetritus)

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: DetFate_fname, isASCIIinputFile, ep_detfate

implicit none

integer, intent(inout) :: detrows, detcols, ndetritus
integer                :: i            ! loop counter
character ( len = 250 ) :: groupname, dummy

! open detritus groups' information file
open(444, file = DetFate_fname, form = 'formatted')

! read number of rows and columns in the file
read(444, *) dummy, detrows
read(444, *) dummy, detcols
read(444, *) dummy, ndetritus
read(444, *)

! allocate size of Ecopath data in variable <m_data>
if (isASCIIinputFile) then
    allocate(ep_detfate(detrows, detcols))
end if

! read Ecopath input data from file
if (isASCIIinputFile) then
    do i = 1, detrows
        read(444,*,end=400) groupname, ep_detfate(i, :)
    end do
end if

400 continue
 close(444)

end subroutine readDetritusFate_io
