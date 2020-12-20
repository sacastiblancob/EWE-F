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

              subroutine readPrimaryProdForcingFunction_io &
          (PrimaryProdForcing_fname, nvars, frows, fcols)

! This subroutine reads forcing time series data
! from tab-delimited ascii file.

use statevartypesecopath, only: ep_data
use statevartypesecosim, only: PrimaryProdForce, producer_no, nproducer

implicit none

!INTENT VARIABLES
integer, INTENT(IN)              :: nvars
INTEGER, INTENT(INOUT)           :: frows, fcols
character(len = 250), INTENT(IN) :: PrimaryProdForcing_fname

! IN SUBROUTINE VARIABLES
integer              :: iec, j          ! loop counter
!integer              :: frows, fcols  ! number of rows and columns
character ( len = 250 ) :: groupname, dummy

! open time-series data file
open(555, file = PrimaryProdForcing_fname, form = 'formatted')

! read number of group, time-series type,
! number of rows and columns in file
read(555, *) dummy, frows
read(555, *) dummy, fcols
read(555, *) dummy, groupname

nproducer = count(ep_data(:)%org_type == 1)
allocate(producer_no(nproducer))
j = 0
do iec = 1, nvars
    if (ep_data(iec)%org_type == 1) then
        j = j + 1
        producer_no(j) = iec
    end if
end do

! allocate the time-series to their corresponding groups
allocate(PrimaryProdForce(frows, fcols))

! read time-series data from file
do iec = 1, frows
    read(555, *, end = 500) dummy, PrimaryProdForce(iec, :)
end do

500 continue
 close(555)

end subroutine readPrimaryProdForcingFunction_io
