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

subroutine readSpatialDistribution_io ()

  use statevartypesecosim, only: nvars
  use statevartypesecospace

  implicit none

  integer :: i, j
  character(len = 30) :: fileinptr

  nvars = 22;
  allocate(spatialhafs(nlat, nlon, nvars))

  do i = 1, nvars
      if (i < 10) then
          write(fileinptr, '(A20, I1, ".csv")') SpatialDistribution_dirname, i
      else
          write(fileinptr, '(A20, I2, ".csv")') SpatialDistribution_dirname, i
      end if
      
      print *, fileinptr
      open(4440444, file = fileinptr, form = 'formatted')
      read(4440444, *)
      read(4440444, *)
      
      do j = 1, nlat
          read(4440444, *, end = 1000) spatialhafs(j, :, i)
      end do

1000  continue
  close(4440444)

  end do

end subroutine readSpatialDistribution_io
