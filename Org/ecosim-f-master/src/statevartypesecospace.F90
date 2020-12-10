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

module statevartypesecospace

#ifdef isWithBFM
  use global_mem
#endif

#ifndef isWithBFM
  use statevartypesecopath, only: RLEN
#endif

  implicit none

  integer                 :: nlat, nlon
  real(RLEN), allocatable :: spatialhafs (:, :, :)
  real(RLEN), allocatable :: BB_spatial (:, :, :)
  real(RLEN), allocatable :: QperB (:, :, :)
  real(RLEN), allocatable :: M2 (:, :, :)
  integer, allocatable    :: grid(:, :)
  integer, allocatable    :: advection(:, :)

  character(len = 250) :: SpatialGrid_fname
  character(len = 250) :: SpatialDistribution_dirname
  character(len = 250) :: Advection_fname

end module statevartypesecospace


