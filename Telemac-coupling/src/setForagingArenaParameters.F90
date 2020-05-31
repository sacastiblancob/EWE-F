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

subroutine setForagingArenaParameters (vrows, vcols)

! This subroutine populates the fields of foraging arena parameters &
! <arena> that are calculated within the Ecosim model. The type of the
! columns of the <arena> matrix is defined in 
! statevartypesEcosim-mod.f95. However, they can be extended by 
! appending new fields to <arena>.

use statevartypesecosim, only: arena

implicit none

! variables inherited from ecopath model
integer, intent(in) :: vrows, vcols

! allocate the fields in memory
allocate(arena%a(vrows, vcols))
allocate(arena%Q_link(vrows, vcols))
allocate(arena%Q_arena(vrows, vcols))
allocate(arena%vul_arena(vrows, vcols))
allocate(arena%vul_biom(vrows, vcols))
allocate(arena%vulrate(vrows, vcols))
allocate(arena%a_link(vrows, vcols))
allocate(arena%base_time_switch(vrows, vcols))
allocate(arena%rela_switch(vrows, vcols))
allocate(arena%a_eff(vrows, vcols))
allocate(arena%v_eff(vrows, vcols))
allocate(arena%v_biom(vrows, vcols))
allocate(arena%v_denom(vrows, vcols))

end subroutine setForagingArenaParameters
