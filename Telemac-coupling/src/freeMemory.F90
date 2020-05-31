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

subroutine freeMemory

  use statevartypesecopath
  use statevartypesecosim
!#ifdef _Ecospace_
    use statevartypesecospace
!#endif

    deallocate(ep_data)
    if ( allocated(ms_data) ) deallocate(ms_data)
    if ( allocated(es_ms_data) ) deallocate (es_ms_data)

!#ifdef _Ecospace_
    deallocate(QperB)
    deallocate(M2)
    deallocate(BB_spatial)
    deallocate(advection)
    deallocate(grid)
    deallocate(spatialhafs)
!#endif

    deallocate(det_export_rate)
    deallocate(detritus_no)
    deallocate(flow2detritus)
    deallocate(BBAvg)
    deallocate(LossAvg)
    deallocate(EatenByAvg)
    deallocate(EatenOfAvg)
    deallocate(PredAvg)
    deallocate(BB)

!#ifdef isWithBFM
!    deallocate(ruHTLc)
!    deallocate(ruHTLn)
!    deallocate(ruHTLp)
!    deallocate(ruHTLl)
!    deallocate(ruHTLi)
!    deallocate(ruHTLs)
!    deallocate(tfluxc)
!    deallocate(tfluxn)
!    deallocate(tfluxp)
!    deallocate(p_qncHTL)
!    deallocate(p_qpcHTL)
!    deallocate(flow2detritusR6c)
!    deallocate(flow2detritusR6n)
!    deallocate(flow2detritusR6p)
!#endif

    deallocate(arena%a)
    deallocate(arena%Q_link)
    deallocate(arena%Q_arena)
    deallocate(arena%vul_arena)
    deallocate(arena%vul_biom)
    deallocate(arena%vulrate)
    deallocate(arena%a_link)
    deallocate(arena%base_time_switch)
    deallocate(arena%rela_switch)
    deallocate(arena%a_eff)
    deallocate(arena%v_eff)
    deallocate(arena%v_biom)
    deallocate(arena%v_denom)

    if ( allocated(producer_no) ) deallocate(producer_no)
    if ( allocated(PrimaryProdForce) ) deallocate(PrimaryProdForce)
    if ( allocated(NutrientForce) ) deallocate(NutrientForce)

    deallocate(es_vul)
    deallocate(force%forcetype)
    deallocate(force%fishforce)
    deallocate(ep_diet)
    deallocate(ep_detfate)
    deallocate(es_data)
    deallocate(NutFreeBase)

  end subroutine freeMemory
