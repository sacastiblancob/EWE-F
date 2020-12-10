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

program ecopath

! This is the main ecopath model. All subroutines are called from
! this model and state equations are calculated here.

use statevartypesecopath
use readHDF5Database

implicit none

integer :: drows, dcols     ! # of rows & columns of diet data
integer :: detrows, detcols ! # of rows & columns of detritus data
integer :: detritus_no      ! detritus number being handled
integer :: var              ! loop counter

! rank of groups in terms of prey-predator relationship
integer, allocatable :: group_rank(:)

! the precedence order of group equations to be solved
integer, allocatable :: order(:)

! the index of detritus groups
integer, allocatable :: detritus_group_no(:)

! numbers of producer & detritus groups
integer              :: nproducer, ndetritus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! get file names from namelist (nml) file
namelist /filenames/ InputData_fname, DietComp_fname, DetFate_fname, &
  GrowthParam_fname, Results_fname, HDF5_fname, isASCIIinputFile
open(1010, file = "filenames.nml", status = 'OLD')
read(1010, nml = filenames)
 close(1010)

!!!!!!!!!! SECTION: READING INPUT DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call readEcopathScenario_io ()

!!!!!!!!!! END SECTION: READING INPUT DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: READING DETRITUS INFORMATION !!!!!!!!!!!!!!!!!!!!!!!!

 call readDetritusFate_io (detrows, detcols, ndetritus)

!!!!!!!!!! END SECTION: READING DETRITUS INFORMATION !!!!!!!!!!!!!!!!!!!!

if (isASCIIinputFile .eqv. .false.) then
    allocate(ep_growth(nstanzas))
    call readHDF5 ()
end if

!!!!!!!!!! SECTION: MULTISTANZA CALCULATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nstanzas > 0) then
  allocate(ms_data(nstanzas))
  call calculateMultistanzaParameters ()
end if
!!!!!!!!!! END SECTION: MULTISTANZA CALCULATIONS !!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: READING DIET COMPOSITION MATRIX !!!!!!!!!!!!!!!!!!!!!

 call readDietComposition_io (drows, dcols)

!!!!!!!!!! END SECTION: READING DIET COMPOSITION MATRIX !!!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: DETERMINE PREY-PREDATOR RANK ORDER !!!!!!!!!!!!!!!!!!

allocate(detritus_group_no(ndetritus))
allocate(group_rank(nvars))
allocate(order(nvars))

 call calculateSolvingOrderOfLinearEquations (nproducer, ndetritus &
             , dcols, order, group_rank, detritus_group_no)

!!!!!!!!!! END SECTION: DETERMINE PREY-PREDATOR RANK ORDER !!!!!!!!!!!!!!

!! estimate Q/B, P/Q or P/B
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then ! if group is a consumer
        if (ep_data(var)%PoB /= -999) then ! if P/B is available
            if (ep_data(var)%PoQ /= -999) then ! if P/Q is available
                ! estimate Q/B
                ep_data(var)%QoB = ep_data(var)%PoB / ep_data(var)%PoQ
            else if (ep_data(var)%QoB /= -999) then ! if Q/B is available
                ! estimate P/Q
                ep_data(var)%PoQ = ep_data(var)%PoB / ep_data(var)%QoB
            else
                print *, "One of Q/B and P/Q should be entered for group", var, "!"
                print *, "The program will halt!"
                print *, "ERROR in model!!!!!!!!!!!!"
                call abort
            end if
        else if (ep_data(var)%QoB /= -999 .and. &
          ep_data(var)%PoQ /= -999) then
            ! otherwise P/B is to be estimated from Q/B and P/Q
            ep_data(var)%PoB = ep_data(var)%QoB * ep_data(var)%PoQ
        else ! then P/B and one of Q/B and P/Q are unknown
            continue 
        end if
    end if
end do

!!!!!!!!!! SECTION: LINEAR EQUATIONS ARE SOLVED !!!!!!!!!!!!!!!!!!!!!!!!!
allocate(m_flows2det(nvars, ndetritus))
allocate(m_consumed(nvars))
m_flows2det = 0
m_consumed = 0
detritus_no = 0

do var = 1, nvars

    if (ep_data(order(var))%PoB /= -999) then
        if (ep_data(order(var))%QoB /= -999 &
          .or. ep_data(order(var))%PoQ /= -999 &
          .or. ep_data(order(var))%org_type == 1) then
            if (ep_data(order(var))%EE /= -999) then
                call calculatePredation (dcols, order, group_rank, var)

                if (ep_data(order(var))%org_type == 2) then

                    if (ep_diet(order(var),order(var)) == 0) then
                        ! no intraguild predation
                        ep_data(order(var))%biomass &
                          = (ep_data(order(var))%landings &
                          + ep_data(order(var))%discards &
                          + m_consumed(order(var))) & 
                          / (ep_data(order(var))%PoB &
                          * ep_data(order(var))%EE)

                    else ! intraguild predation
                        ep_data(order(var))%biomass &
                          = (ep_data(order(var))%landings &
                          + ep_data(order(var))%discards &
                          + m_consumed(order(var))) & 
                          / (ep_data(order(var))%PoB &
                          - ep_data(order(var))%QoB &
                          * ep_diet(order(var),order(var)) &
                          - (1 - ep_data(order(var))%EE) &
                          * ep_data(order(var))%PoB)
                    end if

                else
                    ep_data(order(var))%biomass &
                      = (ep_data(order(var))%landings &
                      + ep_data(order(var))%discards &
                      + m_consumed(order(var))) &
                      / (ep_data(order(var))%PoB &
                      * ep_data(order(var))%EE)
                end if

            end if
        end if
    end if

    if (ep_data(order(var))%biomass /= -999) then
        if (ep_data(order(var))%org_type  /= 0) then ! if not detritus
            if (ep_data(order(var))%PoB /= -999) then
                if (ep_data(order(var))%org_type  /= 1) then
                    ! if not producer
                    if (ep_data(order(var))%QoB /= -999 .or. &
                      ep_data(order(var))%PoQ /= -999) then
                        call calculatePredation (dcols, order, &
                          group_rank, var)

                        ep_data(order(var))%EE &
                          = (ep_data(order(var))%landings &
                          + ep_data(order(var))%discards &
                          + ep_data(order(var))%BA &
                          + m_consumed(order(var))) & 
                          / (ep_data(order(var))%PoB &
                          * ep_data(order(var))%biomass)
                    end if

                else
                    call calculatePredation (dcols, order, group_rank, &
                      var)

                    ep_data(order(var))%EE &
                      = (ep_data(order(var))%landings &
                      + ep_data(order(var))%discards &
                      + ep_data(order(var))%BA &
                      + m_consumed(order(var))) & 
                      / (ep_data(order(var))%PoB &
                      * ep_data(order(var))%biomass)
                end if
            end if
        else ! then detritus
            detritus_no = detritus_no + 1
            call calculatePredation (dcols, order, group_rank, var)
            call calculateDetritalFlows (detritus_no)

        end if
    end if
end do

 call calculateDetritusFate (ndetritus, detritus_group_no)
 call calculateBAofDetritus (ndetritus, detritus_group_no)
 call calculateEEofDetritus (ndetritus, detritus_group_no)


!!!!!!!!!! END SECTION: LINEAR EQUATIONS ARE SOLVED !!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! SECTION: ENERGY BALANCE CALCULATIONS !!!!!!!!!!!!!!!!!!!!!!!!!

 call calculateEnergyBalance ()

!!!!!!!!!! END SECTION: ENERGY BALANCE CALCULATIONS !!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: RESULTS OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1111, file = 'Ecopath_Results.dat', form = 'formatted', &
  status = 'unknown')

write(1111, 13) "Group/Parameter", "Biomass", "P/B", "Q/B", "EE", "P/Q", &
  "Unass. Q", "Detr. Import", "Landings", "Discards", "Org. Type", &
  "Production", "Consumption", "Respiration", "Assimilation", &
  "Flows2Detritus"
13 format(a17,",", 9(a21,","), a12, 5(a21,","))

do var = 1, nvars
    write(1111, 15) groupnames(var), ep_data(var)%biomass, ep_data(var)%PoB, &
      ep_data(var)%QoB, ep_data(var)%EE, ep_data(var)%PoQ, &
      ep_data(var)%unass_Q, ep_data(var)%detritus_import, &
      ep_data(var)%landings, ep_data(var)%discards, ep_data(var)%org_type, &
      ep_data(var)%production, ep_data(var)%consumption, &
      ep_data(var)%respiration, ep_data(var)%assimilation, &
      sum(m_flows2det(var, :))
end do
15 format(A35,",", 9(f27.9,","), i12, 5(f27.9,","))

 close(1111)
 
 ! store results in HDF5 format
 call writeHDF5ResultsFile (drows, dcols, detrows, detcols)

!!!!!!!!!! END SECTION: RESULTS OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: SANITY CHECKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call warnSanityChecks ()

!!!!!!!!!! END SECTION: SANITY CHECKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if ( allocated(ep_growth) ) deallocate(ep_growth)
 if ( allocated(ms_data) ) deallocate(ms_data)
 if ( allocated(ep_diet) ) deallocate(ep_diet)
 if ( allocated(ep_detfate) ) deallocate(ep_detfate)

 deallocate(ep_data)
 deallocate(detritus_group_no)
 deallocate(group_rank)
 deallocate(order)
 deallocate(m_flows2det)
 deallocate(m_consumed)
 
end program ecopath
