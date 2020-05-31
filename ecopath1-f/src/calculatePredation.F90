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

subroutine calculatePredation (dcols, order, group_rank, var)

! This subroutine calculates the predatory consumption
! on each state variable.

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, ep_diet, m_consumed, nvars

implicit none

! variables inherited from Ecopath-FORTRAN model
integer, intent(in)  :: dcols
integer, intent(in)  :: var
integer, intent(in)  :: order(nvars), group_rank(nvars)

! in-subroutine variable declarations
real(RLEN) :: total_consumption            ! total predation on group
real(RLEN) :: consumption                  ! sum of predations
integer :: j                            ! in-subroutine loop index
  

if (group_rank(order(var)) == dcols) then ! if no predators
    m_consumed(order(var)) = 0

else
    if (ep_data(order(var))%org_type == 2) then ! if a consumer

        if (ep_diet(order(var),order(var)) == 0) then
            ! if no intraguild predation
            total_consumption = 0

            do j = 1, dcols
                if (ep_data(j)%biomass /= -999) then
                    consumption = ep_data(j)%biomass * ep_data(j)%QoB &
                      * ep_diet(order(var), j)
                    total_consumption = total_consumption + consumption
                end if
            end do
            
            m_consumed(order(var)) = total_consumption   
        else ! there is intraguild predation            
            if (ep_data(order(var))%biomass == -999 .or. &
              ep_data(order(var))%QoB == -999) then 
                ! either biomass or QoB is to be estimated
                total_consumption = 0

                do j = 1, dcols
                    if (j == order(var)) then
                        ! skip the group itself due to either biomass
                        ! or QoB parameter is to be estimated
                        cycle
                    end if
                    if (ep_data(j)%biomass /= -999) then
                        consumption = ep_data(j)%biomass * ep_data(j)%QoB &
                          * ep_diet(order(var), j)
                        total_consumption = total_consumption + consumption
                    end if
                end do

                m_consumed(order(var)) = total_consumption
            else ! otherwise biomass and QoB of group(i) are both known
                total_consumption = 0

                do j = 1, dcols
                    if (ep_data(j)%biomass /= -999) then
                        consumption = ep_data(j)%biomass * ep_data(j)%QoB &
                          * ep_diet(order(var), j)
                        total_consumption = total_consumption + consumption
                    end if
                end do

                m_consumed(order(var)) = total_consumption
            end if  
        end if
    else ! then this group is either detritus or producer

        total_consumption = 0

        do j = 1, dcols
            if (ep_data(j)%biomass /= -999) then
                consumption = ep_data(j)%biomass * ep_data(j)%QoB &
                  * ep_diet(order(var), j)
                total_consumption = total_consumption + consumption
            end if
        end do

        m_consumed(order(var)) = total_consumption
    end if
                                                      
end if

end subroutine calculatePredation



