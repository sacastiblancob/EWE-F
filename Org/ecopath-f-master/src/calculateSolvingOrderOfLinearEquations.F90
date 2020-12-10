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

subroutine calculateSolvingOrderOfLinearEquations (nproducer, ndetritus &
             , dcols, order, group_rank, detritus_group_no)

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: ep_data, ep_diet, nvars

implicit none

integer, intent(in)    :: dcols, ndetritus
integer, intent(inout) :: nproducer
integer, intent(inout) :: group_rank(nvars)
integer, intent(inout) :: order(nvars)
integer, intent(inout) :: detritus_group_no(ndetritus)

integer, allocatable :: int_matrix(:)      ! trans. of group_rank matrix
integer              :: var
integer              :: ind, m, n, i, j 
integer              :: counter, counter2  ! index to order precedence
integer              :: lastgroups

! determine the number of producer groups
nproducer = count(ep_data(:)%org_type == 1)

! count number of zeros in diet matrix for each group
! more zeros for group(i) means less predators
counter = 0
counter2 = 0
do var = 1, nvars
    if (ep_data(var)%org_type == 2) then
        counter = counter + 1
        group_rank(counter) = count(ep_diet(var, :) == 0)
    elseif (ep_data(var)%org_type == 0) then
        counter2 = counter2 + 1
        detritus_group_no(counter2) = var
    end if
end do

! determine the precedence of the linear equations
! to be solved; first solve equations of groups with 
! less predators rather than
! equations of groups with more predators
allocate(int_matrix(nvars))
int_matrix = group_rank(:)
do var = 1, (nvars - nproducer - ndetritus)
    ind = maxloc(int_matrix(1: (nvars - nproducer - ndetritus)), dim = 1)
    order(var) = ind
    int_matrix(ind) = -999
end do

! check to see if the precedence is correct
! be sure that biomass of predator j of prey i is not unknown
! while calculating parameters of prey i
do var = 1, (nvars - nproducer - ndetritus)
    do j = 1, dcols
        if (ep_diet(order(var), j) /= 0) then
            if (ep_data(j)%biomass == -999) then
                do n = 1, nvars
                    if (order(n) == j) then
                        m = n
                    end if
                end do
                if (m > var) then
                    where (order == j)
                        order = order(var)
                    end where
                    order(var) = j
                end if
            end if
        end if
    end do
end do

! check to align producer groups before detritus 
! so that they are calculated before detritus
lastgroups = nproducer
do var = 1, nvars
    if (ep_data(var)%org_type == 1) then
        order(nvars - ndetritus - lastgroups + 1) = var
        lastgroups = lastgroups - 1
    end if
end do

! check to align detritus groups last 
! so that they are handled last
lastgroups = ndetritus
do i = 1, ndetritus 
    order(nvars -lastgroups + 1) = detritus_group_no(i)
    lastgroups = lastgroups - 1
end do

end subroutine calculateSolvingOrderOfLinearEquations
