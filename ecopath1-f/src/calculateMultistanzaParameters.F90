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

subroutine calculateMultistanzaParameters ()

! This subroutine is run if there are multistanza groups in the model
! The subroutine is called from the ecopath.f90 and run 
! before mass-balance calculations

use statevartypesecopath, only: RLEN
use statevartypesecopath, only: GrowthParam_fname, isASCIIinputFile, &
                          ep_data, ep_growth, ms_data, nvars, nstanzas

implicit none

integer :: counter, i, j        ! loop variables
integer :: paramrows, paramcols ! # of rows and columns of growth data
integer :: stanza               ! loop var. # of multistanza groups
real(RLEN) :: previous_survival    ! survival rate of previous age
real(RLEN) :: survival             ! monthly survival rate of stanza

! sum of biomasses of all ages in the leading stanza
real(RLEN) :: SumB

real(RLEN) :: recruits             ! # of recruits to the leading stanza
real(RLEN) :: Bio_t                ! the biomass of a substanza
real(RLEN) :: K                    ! the real K value from VBGF

! temporary matrix to read growth parameters from file
real(RLEN), allocatable :: m_data(:, :)

! accumulated biomass in each substanza in the multistanza group
real(RLEN), allocatable :: BA_t(:)

! location of calculated stanza parameters in the main ep_data matrix
integer, allocatable :: ep_dataloc(:)

real(RLEN), parameter   :: Fill_Value0 = -999.0
real(RLEN), parameter   :: Fill_Value1 = -999.0
integer, parameter   :: Fill_Value2 = -999
character ( len = 250 ) :: groupname, dummy

!!!!!!!!!! SECTION: READING GROWTH PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!

! open growth parameter data file
open(999, file = GrowthParam_fname, form = 'formatted')

! read number of rows and columns in the file
read(999, *) dummy, paramrows
read(999, *) dummy, paramcols
read(999, *)

! allocate size of multistanza data in variable <m_data>
allocate(m_data(paramrows, paramcols))

! read growth parameters (vbK etc.) input data from file
if (isASCIIinputFile) then
    do i = 1, paramrows
        read(999, *, end = 900) groupname, m_data(i, :)
    end do
end if

900 continue
 close(999)

do i = 1, nstanzas
    if (isASCIIinputFile) then
        ms_data(i)%vbK = m_data(i, 1)
        ms_data(i)%rec_power = m_data(i, 2)
        ms_data(i)%rel_BA  = m_data(i, 3)
        ms_data(i)%Wmat_Winf = m_data(i, 4)
    else
        ms_data(ep_growth(i)%stanza_no)%vbK = ep_growth(i)%K
        ms_data(ep_growth(i)%stanza_no)%rec_power = ep_growth(i)%RecPow
        ms_data(ep_growth(i)%stanza_no)%rel_BA  = ep_growth(i)%BaB
        ms_data(ep_growth(i)%stanza_no)%Wmat_Winf = ep_growth(i)%WmatWinf
    end if
end do

!!!!!!!!!! END SECTION: READING GROWTH PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! SECTION: CALCULATION OF STANZA PARAMETERS !!!!!!!!!!!!!!!!!!!!

do stanza = 1, nstanzas

    ! Calculate the number of substanzas in the multistanza group
    counter = 0
    do i = 1, nvars
        if (ep_data(i)%stanza_no == stanza) then
            counter = counter + 1
        end if
    end do
    ms_data(stanza)%substanzas = counter

    ! allocate the original places of transfered variables in the ep_data
    allocate(ep_dataloc(counter))

    ! fill the rest of the first five arrays of the ms_data matrix 
    ! with respective Fill_Value
    do i = (counter + 1), 10
        ms_data(stanza)%ep_groupno = Fill_Value2
        ms_data(stanza)%biomass(i) = Fill_Value0
        ms_data(stanza)%age_start(i) = Fill_Value2
        ms_data(stanza)%mortality(i) = Fill_Value0
        ms_data(stanza)%QoB(i) = Fill_Value0
        ms_data(stanza)%isleading(i) = Fill_Value2
    end do

    ! then transfer data to in-subroutine variables
    j = 0
    do i = 1, nvars
        if (ep_data(i)%stanza_no == stanza) then
            j = j + 1
            ms_data(stanza)%ep_groupno(j)= i
            ms_data(stanza)%biomass(j) = ep_data(i)%biomass
            ms_data(stanza)%age_start(j) = ep_data(i)%age_start
            ms_data(stanza)%mortality(j) = ep_data(i)%pob
            ms_data(stanza)%QoB(j) = ep_data(i)%QoB
            ms_data(stanza)%isleading(j) = ep_data(i)%isleading
        end if
    end do

    ! now it is time to move forward for calculations

    ms_data(stanza)%age_infinity = nint(40.3978 / ms_data(stanza)%vbK)

    do i = (ms_data(stanza)%age_infinity + 1), 299
        ms_data(stanza)%Wage(i) = Fill_Value1
        ms_data(stanza)%WWa(i) = Fill_Value1
    end do

    ms_data(stanza)%vbM = 1 - 3 * ms_data(stanza)%vbK / 12

    ! calculate weight at age (Wage, monthly) and 
    ! the relative Q (WWa, monthly) depending on this weight
    do i = 0, ms_data(stanza)%age_infinity
        ! Weight based von Bertalanffy equation
        ms_data(stanza)%Wage(i) = real((1 - exp(dble(-ms_data(stanza)%vbK) &
          * i / 12.0D0)) ** 3.0D0, 4)
        ms_data(stanza)%WWa(i) = real(dble(ms_data(stanza)%Wage(i)) ** (2.0D0 / 3.0D0), 4)
    end do

    ! calculate survival rate of each monthly cohort 
    ! in the multistanza group
    do i = (ms_data(stanza)%age_infinity + 1), 299
        ms_data(stanza)%survive(i) = Fill_Value1
    end do

    ms_data(stanza)%survive(0) = 1
    previous_survival = 1

    do i = 1, counter
        survival = real(exp(-(dble(ms_data(stanza)%mortality(i)) &
          + dble(ms_data(stanza)%rel_BA)) / 12.0D0), 4)
  
        if (survival > 0) then
    
            if (ms_data(stanza)%age_start(i) > 0) then
                ms_data(stanza)%survive(ms_data(stanza)%age_start(i)) &
                  = ms_data(stanza)%survive(ms_data(stanza)%age_start(i) &
                  - 1) * previous_survival
            end if
    
            if (i < counter) then
                do j = (ms_data(stanza)%age_start(i) + 1), &
                  ms_data(stanza)%age_start(i + 1)
                    ms_data(stanza)%survive(j) &
                      = ms_data(stanza)%survive(j - 1) * survival
                end do
            else
                do j = (ms_data(stanza)%age_start(i) + 1), &
                  ms_data(stanza)%age_infinity
                    ms_data(stanza)%survive(j) &
                      = ms_data(stanza)%survive(j - 1) * survival
                end do
            end if
            previous_survival = survival
        end if
    end do

    !! now it is time for accumulator calculations for last age
    if (survival < 1) then
        ms_data(stanza)%survive(j - 1) = ms_data(stanza)%survive(j - 1) &
          / (1.0 - survival)
    end if

    ms_data(stanza)%Rhat &
      = (ms_data(stanza)%Wage(ms_data(stanza)%age_infinity) &
      - ms_data(stanza)%Wage(ms_data(stanza)%age_infinity - 1)) &
      / (ms_data(stanza)%Wage(ms_data(stanza)%age_infinity - 1) &
      - ms_data(stanza)%Wage(ms_data(stanza)%age_infinity - 2))

    ms_data(stanza)%Ahat = ms_data(stanza)%Wage(ms_data(stanza)%age_infinity) &
      - ms_data(stanza)%Rhat &
      * ms_data(stanza)%Wage(ms_data(stanza)%age_infinity - 1)

    ms_data(stanza)%Wage(ms_data(stanza)%age_infinity) = (survival &
      * ms_data(stanza)%Ahat + (1 - survival) &
      * ms_data(stanza)%Wage(ms_data(stanza)%age_infinity)) &
      / (1 - ms_data(stanza)%Rhat * survival)

    ms_data(stanza)%RhatC &
      = (ms_data(stanza)%WWa(ms_data(stanza)%age_infinity) &
      - ms_data(stanza)%WWa(ms_data(stanza)%age_infinity - 1)) &
      / (ms_data(stanza)%WWa(ms_data(stanza)%age_infinity - 1) &
      - ms_data(stanza)%WWa(ms_data(stanza)%age_infinity - 2))

    ms_data(stanza)%AhatC &
      = ms_data(stanza)%WWa(ms_data(stanza)%age_infinity) &
      - ms_data(stanza)%Rhat &
      * ms_data(stanza)%WWa(ms_data(stanza)%age_infinity - 1)

    ms_data(stanza)%WWa(ms_data(stanza)%age_infinity) = (survival &
      * ms_data(stanza)%AhatC + (1 - survival) &
      * ms_data(stanza)%WWa(ms_data(stanza)%age_infinity)) &
      / (1 - ms_data(stanza)%RhatC * survival)
    !! end of accumulator calculations

    ! calculate the biomass of the all ages in the leading stanza
    SumB = 0
    do i = ms_data(stanza)%age_start(counter), &
      ms_data(stanza)%age_infinity
        SumB = SumB + ms_data(stanza)%survive(i) &
          * ms_data(stanza)%Wage(i)
    end do

    ! number of recruits from the leading stanza
    recruits = ms_data(stanza)%biomass(counter) / SumB
    ms_data(stanza)%RzeroS   = real(dble(recruits) * exp(dble(ms_data(stanza)%rel_BA) &
      / 12.0D0), 4)

    ! calculate the survivors from each monthly cohort
    do i = (ms_data(stanza)%age_infinity + 1), 299
        ms_data(stanza)%splitno(i) = Fill_Value1
    end do

    do j = 0, ms_data(stanza)%age_infinity
        ms_data(stanza)%splitno(j) = recruits &
          * ms_data(stanza)%survive(j)
    end do

    ! calculate the total biomass of each substanza 
    ! in the multistanza group
    do i = 1, (counter - 1)
        Bio_t = 0
        do j = ms_data(stanza)%age_start(i), &
          (ms_data(stanza)%age_start(i + 1) - 1)
            Bio_t = Bio_t + ms_data(stanza)%splitno(j) &
              * ms_data(stanza)%Wage(j)
        end do
        ms_data(stanza)%biomass(i) = Bio_t
    end do

    !! now calculate Q/B values of stanzas
    ! this is the real K value
    K = 0
    do i = ms_data(stanza)%age_start(counter), &
      ms_data(stanza)%age_infinity
        K = K + ms_data(stanza)%splitno(i) * ms_data(stanza)%WWa(i)
    end do

    if (K > 0) then
        K = ms_data(stanza)%QoB(counter) &
          * ms_data(stanza)%biomass(counter) / K
    end if

    ! accumulated biomass in each substanza in the multistanza group
    allocate(BA_t(counter))
    do i = 1, (counter - 1)

        BA_t(i) = ms_data(stanza)%rel_BA * ms_data(stanza)%biomass(i)
        ms_data(stanza)%QoB(i) = 0
  
        do j = ms_data(stanza)%age_start(i), &
          (ms_data(stanza)%age_start(i + 1) - 1)
            ms_data(stanza)%QoB(i) = ms_data(stanza)%QoB(i) &
            + ms_data(stanza)%splitno(j) * ms_data(stanza)%WWa(j)
        end do

        ms_data(stanza)%QoB(i) = K * ms_data(stanza)%QoB(i)
        if (ms_data(stanza)%biomass(i) > 0) then
            ms_data(stanza)%QoB(i) = ms_data(stanza)%QoB(i) &
            / ms_data(stanza)%biomass(i)
        end if

    end do
    !! end of Q/B calculations

    ! transfer calculated stanza parameters to the ep_data matrix
    do i = 1, counter
        ep_data(ms_data(stanza)%ep_groupno(i))%biomass &
          = ms_data(stanza)%biomass(i)
        ep_data(ms_data(stanza)%ep_groupno(i))%QoB &
          = ms_data(stanza)%QoB(i)
    end do

    ! deallocate variables for use in the next iteration
    deallocate(ep_dataloc)
    deallocate(BA_t)

end do

end subroutine calculateMultistanzaParameters
