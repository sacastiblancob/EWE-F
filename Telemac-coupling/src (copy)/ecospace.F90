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

subroutine ecospace (time, biomass, lat, lon)

!#ifdef isWithBFM
!  use global_mem
!#endif

!#ifndef isWithBFM
  use statevartypesecopath, only: RLEN, ep_data
!#endif
  use statevartypesecosim, only: nvars, es_data
  use statevartypesecospace, only: nlat, nlon, grid, QperB, M2, BB_spatial, advection

  implicit none

  real(RLEN), intent(in)    :: time
  real(RLEN), intent(inout) :: biomass(nvars)
  integer, intent(in)       :: lat, lon

  !in-subroutine variables
  integer    :: var
  real(RLEN) :: risk_ratio_in_grid(nvars), risk_ratio_left(nvars)
  real(RLEN) :: risk_ratio_right(nvars), risk_ratio_below(nvars)
  real(RLEN) :: risk_ratio_above(nvars)

  risk_ratio_in_grid = 0
  risk_ratio_below = 0
  risk_ratio_above = 0
  risk_ratio_left = 0
  risk_ratio_right = 0

  ! First calculate risk ratios in cell and its adjacent cells (no diagnals!)
  if (time == 0) then
      do var = 1, nvars
          if (es_data(var)%isAdvected == 0) then

!!$              risk_ratio_in_grid(var) = M2(lat, lon, var) / QperB(lat, lon, var)
              risk_ratio_in_grid(var) = bb_spatial(lat, lon, var) &
                   * ep_data(var)%eatenof / ep_data(var)%eatenby

              if (grid(lat-1, lon) == 1) then
                  risk_ratio_below(var) = BB_spatial(lat-1, lon, var) &
                       * (ep_data(var)%EatenOf / ep_data(var)%EatenBy)
              else
                  risk_ratio_below(var) = -999
              end if
              if (grid(lat+1, lon) == 1) then
                  risk_ratio_above(var) = BB_spatial(lat+1, lon, var) &
                       * (ep_data(var)%EatenOf / ep_data(var)%EatenBy)
              else
                  risk_ratio_above(var) = -999
              end if
              if (grid(lat, lon-1) == 1) then
                  risk_ratio_left(var) = BB_spatial(lat, lon-1, var) &
                       * (ep_data(var)%EatenOf / ep_data(var)%EatenBy)
              else
                  risk_ratio_left(var) = -999
              end if
              if (grid(lat, lon+1) == 1) then
                  risk_ratio_right(var) = BB_spatial(lat, lon+1, var) &
                       * (ep_data(var)%EatenOf / ep_data(var)%EatenBy)
              else
                  risk_ratio_right(var) = -999
              end if
          end if
      end do
  else
      do var = 1, nvars
          if (es_data(var)%isAdvected == 0) then

              risk_ratio_in_grid(var) = M2(lat, lon, var) &
                   / QperB(lat, lon, var)

              if (grid(lat-1, lon) == 1) then
                  risk_ratio_below(var) = M2(lat-1, lon, var) &
                       / QperB(lat-1, lon, var)
              else
                  risk_ratio_below(var) = -999
              end if
              if (grid(lat+1, lon) == 1) then
                  risk_ratio_above(var) = M2(lat+1, lon, var) &
                       / QperB(lat+1, lon, var)
              else
                  risk_ratio_above(var) = -999
              end if
              if (grid(lat, lon-1) == 1) then
                  risk_ratio_left(var) = M2(lat, lon-1, var) &
                       / QperB(lat, lon-1, var)
              else
                  risk_ratio_left(var) = -999
              end if
              if (grid(lat, lon+1) == 1) then
                  risk_ratio_right(var) = M2(lat, lon+1, var) &
                       / QperB(lat, lon+1, var)
              else
                  risk_ratio_right(var) = -999
              end if
          end if
      end do
  end if

  ! now calculate and apply biomass fluxes
  do var = 1, nvars

      if (es_data(var)%isAdvected == 0) then

          if (risk_ratio_below(var) > 0 .and. risk_ratio_in_grid(var) > risk_ratio_below(var)) then
              BB_spatial(lat-1, lon, var) = BB_spatial(lat-1, lon, var) &
                   + (1 - (risk_ratio_below(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_below(var)))) * biomass(var) * 0.1
              biomass(var) = biomass(var) - (1 - (risk_ratio_below(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_below(var)))) * biomass(var) * 0.1
          end if
          if (risk_ratio_above(var) > 0 .and. risk_ratio_in_grid(var) > risk_ratio_above(var)) then
              BB_spatial(lat+1, lon, var) = BB_spatial(lat+1, lon, var) &
                   + (1 - (risk_ratio_above(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_above(var)))) * biomass(var) * 0.1
              biomass(var) = biomass(var) - (1 - (risk_ratio_above(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_above(var)))) * biomass(var) * 0.1
          end if
          if (risk_ratio_left(var) > 0 .and. risk_ratio_in_grid(var) > risk_ratio_left(var)) then
              BB_spatial(lat, lon-1, var) = BB_spatial(lat, lon-1, var) &
                   + (1 - (risk_ratio_left(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_left(var)))) * biomass(var) * 0.1
              biomass(var) = biomass(var) - (1 - (risk_ratio_left(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_left(var)))) * biomass(var) * 0.1
          end if
          if (risk_ratio_right(var) > 0 .and. risk_ratio_in_grid(var) > risk_ratio_right(var)) then
              BB_spatial(lat, lon+1, var) = BB_spatial(lat, lon+1, var) &
                   + (1 - (risk_ratio_right(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_right(var)))) * biomass(var) * 0.1
              biomass(var) = biomass(var) - (1 - (risk_ratio_right(var) / (risk_ratio_in_grid(var) &
                   + risk_ratio_right(var)))) * biomass(var) * 0.1
          end if

      end if

  end do


  ! now calculate and apply biomass fluxes due to advection
  do var = 1, nvars

      if (es_data(var)%isAdvected == 1) then

          if (advection(lat, lon) == 180) then
              ! Advect eastward
              BB_spatial(lat, lon+1, var) = BB_spatial(lat, lon+1, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == -180) then
              ! Advect westward
              BB_spatial(lat, lon-1, var) = BB_spatial(lat, lon-1, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == -90) then
              ! Advect southward
              BB_spatial(lat-1, lon, var) = BB_spatial(lat-1, lon, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == 90) then
              ! Advect northward
              BB_spatial(lat+1, lon, var) = BB_spatial(lat+1, lon, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == -45) then
              ! Advect southeastward
              BB_spatial(lat-1, lon+1, var) = BB_spatial(lat-1, lon+1, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == 45) then
              ! Advect northeastward
              BB_spatial(lat+1, lon+1, var) = BB_spatial(lat+1, lon+1, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == -135) then
              ! Advect southwestward
              BB_spatial(lat-1, lon-1, var) = BB_spatial(lat-1, lon-1, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if
          if (advection(lat, lon) == 135) then
              ! Advect northwestward
              BB_spatial(lat+1, lon-1, var) = BB_spatial(lat+1, lon-1, var) + biomass(var) * 0.1
              biomass(var) = biomass(var) - biomass(var) * 0.1
          end if

      end if

  end do


end subroutine ecospace
