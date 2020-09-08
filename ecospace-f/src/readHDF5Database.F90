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

module readHDF5Database

contains

subroutine readHDF5 (nstanzas)

use statevartypesecopath, only: HDF5_fname, ep_data, ms_data, &
                          ep_diet, ep_detfate
use hdf5

implicit none

! Ecopath basic parameters dataset
character(len = 7), parameter   :: ep_setname = "ep_data"
! Ecopath multistanza parameter dataset
character(len = 7), parameter   :: ms_setname = "ms_data"
! Ecopath diet composition matrix
character(len = 7), parameter   :: dc_setname = "ep_diet"
! Ecopath detritus fate matrix
character(len = 10), parameter  :: df_setname = "ep_detfate"

integer, intent(in)              :: nstanzas
integer(hsize_t), dimension(1:2) :: adims    ! dimensions of matrices

integer        :: hdferr
integer(hid_t) :: file_id, dset_id, memtype, filetype       ! handles
type(c_ptr)    :: f_ptr_ep, f_ptr_ms, f_ptr_dc, f_ptr_df

! initialize the Fortran interface
 call h5open_f(hdferr)

! open file
 call h5fopen_f(HDF5_fname, H5F_ACC_RDONLY_F, file_id, hdferr)

!!!!!!!!!! READ ECOPATH BASIC PARAMETERS DATASET !!!!!!!!!!!!!!!!!!!!!!!!
! open dataset
 call h5dopen_f(file_id, ep_setname, dset_id, hdferr)

! get the datatype and dimensions
 call h5dget_type_f(dset_id, memtype, hdferr)

! read the data
f_ptr_ep = c_loc(ep_data(1))
 call h5dread_f(dset_id, memtype, f_ptr_ep, hdferr)

! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5tclose_f(memtype, hdferr)

!!!!!!!!!! READ ECOPATH MULTISTANZA PARAMETER DATASET !!!!!!!!!!!!!!!!!!!
if ( nstanzas /= 0 ) then

    ! open dataset
    call h5dopen_f(file_id, ms_setname, dset_id, hdferr)

    ! get the datatype and dimensions
    call h5dget_type_f(dset_id, memtype, hdferr)

    ! read the data
    f_ptr_ms = c_loc(ms_data(1))
    call h5dread_f(dset_id, memtype, f_ptr_ms, hdferr)

    ! close and release resources
    call h5dclose_f(dset_id, hdferr)
    call h5tclose_f(memtype, hdferr)
end if

!!!!!!!!!! READ DIET COMPOSITION MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataset
 call h5dopen_f(file_id, dc_setname, dset_id, hdferr)

! get the datatype and dimensions
 call h5dget_type_f(dset_id, filetype, hdferr)
 call h5tget_array_dims_f(filetype, adims, hdferr)

! allocate memory for read buffer.
allocate(ep_diet(adims(2), adims(1)))

! read the data
f_ptr_dc = c_loc(ep_diet)
 call h5dread_f(dset_id, filetype, f_ptr_dc, hdferr)

! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5tclose_f(filetype, hdferr)
 
!!!!!!!!!! READ DETRITUS FATE MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataset
 call h5dopen_f(file_id, df_setname, dset_id, hdferr)

! get the datatype and dimensions
 call h5dget_type_f(dset_id, filetype, hdferr)
 call h5tget_array_dims_f(filetype, adims, hdferr)

! allocate memory for read buffer.
allocate(ep_detfate(adims(2), adims(1)))

! read the data
f_ptr_df = c_loc(ep_detfate)
 call h5dread_f(dset_id, filetype, f_ptr_df, hdferr)

! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5tclose_f(filetype, hdferr)
 call h5fclose_f(file_id, hdferr)

 call h5close_f(hdferr)

end subroutine readHDF5

end module readHDF5Database
