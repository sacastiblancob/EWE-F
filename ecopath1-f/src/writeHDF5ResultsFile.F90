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

subroutine writeHDF5ResultsFile (drows, dcols, detrows, detcols)

use statevartypesecopath, only: Results_fname, ep_data, ms_data, &
                          ep_diet, ep_detfate, nvars, nstanzas
use hdf5

implicit none

! This is the name of the data file we will read. 
! name of the Ecopath base dataset
character (len = 7), parameter  :: ep_setname = "ep_data"
! name of the Ecopath multistanza dataset
character (len = 7), parameter  :: ms_setname = "ms_data"
! name of the Ecopath diet dataset
character (len = 7), parameter  :: dc_setname = "ep_diet"
! name of the Ecopath detritus fate dataset
character (len = 10), parameter :: df_setname = "ep_detfate"

! variables inherited from Ecopath-FORTRAN model
integer, intent(in)               :: drows, dcols, detrows, detcols

! in-subroutine variable declarations for ep_data
integer                        :: ep_dim0, ms_dim0, dim0
integer, parameter             :: rank   = 1
integer                        :: hdferr
integer(hid_t)                 :: file_id, dset_id, dspace_id  ! handles
integer(hid_t)                 :: memtype, filetype, plist_id  ! handles
integer, parameter             :: ardim0 = 10
integer, parameter             :: ardim1 = 299
integer(hsize_t), dimension(1) :: ep_dims, ms_dims, dims
integer(hsize_t), dimension(1) :: f5len, l4len
integer(hsize_t), dimension(2) :: dc_dims, df_dims
type(c_ptr)                    :: f_ptr_ep, f_ptr_ms, f_ptr_dc, f_ptr_df

integer(hid_t) :: s1_tid, s2_tid, s3_tid, s4_tid, s5_tid, s6_tid, s7_tid
integer(hid_t) :: s8_tid, s9_tid, s10_tid, s11_tid, s12_tid, s13_tid
integer(hid_t) :: s14_tid, s15_tid, s16_tid, s17_tid, s18_tid, s19_tid
integer(hid_t) :: s20_tid, s21_tid, s22_tid, s_tid
integer(hid_t) :: RHDF5

 ! initialize FORTRAN interface
 call h5open_f(hdferr)

#ifdef _isDP_
  call h5tcopy_f (H5T_NATIVE_DOUBLE, RHDF5, hdferr);
#endif
#ifndef _isDP_
  call h5tcopy_f (H5T_NATIVE_REAL, RHDF5, hdferr);
#endif

 ! transfer properties
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 call h5pset_preserve_f(plist_id, .TRUE., hdferr)

 ! create file, if it already exists overwrite (H5F_ACC_TRUNC_F)
 call h5fcreate_f(Results_fname, H5F_ACC_TRUNC_F, file_id, hdferr)

!!!!!!!!!! ECOPATH DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! dimensions of (n) ecopath data matrix
 ep_dim0 = nvars
 ep_dims = (/ ep_dim0 /)

 ! create the dataspace
 call h5screate_simple_f(rank, ep_dims, dspace_id, hdferr)
  
 ! create the compound datatype for memory
 call h5tcreate_f(H5T_COMPOUND_F,  H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(2))), memtype, hdferr)

   ! insert members
 call h5tinsert_f(memtype, "biomass", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%biomass)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "PoB", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%PoB)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "QoB", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%QoB)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "EE", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%EE)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "PoQ", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%PoQ)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "unass_Q", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%unass_Q)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "detritus_import", &
                H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%detritus_import)), &
                RHDF5, hdferr)
 call h5tinsert_f(memtype, "landings", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%landings)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "discards", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%discards)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "org_type", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%org_type)), H5T_NATIVE_INTEGER, hdferr)
 call h5tinsert_f(memtype, "isstanza", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%isstanza)), H5T_NATIVE_INTEGER, hdferr)
 call h5tinsert_f(memtype, "stanza_no", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%stanza_no)), H5T_NATIVE_INTEGER, hdferr)
 call h5tinsert_f(memtype, "age_start", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%age_start)), H5T_NATIVE_INTEGER, hdferr)
 call h5tinsert_f(memtype, "isleading", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%isleading)), H5T_NATIVE_INTEGER, hdferr)
 call h5tinsert_f(memtype, "production", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%production)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "consumption", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%consumption)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "respiration", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%respiration)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "assimilation", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%assimilation)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "EatenOf", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%EatenOf)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "EatenBy", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%EatenBy)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "DetPassedProp", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%DetPassedProp)), RHDF5, hdferr)
 call h5tinsert_f(memtype, "BA", H5OFFSETOF(c_loc(ep_data(1)), &
                c_loc(ep_data(1)%BA)), RHDF5, hdferr)

 ! Create the dataset and write the array data to it.
 call h5dcreate_f(file_id, ep_setname, memtype, dspace_id, dset_id, &
                hdferr)
 
 ! write data
 f_ptr_ep = C_LOC(ep_data(1))
 call h5dwrite_f(dset_id, memtype, f_ptr_ep, hdferr, xfer_prp = plist_id)
 
 ! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5tclose_f(memtype, hdferr)
 call h5sclose_f(dspace_id, hdferr)

!!!!!!!!!! MULTISTANZA DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nstanzas /= 0) then

 ! dimension (n) of multistanza data matrix
 ms_dim0 = nstanzas
 ms_dims = (/ ms_dim0 /)
 
 ! dimesion (n) of subarrays in multistanza data matrix
 f5len = (/ ardim0 /)
 l4len = (/ ardim1 /)

 ! create the dataspace
 call h5screate_simple_f(rank, ms_dims, dspace_id, hdferr)

 ! create arrays
 call H5Tarray_create_f(H5T_NATIVE_INTEGER, rank, f5len, s1_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, f5len, s2_tid, hdferr)
 call H5Tarray_create_f(H5T_NATIVE_INTEGER, rank, f5len, s3_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, f5len, s4_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, f5len, s5_tid, hdferr)
 call H5Tarray_create_f(H5T_NATIVE_INTEGER, rank, f5len, s6_tid, hdferr)
 call H5Tcopy_f(H5T_NATIVE_INTEGER , s7_tid, hdferr)
 call H5Tcopy_f(H5T_NATIVE_INTEGER , s8_tid, hdferr)
 call H5Tcopy_f(RHDF5, s9_tid, hdferr)
 call H5Tcopy_f(RHDF5, s10_tid, hdferr)
 call H5Tcopy_f(RHDF5, s11_tid, hdferr)
 call H5Tcopy_f(RHDF5, s12_tid, hdferr)
 call H5Tcopy_f(RHDF5, s13_tid, hdferr)
 call H5Tcopy_f(RHDF5, s14_tid, hdferr)
 call H5Tcopy_f(RHDF5, s15_tid, hdferr)
 call H5Tcopy_f(RHDF5, s16_tid, hdferr)
 call H5Tcopy_f(RHDF5, s17_tid, hdferr)
 call H5Tcopy_f(RHDF5, s18_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, l4len, s19_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, l4len, s20_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, l4len, s21_tid, hdferr)
 call H5Tarray_create_f(RHDF5, rank, l4len, s22_tid, hdferr)

 ! create the compound datatype "s_tid"  which will embrace the other arrays
 call H5Tcreate_f(H5T_COMPOUND_F, sizeof(ms_data(1)), s_tid, hdferr)

 ! insert arrays into the compound datatype
 call H5Tinsert_f(s_tid, "ep_groupno", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%ep_groupno)), s1_tid, hdferr)
 call H5Tinsert_f(s_tid, "biomass", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%biomass)), s2_tid, hdferr)
 call H5Tinsert_f(s_tid, "age_start", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%age_start)), s3_tid, hdferr)
 call H5Tinsert_f(s_tid, "mortality", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%mortality)), s4_tid, hdferr)
 call H5Tinsert_f(s_tid, "QoB", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%QoB)), s5_tid, hdferr)
 call H5Tinsert_f(s_tid, "isleading", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%isleading)), s6_tid, hdferr)
 call H5Tinsert_f(s_tid, "age_infinity", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%age_infinity)), s7_tid, hdferr)
 call H5Tinsert_f(s_tid, "substanzas", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%substanzas)), s8_tid, hdferr)
 call H5Tinsert_f(s_tid, "vbK", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%vbK)), s9_tid, hdferr)
 call H5Tinsert_f(s_tid, "vbM", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%vbM)), s10_tid, hdferr)
 call H5Tinsert_f(s_tid, "rec_power", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%rec_power)), s11_tid, hdferr)
 call H5Tinsert_f(s_tid, "rel_BA", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%rel_BA)), s12_tid, hdferr)
 call H5Tinsert_f(s_tid, "Wmat_Winf", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%Wmat_Winf)), s13_tid, hdferr)
 call H5Tinsert_f(s_tid, "RzeroS", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%RzeroS)), s14_tid, hdferr)
 call H5Tinsert_f(s_tid, "Rhat", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%Rhat)), s15_tid, hdferr)
 call H5Tinsert_f(s_tid, "Ahat", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%Ahat)), s16_tid, hdferr)
 call H5Tinsert_f(s_tid, "RhatC", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%RhatC)), s17_tid, hdferr)
 call H5Tinsert_f(s_tid, "AhatC", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%AhatC)), s18_tid, hdferr)
 call H5Tinsert_f(s_tid, "Wage", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%Wage)), s19_tid, hdferr)
 call H5Tinsert_f(s_tid, "WWa", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%WWa)), s20_tid, hdferr)
 call H5Tinsert_f(s_tid, "survive", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%survive)), s21_tid, hdferr)
 call H5Tinsert_f(s_tid, "splitno", H5OFFSETOF(c_loc(ms_data(1)), &
                c_loc(ms_data(1)%splitno)), s22_tid, hdferr)

 ! Create the dataset and write the array data to it.
 call h5dcreate_f(file_id, ms_setname, s_tid, dspace_id, dset_id, hdferr)

 ! write data to file
 f_ptr_ms = c_loc(ms_data(1))
 call h5dwrite_f(dset_id, s_tid, f_ptr_ms, hdferr)

 ! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5tclose_f(s_tid, hdferr)
 call h5sclose_f(dspace_id, hdferr)

endif
!!!!!!!!!! DIET COMPOSITION MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! dimensions (n x m) of diet composition matrix
 dim0 = 1
 dims = (/ dim0 /)
 dc_dims = (/ dcols, drows /)
 
 ! create diet composition matrix datatype for file and memory
 call H5Tarray_create_f(INT(RHDF5, HID_T), 2, dc_dims, &
                      filetype, hdferr)
 call H5Tarray_create_f(RHDF5, 2, dc_dims, memtype, hdferr)

 ! create the dataspace
 call h5screate_simple_f(rank, dims, dspace_id, hdferr)

 ! create the dataset
 call h5dcreate_f(file_id, dc_setname, filetype, dspace_id, dset_id, &
                hdferr)

 ! write data
 f_ptr_dc = c_loc(ep_diet)
 call h5dwrite_f(dset_id, memtype, f_ptr_dc, hdferr)

 ! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(dspace_id, hdferr)
 call h5tclose_f(filetype, hdferr)
 call h5tclose_f(memtype, hdferr)

!!!!!!!!!! DETRITUS FATE MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! dimensions (n x m) of detritus fate matrix
 df_dims = (/ detcols, detrows /)
 
 ! create diet composition matrix datatype for file and memory
 call H5Tarray_create_f(INT(RHDF5, HID_T), 2, df_dims, &
                      filetype, hdferr)
 call H5Tarray_create_f(RHDF5, 2, df_dims, memtype, hdferr)

 ! create the dataspace
 call h5screate_simple_f(rank, dims, dspace_id, hdferr)

 ! create the dataset
 call h5dcreate_f(file_id, df_setname, filetype, dspace_id, dset_id, &
                hdferr)

 ! write data
 f_ptr_df = c_loc(ep_detfate)
 call h5dwrite_f(dset_id, memtype, f_ptr_df, hdferr)

 ! close and release resources
 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(dspace_id, hdferr)
 call h5tclose_f(filetype, hdferr)
 call h5tclose_f(memtype, hdferr)
 call h5pclose_f(plist_id, hdferr)
 call h5fclose_f(file_id, hdferr)

 ! close FORTRAN interface
 call h5close_f(hdferr)

end subroutine writeHDF5ResultsFile

