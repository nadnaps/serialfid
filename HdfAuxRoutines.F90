!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : HdfAuxRoutines.F90
!     Contains : Routines for h5 files open and close
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfStart
      use hdf5
      implicit none
      integer :: hdf_error

      call h5open_f(hdf_error)

      end subroutine HdfStart

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfClose
      use hdf5
      implicit none
      integer :: hdf_error

      call h5close_f(hdf_error)

      end subroutine HdfClose

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfSerialWriteReal1D(dsetname,filename,var,sz)
      use constants
      use hdf5
      implicit none
      character*50,intent(in) :: dsetname,filename
      integer, intent(in) :: sz
      real(DP), dimension(1:sz), intent(in) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)
      logical :: dsetexists

      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sz

      call h5screate_simple_f(1, dims, filespace, hdf_error)

      call h5lexists_f(file_id,dsetname,dsetexists,hdf_error)

      if(dsetexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
                      filespace, dset, hdf_error)


       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
         var(1:sz), dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfSerialWriteReal1D

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfSerialWriteReal2D(dsetname,filename,var,sx,sy)
      use constants
      use hdf5
      implicit none
      character*50,intent(in) :: dsetname,filename
      integer, intent(in) :: sx,sy
      real(DP), dimension(sx,sy), intent(in) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error,ndims
      integer(HSIZE_T) :: dims(2)
      logical :: dsetexists


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      ndims = 2
      dims(1)=sx
      dims(2)=sy

      call h5screate_simple_f(ndims, dims, &
                              filespace, hdf_error)

      call h5lexists_f(file_id,dsetname,dsetexists,hdf_error)

      if(dsetexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
                      filespace, dset, hdf_error)


       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
         var(1:sx,1:sy), dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfSerialWriteReal2D

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfSerialWriteReal3D(dsetname,filename,var,sx,sy,sz)
      use constants
      use hdf5
      implicit none
      character*50,intent(in) :: dsetname,filename
      integer, intent(in) :: sx,sy,sz
      real(DP), dimension(sx,sy,sz), intent(in) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset, filespace
      integer :: hdf_error,ndims
      integer(HSIZE_T) :: dims(3)
      logical :: dsetexists


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      ndims = 3
      dims(1)=sx
      dims(2)=sy
      dims(3)=sz

      call h5screate_simple_f(ndims, dims, &
                              filespace, hdf_error)

      call h5lexists_f(file_id,dsetname,dsetexists,hdf_error)

      if(dsetexists) call h5ldelete_f(file_id,dsetname,hdf_error)

      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
                      filespace, dset, hdf_error)


       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
         var(1:sx,1:sy,1:sz), dims, hdf_error)


      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfSerialWriteReal3D

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfCreateBlankFile(filename)
      use constants
      use hdf5
      implicit none
      character*50,intent(in) :: filename
      integer(HID_T) :: file_id
      integer :: hdf_error

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfCreateBlankFile

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfSerialReadReal1D(dsetname,filename,var,sz)
      use constants
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sz
      real(DP), dimension(sz), intent(out) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sz

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
         var, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfSerialReadReal1D

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfSerialReadReal2D(dsetname,filename,var,sx,sy)
      use constants
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sx,sy
      real(DP), dimension(sx,sy), intent(out) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(2)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sx
      dims(2)=sy

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
         var, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfSerialReadReal2D

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine HdfSerialReadReal3D(dsetname,filename,var,sx,sy,sz)
      use constants
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: sx,sy,sz
      real(DP), dimension(sx,sy,sz), intent(out) :: var
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(3)


      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

      dims(1)=sx
      dims(2)=sy
      dims(3)=sz

      call h5dopen_f(file_id, dsetname, dset, hdf_error)

      call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
         var, dims, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfSerialReadReal3D

!     -----------------------------------------------------
!     -----------------------------------------------------
