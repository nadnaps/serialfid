!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : MemRoutines.F90
!     Contains : Routines for allocating and deallocating
!                memory
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine MemoryAlloc
      use param, only: n1,n2,n3
      use local_arrays
      use aux_routines
      implicit none

      call AllocateReal3DArray(vx,1,n1,1,n2,1,n3)
      call AllocateReal3DArray(vy,1,n1,1,n2,1,n3)
      call AllocateReal3DArray(vz,1,n1,1,n2,1,n3)
      call AllocateReal3DArray(pr,1,n1,1,n2,1,n3)

      call AllocateReal3DArray(dph,1,n1,1,n2+2,1,n3)

      call AllocateReal3DArray(rhs,1,n1,1,n2, 1,n3)
      call AllocateReal3DArray(dq,1,n1,1,n2,  1,n3)
      call AllocateReal3DArray(qcap,1,n1,1,n2,1,n3)
      call AllocateReal3DArray(ru1,1,n1,1,n2, 1,n3)
      call AllocateReal3DArray(ru2,1,n1,1,n2, 1,n3)
      call AllocateReal3DArray(ru3,1,n1,1,n2, 1,n3)


      return
      end subroutine MemoryAlloc

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine MemoryDeAlloc
      use local_arrays
      implicit none

      if(allocated(vx)) deallocate(vx)
      if(allocated(vy)) deallocate(vy)
      if(allocated(vz)) deallocate(vz)
      if(allocated(pr)) deallocate(pr)

      if(allocated(dph)) deallocate(dph)

      if(allocated(rhs)) deallocate(rhs)
      if(allocated(dq)) deallocate(dq)
      if(allocated(qcap)) deallocate(qcap)

      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)


      return
      end subroutine MemoryDeAlloc

!     -----------------------------------------------------
!     -----------------------------------------------------


