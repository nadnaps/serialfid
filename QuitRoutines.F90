!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : QuitRoutines.F90
!     Contains : routines for code exit
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine QuitRoutine(tin,normalexit,errorcode)
      use hdf5
      use param
      use local_arrays, only: vx,vy,vz,pr
      implicit none
      logical, intent(in) :: normalexit
      integer :: errorcode
      real(DP) :: tin(3)
      character*50 :: dsetname, filename

      if(errorcode.ne.100) then ! skip if already finalized

      !write(6,'(a,f10.2,a)') 'Total Iteration Time = ',tin(3) -tin(2),' sec.'

      call NotifyError(errorcode)
      

      if(normalexit) then
        dsetname = trim('vx') ; filename = trim('contin_vx.h5')
        call HdfCreateBlankFile(filename)
        call HdfSerialWriteReal3D(dsetname,filename,vx,n1,n2,n3)
        dsetname = trim('vy') ; filename = trim('contin_vy.h5')
        call HdfCreateBlankFile(filename)
        call HdfSerialWriteReal3D(dsetname,filename,vy,n1,n2,n3)
        dsetname = trim('vz') ; filename = trim('contin_vz.h5')
        call HdfCreateBlankFile(filename)
        call HdfSerialWriteReal3D(dsetname,filename,vz,n1,n2,n3)
        dsetname = trim('pr') ; filename = trim('contin_pr.h5')
        call HdfCreateBlankFile(filename)
        call HdfSerialWriteReal3D(dsetname,filename,pr,n1,n2,n3)
      else
       stop
      endif

      call dfftw_destroy_plan(fwd_plan)
      call dfftw_destroy_plan(bck_plan)


      call HdfClose

      call CloseLogs

      call MemoryDeAlloc

      endif

      end subroutine QuitRoutine

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine NotifyError(errorcode)
      use param
      implicit none
      integer, intent(in) :: errorcode

      if(errorcode.eq.166) then
        write(6,168) dt
168     format(10x,'dt too small, DT= ',e14.7)
      else if(errorcode.eq.165) then
        write(6,164)
164     format(10x,'cfl too large  ')
      else if(errorcode.eq.266) then
        write(6,268)
268     format(10x,'velocities diverged')
      else if(errorcode.eq.169) then
        write(6,178)
        write(6,179)
        write(6,180)
178     format(10x,'too large local residue for mass conservation at:')
179     format(10x,'Probably the matrix in SolvePressureCorrection')
180     format(10x,'is singular. Try changing nxm or str3')
        call LocateLargeDivergence
      else if(errorcode.eq.333) then
         write(*,*) "time greater than tmax"
         write(*,*) "statistics and continuation updated"
      else if(errorcode.eq.334) then
         write(*,*) "walltime greater than walltimemax"
         write(*,*) "statistics and continuation updated"
      else if(errorcode.eq.444) then
         write(*,*) "FFT size in ny or nz is not efficient"
      else
         write(*,*) "Maximum number of timesteps reached"
      end if

      return
      end

!     -----------------------------------------------------
!     -----------------------------------------------------
