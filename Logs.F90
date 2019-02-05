!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : Logs.F90
!     Contains : Routines for opening and closing simple 
!                output files
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine OpenLogs
      use param, only : ireset
      implicit none

      open(210,file='force.out',status='unknown',position='append',access='sequential')


      if(ireset.eq.1)then
            rewind(210)
      end if

      return
      end subroutine OpenLogs

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine CloseLogs
      implicit none

      close(210)

      return
      end subroutine CloseLogs

!     -----------------------------------------------------
!     -----------------------------------------------------


