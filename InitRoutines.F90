!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : InitRoutines.F90
!     Contains : Routines for initialising different things
!
!     -----------------------------------------------------
!     -----------------------------------------------------


      subroutine InitTimeMarchScheme
      use param
      implicit none
      
      integer :: ns

      if(nsst.gt.1) then
      
        gam(1) = 8.d0/15.d0
        gam(2) = 5.d0/12.d0
        gam(3) = 3.d0/4.d0

        rom(1) = 0.d0
        rom(2) =-17.d0/60.d0
        rom(3) =-5.d0/12.d0

!     output to screen      
      write(6,121)(gam(ns),ns=1,nsst),(rom(ns),ns=1,nsst)
  121 format(/,10x,'Time Scheme is III order Runge-Kutta', &
             4x,'gam= ',3f8.3,4x, 'ro= ',3f8.3)
      
      else

        gam(1) = 1.5d0 
        gam(2) = 0.d0 
        gam(3) = 0.d0 

        rom(1) = -0.5d0
        rom(2) = 0.d0
        rom(3) = 0.d0

!     output to screen      
      write(6,122)gam(1),rom(1)
  122 format(/,10x,'Time Scheme is Adams-Bashfort', &
             4x,'gam= ',f8.3,4x, 'ro= ',f8.3)

      end if

      do ns=1,nsst
        alm(ns) = gam(ns)+rom(ns)
      end do

      return
      end subroutine InitTimeMarchScheme

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine InitArrays
      use param
      use local_arrays
      implicit none
      integer :: j,k,i
      
      do k=1,n3
      do j=1,n2
      do i=1,n1
      
      vx(i,j,k)=0.d0
      vy(i,j,k)=0.d0
      vz(i,j,k)=0.d0

      pr(i,j,k)=0.d0
      dph(i,j,k)=0.d0
      dq(i,j,k)=0.d0
      rhs(i,j,k)=0.d0
      ru1(i,j,k)=0.d0
      ru2(i,j,k)=0.d0
      ru3(i,j,k)=0.d0
      qcap(i,j,k)=0.d0
      
      
      enddo
      enddo
      enddo

      return 
      end subroutine InitArrays

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine CreateInitCond
      use param
      use local_arrays, only: vx,vy,vz
      implicit none
      integer :: j,k,i
      
      do k=1,n3
      do j=1,n2
      do i=1,n1
      
      vx(i,j,k)=0.d0
      vy(i,j,k)=0.d0
      vz(i,j,k)=0.d0
      
      enddo
      enddo
      enddo

      return 
      end subroutine CreateInitCond

!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine ReadInitCond
      use param
      use local_arrays, only: vx,vy,vz
      IMPLICIT NONE
      character*70 :: filcn
      integer :: n1o,n2o,n3o   
      
!    Reading old grid information
      filcn = 'contin_grid.dat'
      open(13,file=filcn,status='unknown')
      rewind(13)                                                      
      read(13,*) n1o,n2o,n3o,time
      close(13)
      
      if(n1o.ne.n1.or.n2o.ne.n2.or.n3o.ne.n3) then
        write(6,*)'Grids dont match in ReadInitCond ',n1,n2,n3,' <-- NE --> ',n1o,n2o,n3o
      end if
        
!    One to one HDF read
      call HdfSerialReadReal3D('vx','contin_vx.h5',vx,n1,n2,n3)
      call HdfSerialReadReal3D('vy','contin_vy.h5',vy,n1,n2,n3)
      call HdfSerialReadReal3D('vz','contin_vz.h5',vz,n1,n2,n3)

       vx(1:n1,1:n2,n3) = 0.0
       vy(1:n1,1:n2,n3) = 0.0
       vz(1:n1,1:n2,n3) = 0.0

      if (ireset.eq.1) then                                             
       time=0.
      endif                                                             

      return                                    
      end                                                               

