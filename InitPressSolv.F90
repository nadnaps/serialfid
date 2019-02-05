!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : InitPressSolv.F90
!     Contains : Routines for initialising pressure solver
!
!     -----------------------------------------------------
!     ----------------------------------------------------- 

!   Initialise pressure solver and fftw plans

      subroutine InitPressSolv
      use param
      implicit none
      integer FFTW_EXHAUSTIVE
      parameter(FFTW_EXHAUSTIVE=64)
      real(DP), allocatable, dimension(:,:) :: xr
      complex(DP), allocatable, dimension(:,:) :: xa

      allocate(xr(n2m,n1m))
      allocate(xa(n2mh,n1m))
    
      call InitFourierMetric

      call TridiagMatrices

      call dfftw_plan_dft_r2c_2d(fwd_plan,n2m,n1m,xr,xa,FFTW_EXHAUSTIVE)
      call dfftw_plan_dft_c2r_2d(bck_plan,n2m,n1m,xa,xr,FFTW_EXHAUSTIVE)

      if(allocated(xr)) deallocate(xr)
      if(allocated(xa)) deallocate(xa)


      return
      end subroutine InitPressSolv
      
!     -----------------------------------------------------
!     -----------------------------------------------------

!     Perform the calculation of trigz for temperton fft

      subroutine InitFourierMetric
      use param
      implicit none
      integer :: n2mp,j,i,n1mp
      integer :: k, n3mh, n3mp


      n1mp=n1mh+1
      n2mp=n2mh+1

      allocate(ao(1:n1))
      allocate(ap(1:n2))

      allocate(ak1(1:n1))
      allocate(ak2(1:n2))
!
!     wave number definition
!
      do i=1,n1mh
        ao(i)=(i-1)*2.d0*pi
      enddo
      do i=n1mp,n1m
        ao(i)=-(n1m-i+1)*2.d0*pi
      enddo
      do i=1,n1m
        ak1(i)=2.d0*(1.d0-dcos(ao(i)/n1m))*(float(n1m)/xlen)**2
      enddo

      do j=1,n2mh
        ap(j)=(j-1)*2.d0*pi
      enddo
      do j=n2mp,n2m
        ap(j)=-(n2m-j+1)*2.d0*pi
      enddo
      do j=1,n2m
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/n2m))*(float(n2m)/ylen)**2
      enddo

      return
      end
      
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine TridiagMatrices
      use param
      implicit none
      integer  :: kc,km,kp
      real(DP) :: ugmmm,a33icc,a33icp

      allocate(amphk(1:n3))
      allocate(acphk(1:n3))
      allocate(apphk(1:n3))


!   tridiagonal matrix coefficients at each k and i
!   x1 and x3 cartesian coordinates

      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo

      end subroutine TridiagMatrices

!     -----------------------------------------------------
!     -----------------------------------------------------
