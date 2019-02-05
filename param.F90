!     ===========================================================
!     Declaration of global variables
!     ===========================================================

!     -----------------------------------------------------
      module constants
        integer, parameter :: DP = kind(1.d0)
      end module

!     -----------------------------------------------------
      module param
        use constants
        implicit none

!     Read from input file
        integer :: n1m,n2m,n3m
        integer :: nsst,nwrit,nread
        integer :: ntst,ireset
        integer :: idtv
        integer :: starea
        integer :: istr3,str3
        integer :: ubcbot, ubctop

        real(DP):: tframe,tpin,tmax,walltimemax
        real(DP):: xlen,ylen,zlen 
        real(DP):: ren
        real(DP):: dt,dtmax,cflmax,cfllim,resid
        real(DP)::tsta


!     Other variables

!     Grid parameters
        integer  :: n1,n2,n3,n1mh,n2mh
        real(DP) :: dx1,dx2,dx3
        real(DP) :: dx1q,dx2q,dx3q

        real(DP),allocatable,dimension(:) :: xc,xm
        real(DP),allocatable,dimension(:) :: yc,ym
        real(DP),allocatable,dimension(:) :: zc,zm

        integer, allocatable, dimension(:) :: imv,ipv
        integer, allocatable, dimension(:) :: jmv,jpv
        integer, allocatable, dimension(:) :: kmv,kpv
        integer, allocatable, dimension(:) :: jmhv
        integer, allocatable, dimension(:) :: kmc,kpc,kup,kum

        real(DP), allocatable, dimension(:) :: udx3c, udx3m
        real(DP), allocatable, dimension(:) :: g3rc, g3rm

!     Metric coefficients
        real(DP), allocatable, dimension(:) :: ap3j,ac3j,am3j
        real(DP), allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
        real(DP), allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
        real(DP), allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   



!     For FFTW and Poisson solver
        real(DP),dimension(13) :: ifx1
        integer*8 :: fwd_plan, bck_plan
        integer*8 :: fwdplan_1d, bckplan_1d
        real(DP), allocatable, dimension(:) :: ao,ap
        real(DP), allocatable, dimension(:) :: ak1,ak2
        real(DP), allocatable, dimension(:) :: amphk, acphk, apphk


!     Auxillary variables
        real(DP):: time
        integer :: ntime
        
        real(DP):: al, ga, ro, beta
        real(DP):: pi

        real(DP), dimension(3) :: gam,rom,alm
        real(DP), dimension(3) :: vmax

      end module param

!     -----------------------------------------------------

      module local_arrays
      use param
      implicit none

      real,allocatable,dimension(:,:,:) :: vx,vy,vz
      real,allocatable,dimension(:,:,:) :: pr,rhs
      real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3
      real,allocatable,dimension(:,:,:) :: qcap,dph,dq
      
      end module local_arrays

!     -----------------------------------------------------
