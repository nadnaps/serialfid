!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : ExplicitTermsVY.F90
!     Contains : explicit non-linear terms  
!
!     -----------------------------------------------------
!     -----------------------------------------------------


      subroutine ExplicitTermsVY
      use param
      use local_arrays, only: vx,vy,vz,dph
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real(DP)    :: h22,h23,udx1,udx2,h21

      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

!   h term for the vy momentum equation at i+1/2,j,k+1/2

      do kc=1,n3m
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)


!     h21 term here for a 3D case
      
!     vy vy term

            h22=( (vy(ic,jp,kc)+vy(ic,jc,kc)) &
                 *(vy(ic,jp,kc)+vy(ic,jc,kc)) &
                 -(vy(ic,jm,kc)+vy(ic,jc,kc)) &
                 *(vy(ic,jm,kc)+vy(ic,jc,kc)) &
                )*udx2

!     vy vz term


            h23=( (vz(ic,jc,kp)+vz(ic,jm,kp))   &
                  *(vy(ic,jc,kpp)+vy(ic,jc,kc)) &
                 -(vz(ic,jc,kc)+vz(ic,jm,kc))   &
                  *(vy(ic,jc,kc)+vy(ic,jc,kmm)) &
                )*udx3m(kc)*0.25d0


            dph(ic,jc,kc)=-(h22+h23) +1.0   !+ h21 for 3D + any explicit body force

          enddo
        enddo
      enddo
      
      return
      end subroutine ExplicitTermsVY


