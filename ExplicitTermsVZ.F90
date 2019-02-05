!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : ExplicitTermsVZ.F90
!     Contains : explicit non-linear terms  
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine ExplicitTermsVZ
      use param
      use local_arrays, only: vz,vy,vz,qcap
      implicit none
      integer :: ic,jc,kc
      integer :: km,kp,jm,jp,im,ip
      real(DP)    :: h32,h33,h31
      real(DP)    :: udx1,udx2


      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

!   h term for the vz momentum equation at i+1/2,j+1/2,k

      do kc=2,n3m
        km=kmv(kc)
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

!   	h31 term here for a 3D case


!		vy vz term

            h32=( (vy(ic,jp,kc)+vy(ic,jp,km))  &
                 *(vz(ic,jp,kc)+vz(ic,jc,kc))  &
                 -(vy(ic,jc,kc)+vy(ic,jc,km))  &
                 *(vz(ic,jc,kc)+vz(ic,jm,kc)))*udx2


!		vz vz term

            h33=( (vz(ic,jc,kp)+vz(ic,jc,kc)) &
                 *(vz(ic,jc,kp)+vz(ic,jc,kc)) &
                 -(vz(ic,jc,kc)+vz(ic,jc,km)) &
                 *(vz(ic,jc,kc)+vz(ic,jc,km)) &
                )*udx3c(kc)*0.25d0



            qcap(ic,jc,kc) = -(h32+h33) !+ h31 for 3D + any explicit body force



          enddo
        enddo
      enddo

      return
      end subroutine ExplicitTermsVZ