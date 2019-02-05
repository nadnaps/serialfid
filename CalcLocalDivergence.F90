!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : CalcLocalDivergence.F90
!     Contains : calculate local divergence of the inter -
!                mediate velocity field for correction step
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine CalcLocalDivergence
      use param
      use local_arrays, only: vx,vy,vz,dph
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real(DP):: usdtal,dqcap   

      usdtal = 1.d0/(dt*al)

      do kc=1,n3m
        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
                    +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
                    +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
              dph(ic,jc,kc)=dqcap*usdtal
            enddo
         enddo

      enddo
      
      return
      end subroutine CalcLocalDivergence


      
