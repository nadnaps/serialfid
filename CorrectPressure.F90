!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : CorrectPressure.F90
!     Contains : calculates the pressure. evaluated at 
!                centre of box  
!               p^{n+1}=p^{n}+phi^{n+1}-b*Nabla^2 phi^{n+1}
!     -----------------------------------------------------
!     -----------------------------------------------------


      subroutine CorrectPressure
      use param
      use local_arrays, only: pr,dph
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real(DP):: be,amm,acc,app

      be=al*beta
      do kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        amm=amphk(kc)
        acc=acphk(kc)
        app=apphk(kc)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
            pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*(              &
              (dph(ip,jc,kc)-2.d0*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q   &
             +(dph(ic,jp,kc)-2.d0*dph(ic,jc,kc)+dph(ic,jm,kc))*dx2q   &
             +(dph(ic,jc,kp)*app+dph(ic,jc,kc)*acc+dph(ic,jc,km)*amm))
          enddo
        enddo
      enddo

      return
      end subroutine CorrectPressure
