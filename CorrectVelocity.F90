!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : CorrectVelocity.F90
!     Contains : calculate the solenoidal velocity field
!                q(n+1)=qhat-grad(dph)*dt ,  pr=dph
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine CorrectVelocity
      use param
      use local_arrays, only: vx,vy,vz,dph
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real(DP)    :: usukm,udx2,udx1,locdph

      udx1 = al*dt*dx1
      udx2 = al*dt*dx2
      do kc=1,n3m
        km=kmv(kc)
        usukm = al*dt*udx3c(kc)

        do jc=1,n2m
          jm=jmv(jc)
          do ic=1,n1m
          im=imv(ic)
          locdph=dph(ic,jc,kc)

!         correct vx for 3d case
          vy(ic,jc,kc)=vy(ic,jc,kc)- &
            (locdph-dph(ic,jm,kc))*udx2
          vz(ic,jc,kc)=vz(ic,jc,kc)- &
            (locdph-dph(ic,jc,km))*usukm

        enddo 
       enddo
      enddo

      return
      end subroutine CorrectVelocity
