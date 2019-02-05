!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : ImplicitAndUpdateVZ.F90
!     Contains : implicit terms and solve the approx. fact.
!                equation  
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine ImplicitAndUpdateVZ
      use param
      use local_arrays, only: vz,qcap,pr,ru3,rhs
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real(DP)    :: udx3
      real(DP)    :: dq32,dq33,dcq3,dpx33,dq31
      real(DP)    :: app,acc,amm
      real(DP)    :: alre,udx1q,udx2q
      
      
      alre=al/ren
      

      udx1q=dx1q
      udx2q=dx2q

!  compute the rhs of the factored equation
!  everything at i+1/2,j+1/2,k

      do kc=2,n3m
        km=kmv(kc)
        kp=kc+1
        udx3 = al*udx3c(kc)
        amm=am3ck(kc)
        acc=ac3ck(kc)
        app=ap3ck(kc)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do ic=1,n1m
              im=imv(ic)
              ip=ipv(ic)

!   11 second derivatives of q3 [for 3d case]

!   22 second derivatives of q3

            dq32=(vz(ic,jm,kc) &
                 -2.d0*vz(ic,jc,kc) &
                 +vz(ic,jp,kc))*udx2q 

!   33 second derivatives of q3

            dq33=vz(ic,jc,kp)*app &
                +vz(ic,jc,kc)*acc &
                +vz(ic,jc,km)*amm 

!   viscous terms

            dcq3=dq32+dq33! +dq31 for 3d case

!  component of grad(pr) along x3 direction

            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3

            rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc) &
                          +alre*dcq3-dpx33)*dt 

!  updating of the non-linear terms

            ru3(ic,jc,kc)=qcap(ic,jc,kc)

         enddo
       enddo
      enddo
      
      call SolveTridY (beta*al*dx2q)
      call SolveTridZ

      vz(:,:,1) = 0.0
      vz(:,:,n3)= 0.0

      return
      end subroutine ImplicitAndUpdateVZ
