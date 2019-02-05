!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : ImplicitAndUpdateVY.F90
!     Contains : implicit terms and solve the approx. fact.
!                equation  
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine ImplicitAndUpdateVY
      use param
      use local_arrays, only: vy,pr,rhs,dph,ru2
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real(DP)    :: udx2,amm,app,acc
      real(DP)    :: dcq2,dpx22
      real(DP)    :: d22q2,d33q2,d11q2
      real(DP)    :: alre,udx1q,udx2q

      
      alre=al/ren

      udx2=dx2*al
      udx1q=dx1q
      udx2q=dx2q

!  compute the rhs of the factored equation
!  everything at i+1/2,j,k+1/2
!


        do kc=1,n3m
          km=kmv(kc)
          kp=kpv(kc)
          amm=am3sk(kc)
          acc=ac3sk(kc)
          app=ap3sk(kc)
          do jc=1,n2m
           jm=jmv(jc)
           jp=jpv(jc)
            do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

!   11 second derivative of vy [for 3d case]

!
!   22 second derivative of vy
!
            d22q2=(vy(ic,jp,kc)         &
                  -2.d0*vy(ic,jc,kc)    &
                  +vy(ic,jm,kc))*udx2q   
!
!   33 second derivative of q2
!
            d33q2=vy(ic,jc,kp)*app      &
                 +vy(ic,jc,kc)*acc      &  
                 +vy(ic,jc,km)*amm      

!
!    viscid terms
!
            dcq2=d22q2+d33q2      !+d11q2 for a 3d case

!
!    component of grad(pr) along 2 direction
!    you can add an ext. pr. gradient here or also throttle it using a PI controller

            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2 !  + 1.0


!     compute the total rhs

            rhs(ic,jc,kc)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc) &
                         +alre*dcq2-dpx22)*dt


            ru2(ic,jc,kc)=dph(ic,jc,kc)
         enddo
       enddo
      enddo

      call SolveTridY (beta*al*dx2q)
      
      call SolveTridXYZ(vy(1:n1,1:n2,1:n3),2)
     
      vy(:,n2,:) = vy(:,1,:)

      return
      end subroutine ImplicitAndUpdateVY
