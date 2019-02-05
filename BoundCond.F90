!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : BoundCond.F90
!     Contains : impose velocity boundary conditions
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine BoundCond
      use param
      use local_arrays, only: vx,vy,vz
      implicit none
      integer :: jc,kc,ic

      real dl1q, dl2q, dlf1, dlf2

!  periodic boundary condition
  
      do kc=1,n3m
        do jc=1,n2m
          q1(n1,jc,kc) = q1(1,jc,kc)
          q2(n1,jc,kc) = q2(1,jc,kc)
          q3(n1,jc,kc) = q3(1,jc,kc)
        enddo

        do ic=1,n1
          q1(ic,n2,kc) = q1(ic,1,kc)
          q2(ic,n2,kc) = q2(ic,1,kc)
          q3(ic,n2,kc) = q3(ic,1,kc)
        enddo
      enddo

!    the boundary condition at upper and lower wall
!    for free-slip boundary 

        if(ubcbot.eq.0)then
          dl1q = (zm(1)-zc(1))**2
          dl2q = (zm(2)-zc(1))**2
          dlf1 = dl2q/(dl2q-dl1q)
          dlf2 = dl1q/(dl2q-dl1q)

          do jc=1,n2
            do ic=1,n1
              q1(ic,jc,0) = dlf1*q1(ic,jc,1)-dlf2*q1(ic,jc,2)
              q2(ic,jc,0) = dlf1*q2(ic,jc,1)-dlf2*q2(ic,jc,2)
              q3(ic,jc,0) = 0.d0
              q3(ic,jc,1) = 0.d0
            enddo
          enddo

        else

          do jc=1,n2
            do ic=1,n1
              q1(ic,jc,0) = 0.d0
              q2(ic,jc,0) = 0.d0
              q3(ic,jc,0) = 0.d0
              q3(ic,jc,1) = 0.d0
            enddo
          enddo

        endif


        if(ubctop.eq.0)then
          dl1q = (zm(n3m)-zc(n3))**2
          dl2q = (zm(n3m-1)-zc(n3))**2
          dlf1 = dl2q/(dl2q-dl1q)
          dlf2 = dl1q/(dl2q-dl1q)

          do jc=1,n2
            do ic=1,n1
              q1(ic,jc,n3) = dlf1*q1(ic,jc,n3m) - dlf2*q1(ic,jc,n3m-1)
              q2(ic,jc,n3) = dlf1*q2(ic,jc,n3m) - dlf2*q2(ic,jc,n3m-1)
              q3(ic,jc,n3) = 0.d0
            enddo
          enddo

        else

          do jc=1,n2
            do ic=1,n1
              q1(ic,jc,n3) = 0.d0
              q2(ic,jc,n3) = 0.d0
              q3(ic,jc,n3) = 0.d0
            enddo
          enddo

        endif

      return
      end subroutine BoundCond

