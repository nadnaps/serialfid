!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : SolveTridZ.F90
!     Contains : inversion of q3 momentum equation  
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine SolveTridZ
      use param
      use local_arrays, only : vz,rhs
      implicit none
      real(DP), dimension(n3) :: amkl,apkl,ackl, fkl
      real(DP) :: amkT(n3-1),apkT(n3-1)
      real(DP) :: appk(n3-2)
      real(DP) :: ackT(n3)
      integer(DP) :: jc,kc,info,ic,n
      integer(DP) :: ipkv(n3)
      real(DP) :: betadx,ackl_b
      real(DP), allocatable, dimension(:,:,:) :: rhst
      

      allocate(rhst(1:n3,1:n1,1:n2m))

!	routine to transpose while running in parallel
      do kc=1,n3
       do jc=1,n2m
        do ic=1,n1

          rhst(kc,ic,jc) = rhs(ic,jc,kc)
        end do 
       end do
      end do
      

      betadx=beta*al

            fkl(1)= 0.d0
            amkl(1)=0.d0
            apkl(1)=0.d0
            ackl(1)=1.d0


            amkl(n3)=0.d0
            apkl(n3)=0.d0
            ackl(n3)=1.d0
            fkl(n3)= 0.d0

      do jc=1,n2m
          do ic=1,n1m
            do kc=2,n3m
             ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
             amkl(kc)=-am3ck(kc)*betadx*ackl_b
             ackl(kc)=1.0d0
             apkl(kc)=-ap3ck(kc)*betadx*ackl_b
             fkl(kc) = rhst(kc,ic,jc)*ackl_b
            enddo

           amkT=amkl(2:n3)
           apkT=apkl(1:(n3-1))
           ackT=ackl(1:n3)
     
           call dgttrf(n3,amkT,ackT,apkT,appk,ipkv,info)
          
           call dgttrs('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)

            do kc=1,n3
              rhst(kc,ic,jc)= fkl(kc)
            enddo
          enddo
      end do


!	transpose back
      do kc=1,n3m
       do jc=1,n2m
        do ic=1,n1m

          rhs(ic,jc,kc) = rhst(kc,ic,jc)

        end do 
       end do
      end do

      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      vz(ic,jc,kc) = vz(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)


      return
      end subroutine SolveTridZ
