!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : SolveTrid2Z.F90
!     Contains : inversion of q2 momentum equation  
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine SolveTridXYZ(q,ind)
      use param
      use local_arrays, only : rhs
      implicit none
      real(DP), intent(inout) :: q(1:n1,1:n2,1:n3)
      real(DP), dimension(n3) :: amkl,apkl,ackl,fkl
      integer, intent(in) :: ind
      integer :: jc,kc,info,ipkv(n3m),ic,n
      real(DP) :: betadx,ackl_b
      real(DP) :: amkT(n3m-1),ackT(n3m),apkT(n3m-1),appk(n3-3)
      real(DP),allocatable :: rhst(:,:,:)


      allocate(rhst(1:n3m,1:n1m,1:n2m))


!	routine to transpose while running in parallel
      do kc=1,n3m
       do jc=1,n2m
        do ic=1,n1m

          rhst(kc,ic,jc) = rhs(ic,jc,kc)

        end do 
       end do
      end do

      betadx=beta*al
       

      do jc=1,n2m
         do ic=1,n1m
          do kc=1,n3m
            ackl_b=1.0d0/(1.-ac3sk(kc)*betadx)
            amkl(kc)=-am3sk(kc)*betadx*ackl_b
            ackl(kc)=1.0d0
            apkl(kc)=-ap3sk(kc)*betadx*ackl_b
            fkl(kc)=rhst(kc,ic,jc)*ackl_b
          end do

          amkT=amkl(2:n3m)
          apkT=apkl(1:(n3m-1))
          ackT=ackl(1:n3m)

          call dgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)

          call dgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)

          do kc=1,n3m
            rhst(kc,ic,jc)=fkl(kc)
          end do
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
      q(ic,jc,kc) = q(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)


      return
      end subroutine SolveTridXYZ