!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : SolvePressureCorrection.F90
!     Contains : perform the calculation of dph, periodic 
!                direction use the fourier transform
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine SolvePressureCorrection
      use param
      use local_arrays, only: dph
      implicit none
      integer :: i,j,k,info
      real(DP) :: coefnorm,acphT_b
      real(DP) :: xr(n2m,n1m)
      complex(DP) :: xa(n2mh,n1m)
      real(DP) :: appph(n3m-2)
      real(DP) :: amphT(n3m-1), apphT(n3m-1)
      real(DP), dimension(n3m) :: acphT,drhs,apph,amph
      integer(DP) :: phpiv(n3m)      
      integer :: jmh
      real(DP),allocatable,dimension(:,:,:) :: dpht,dpho

      allocate(dpht(1:n3m,1:n1m,1:n2m+2))
      allocate(dpho(1:n1m,1:n2m+2,1:n3m))



      coefnorm = 1.d0/(dble(n1m)*dble(n2m))

!   fft applied to the x2 direction to the
!   complex coeff. from cos fft
!   from physical to wave number space

      do k=1,n3m

        do j=1,n2m
          do i=1,n1m
           xr(j,i)=dph(i,j,k)
          enddo
        enddo
        
        call dfftw_execute_dft_r2c(fwd_plan,xr,xa)

        do j=1,n2mh
         do i=1,n1m
         dpho(i,j,k)=dreal(xa(j,i))*coefnorm
         dpho(i,j+n2mh,k)=dimag(xa(j,i))*coefnorm
        enddo
        enddo

      end do

!   transposing elements
!   a separate routine required while running across processors
      do k=1,n3m
        do j=1,n2m+2
          do i=1,n1m
      
            dpht(k,i,j) = dpho(i,j,k)
      
          end do
        end do 
      end do

!     inversion of the matrix in the x2 and x3 directions (real part)

      do j=1,n2m+2
        jmh=jmhv(j)
        do i=1,n1m
         do k = 1,n3m
          acphT_b=1.0/(acphk(k)-ak2(jmh)-ak1(i))
          drhs(k)=dpht(k,i,j)*acphT_b
          apph(k)=apphk(k)*acphT_b
          amph(k)=amphk(k)*acphT_b
          acphT(k)=1.0d0
         enddo
  
         amphT=amph(2:n3m)
         apphT=apph(1:(n3m-1))
  
         call dgttrf(n3m, amphT, acphT, apphT, appph, phpiv, info)
  
         call dgttrs('N',n3m,1,amphT,acphT,apphT,appph,phpiv,drhs, &
                      n3m, info)
  
          do k=1,n3m
            dpht(k,i,j) = drhs(k)
          enddo
        enddo
      enddo

!   transpose back again
!   again, a separate routine while running in parallel
      do k=1,n3m
        do j=1,n2m+2
          do i=1,n1m
      
            dpho(i,j,k) = dpht(k,i,j)
      
          end do
        end do 
      end do

!   inverse fft applied to the phi x1 direction
!   from wave number space to physical space

      do k=1,n3m

       do j=1,n2mh
        do i=1,n1m
          xa(j,i)=dcmplx(dpho(i,j,k),dpho(i,j+n2mh,k))
        enddo
       end do

      call dfftw_execute_dft_c2r(bck_plan,xa,xr)

       do j=1,n2m
         do i=1,n1m
           dph(i,j,k)=xr(j,i)
         enddo
       end do
      end do

      if(allocated(dpht)) deallocate(dpht)
      if(allocated(dpho)) deallocate(dpho)

      return
      end subroutine SolvePressureCorrection
