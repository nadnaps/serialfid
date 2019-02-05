!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : CompCoeff.F90
!     Contains : Coefficient for integration in wall-normal
!				 direction with non-uniform grid spacing
!
!     -----------------------------------------------------
!     ----------------------------------------------------- 

      subroutine CompCoeff
      use param
      implicit none

      integer :: k,kc,km,kp
      real(DP):: a33, a33m, a33p


      allocate(ap3ck(1:n3))
      allocate(ac3ck(1:n3))
      allocate(am3ck(1:n3))

      allocate(ap3sk(1:n3))
      allocate(ac3sk(1:n3))
      allocate(am3sk(1:n3))

      allocate(ap3ssk(1:n3))
      allocate(ac3ssk(1:n3))
      allocate(am3ssk(1:n3))

!	coefficients for q3   for x3 differentation
!   c means centered that is at k location

      am3ck(1)=0.d0
      ap3ck(1)=0.d0
      ac3ck(1)=1.d0
      am3ck(n3)=0.d0
      ap3ck(n3)=0.d0
      ac3ck(n3)=1.d0
      do kc=2,n3m
        km=kc-1
        kp=kc+1
        a33=dx3q/g3rc(kc)
        a33p=1.d0/g3rm(kc)
        a33m=1.d0/g3rm(km)
        ap3ck(kc)=a33*a33p
        am3ck(kc)=a33*a33m
        ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
      enddo


!  	coefficients for q1, q2 ap3sk,am3sk,ac3sk
!  	ap3ssk,am3ssk,ac3ssk, for x3 differentation
!   s means staggered that is at k+1/2 location

      do kc=2,n3m-1
        kp=kc+1
        km=kc-1
        a33=dx3q/g3rm(kc)
        a33p= +a33/g3rc(kp)
        a33m= +a33/g3rc(kc)
        ap3sk(kc)=a33p
        am3sk(kc)=a33m
        ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
        ap3ssk(kc)=ap3sk(kc)
        am3ssk(kc)=am3sk(kc)
        ac3ssk(kc)=ac3sk(kc)
      enddo
    
!	lower wall  bound. conditions  indicated by ubcbot
!   differentiation of sec. derivative at 1+1/2
    
      kc=1
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=0.d0
      ac3sk(kc)=-(a33p+ubcbot*a33m*2.d0)
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=-(a33p+a33m*2.d0)

    
!   upper wall  bound. conditions  indicated by ubctop
!   differentiation of sec. derivative at n3-1/2
    
      kc=n3m
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=0.d0
      am3sk(kc)=a33m
      ac3sk(kc)=-(a33m+ubctop*a33p*2.d0)
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=-(a33m+a33p*2.d0)



      return
      end subroutine CompCoeff


