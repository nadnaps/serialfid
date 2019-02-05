!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : TimeMarcher.F90
!     Contains : Integration routine 
!
!     -----------------------------------------------------
!     -----------------------------------------------------      


      subroutine TimeMarcher
      use param
      use local_arrays
      implicit none
      integer :: ns, inp, ntr, nsub


      beta=dt/ren*0.5d0

      vx = 0.0

      do ns=1,nsst

            al=alm(ns)
            ga=gam(ns)
            ro=rom(ns)

!            call ExplicitTermsVX
            call ExplicitTermsVY
            !print*,'h2 ',sum(vx),sum(vy),sum(vz)
            call ExplicitTermsVZ
            !print*,'h3 ',sum(vx),sum(vy),sum(vz)

!            call ImplicitAndUpdateVX
            call ImplicitAndUpdateVY  
                        !print*,'i2 ',sum(vx),sum(vy),sum(vz)
   
            call ImplicitAndUpdateVZ  
                        !print*,'i3 ',sum(vx),sum(vy),sum(vz)
  

            call CalcLocalDivergence 
            call SolvePressureCorrection

        
            call CorrectVelocity                 !! SOLENOIDAL VEL FIELD
            call CorrectPressure  




      end do

      return
      end subroutine TimeMarcher



            
