!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : CheckRoutines.F90
!     Contains : check maximum divergence in velocity field
!                and compute cfl values  
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine CheckDivergence(qmax)
      use param
      use local_arrays, only: vy,vz,vx
      implicit none
      real(DP),intent(out) :: qmax
      integer :: jc,kc,kp,jp,ic,ip
      real(DP):: dqcap,dvol,my_qmax
        

      dvol=1.d0/(dx1*dx2*dx3)
      qmax=0.d0                                                     

      do kc=1,n3m
        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
                    +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
                    +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
              qmax = max(abs(dqcap),qmax)          
      enddo
      enddo
      enddo
     
      
      return     
      end subroutine CheckDivergence


!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine CalcMaxCFL(cflm)
      use param
      use local_arrays, only: vx,vy,vz
      implicit none
      real(DP),intent(inout)    :: cflm
      integer :: j,k,jp,kp,i,ip
      real(DP) :: qcf
      
      cflm=0.00000001d0
                                                                       
      do k=1,n3m
        kp=k+1
        do j=1,n2m
          jp=j+1
          do i=1,n1m
            ip=i+1
            qcf=( abs((vx(i,j,k)+vx(ip,j,k))*0.5d0*dx1)  &
                 +abs((vy(i,j,k)+vy(i,jp,k))*0.5d0*dx2)  &
                 +abs((vz(i,j,k)+vz(i,j,kp))*0.5d0*dx3))

            cflm = max(cflm,qcf)
      enddo
      enddo
      enddo
          
      return  
      end subroutine CalcMaxCFL

!     -----------------------------------------------------
!     -----------------------------------------------------


      subroutine CheckMaxVel
      use param
      use local_arrays, only: vy,vz,vx
      implicit none
      integer :: jc,kc,ic

       vmax=-100.d0

       do kc=1,n3
        do jc=1,n2m
          do ic=1,n1m
           vmax(1) = max(vmax(1),abs(vx(ic,jc,kc)))
           vmax(2) = max(vmax(2),abs(vy(ic,jc,kc)))
           vmax(3) = max(vmax(3),abs(vz(ic,jc,kc)))
         enddo
        enddo
       enddo

      return   
      end subroutine CheckMaxVel

!     -----------------------------------------------------
!     -----------------------------------------------------


      subroutine LocateLargeDivergence
      use param
      use local_arrays, only: vy,vz,vx

      implicit none
      integer :: jc,kc,kp,jp,ic,ip
      real(DP):: dqcap
        
      do kc=1,n3m
        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
                    +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
                    +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
              if (abs(dqcap).gt.resid) then
                 write(6,*) 'Divergence exceeds resid ',ic,jc,kc
                 stop
              endif
      enddo
      enddo
      enddo
      
      return     
      end subroutine LocateLargeDivergence

!     -----------------------------------------------------
!     -----------------------------------------------------

