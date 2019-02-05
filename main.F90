!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : main.F90
!     Contains : Startup program
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      program main
      use param

      implicit none
      integer :: ntstf, errorcode
      real(DP) :: cflm, dmax
      real(DP), dimension(2) :: ti
      real(DP), dimension(3) :: tin


      call ReadInputFile

      call OpenLogs

      pi=2.d0*dasin(1.d0)


!   output to screen
      write(6,110)xlen,ylen,zlen
  110 format(//,10x,'Rectangular Channel',//,10x, &
             'X-Length = ',f5.2,' Y-Length = ',f5.2,' Z-Length = ',f5.2)
      write(6,111)ren
  111 format(/,10x,'Parameters:  Re= ',e10.3)
      if(idtv.eq.1) then
      write(6,113) cflmax
  113 format(/,10x,'Variable dt and fixed CFL= ',e11.4,/)
      else
      write(6,114) dtmax,cfllim
  114 format(/,10x,'Fixed dt= ',e11.4,' max. CFL= ',e11.4,/)
      end if
!     ---




      call InitTimeMarchScheme

      call MemoryAlloc

      call InitArrays

      call HdfStart
!      
!      call InitStats
!
      call CreateGrid

      call CompCoeff
!
      call InitPressSolv

      time = 0.d0
      vmax = 0.d0 

      write(6,115)n1m,n2m,n3m
  115 format(/,10x,'Grid resolution: ','n1 = ',i3,'  n2 = ',i3, '  n3 = ',i3/)
      write(6,116)dx1,dx2,dx3
  116 format(/,10x,' dx1 = ',e10.3,' dx2 = ',e10.3,' dx3 = ',e10.3,/)


!   Initial conditions
      if(nread.eq.0) then
        write(6, '(10x, "nread set to 0: creating initial conditions",//)')
        ntime = 0
        time = 0.d0
        vmax = 0.d0 
        cflm = 0.d0
        call CreateInitCond
      else
        write(6,*) 'nread set to 1: reading in initial conditions'
        tmax = tmax + time
        call ReadInitCond
      end if


      call CheckDivergence(dmax)

      ntstf = ntst

      write(6,117) tframe,ntstf,tpin
  117 format(10x,'Time params : tframe =',f10.1, &
            '  ntstf =',i8,2x,'tpin =',f10.1//)


!   output to screen about time-steps
!      if(idtv.eq.1) then
!        write(6,*) 'ntime - time - vmax(1) - vmax(2) - vmax(3) - dt - dmax'
!        write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),dt,dmax
!      else
!        write(6,*) 'ntime - time - vmax(1) - vmax(2) - vmax(3)'
!        write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),cflm,dmax
!      endif
      cflm=cflm*dt


! --------------------------------------------
!   Start the time dependent calculation 
! --------------------------------------------

      do ntime=0,ntstf

        call CalcMaxCFL(cflm)

!    adjust timestep or exit if too small 
      if(idtv.eq.1) then
         if(ntime.ne.1) then
            dt=cflmax/cflm
            if(dt.gt.dtmax) dt=dtmax
         endif
         if(dt.lt. 0.00000001d0) errorcode=166 
      else
         cflm=cflm*dt
         if(cflm.gt.cfllim) errorcode=165 
      endif

!     March a step in time
      call TimeMarcher

      time=time+dt

!    Enter checks either on first time-step or every tpin
!    perform routine checks and output

      if((ntime.eq.1).or.(mod(time,tpin).lt.dt)) then 

        call CheckMaxVel !Make sure velocities are not exceeding maximum
        if(vmax(1).gt.1000.d0.and.vmax(2).gt.1000.d0) errorcode=266 

        call CalcMaxCFL(cflm) !Recalculate CFL 
        if(idtv.ne.1) cflm=cflm*dt 

        call CheckDivergence(dmax) !Make sure velocity is solenoidal
        if(dmax.gt.resid) errorcode=169 !If velocity field not solenoidal, stop

!    Calculate statistics if time > TSTA
        if(time.gt.tsta) then
          !call UpdateStats
          !call CalcDissipation
        endif
      
      end if

      if(time.gt.tmax) errorcode=333

      if(mod(time,tpin).lt.dt) then
          write(6,*) '------------------------------------------'
          write(6,*) 'Maximum divergence = ', dmax
          write(6,*) 'ntime - time - vmax(1) - vmax(2) - vmax(3)'
          write(6,*) ntime,time,vmax(1),vmax(2),vmax(3)
        endif

      if( ntime .eq. ntst ) errorcode = 555


      if(errorcode.ne.0) then

!   dt too small
      if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)
!   cfl too high    
      if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)
!   velocities diverged
      if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)
!   mass not conserved
      if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)
!   Physical time exceeded tmax, no error; normal quit
      if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)
!   walltime exceeded walltimemax, no error; normal quit
      if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)
!   FFT input not correct
      if(errorcode.eq.444) call QuitRoutine(tin,.false.,errorcode)
!   maximum number of timesteps reached, no error; normal quit
      if(errorcode.eq.555) call QuitRoutine(tin,.true.,errorcode)

        errorcode = 100 ! already finalized

        exit

      endif
 
      end do !end ntime loop

      call QuitRoutine(tin,.true.,errorcode)
! --------------------------------------------
!   End the time dependent calculation 
! --------------------------------------------

      stop
      end

!     -----------------------------------------------------
!     -----------------------------------------------------
