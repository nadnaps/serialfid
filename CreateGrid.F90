!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : CreateGrid.F90
!     Contains : routine to create grid index & coordinates
!
!     -----------------------------------------------------
!     -----------------------------------------------------      

      subroutine CreateGrid
      use param
      implicit none
      integer :: ic,jc,kc
      real(DP) :: x1,x2,x3
      real(DP) :: delet
      real(DP) :: tstr3, z2dp
      real(DP),dimension(n3) :: etaz
      real(DP),dimension(n3*2) :: etazm
      integer :: nclip, n3mo

!     Create indexing system

      allocate(imv(1:n1))
      allocate(ipv(1:n1))
      allocate(jmv(1:n2))
      allocate(jpv(1:n2))
      allocate(kmv(1:n3))
      allocate(kpv(1:n3))
      allocate(kmc(1:n3))
      allocate(kpc(1:n3))
      allocate(kup(1:n3))
      allocate(kum(1:n3))


      allocate(jmhv(1:n2+1))
!
!     direction normal to non-slip walls: Z-direction
!
      do kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=n3m
        if(kc.eq.n3m) kpv(kc)=1
      end do

      do kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      enddo
!
!   Indices for horizontal direction: Y-direction
!

      do jc=1,n2m
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=n2m
        if(jc.eq.n2m) jpv(jc)=1
      enddo

      do jc = 1,n2+1
       jmhv(jc) = mod(jc,n2m/2+1)
       if(jmhv(jc).eq.0) jmhv(jc) = n2m/2 + 1
      enddo

!
!   Indices for horizontal direction: X-direction
!
      do ic=1,n1m
        imv(ic)=ic-1
        ipv(ic)=ic+1
        if(ic.eq.1) imv(ic)=n1m
        if(ic.eq.n1m) ipv(ic)=1
      enddo


!     Create coordinate system

      allocate(xc(1:n1))
      allocate(yc(1:n2))
      allocate(zc(1:n3))

      allocate(xm(1:n1))
      allocate(ym(1:n2))
      allocate(zm(1:n3))

      allocate(g3rm(1:n3))
      allocate(g3rc(1:n3))
      allocate(udx3c(1:n3))
      allocate(udx3m(1:n3))


      dx1=real(n1m)/xlen   ! inverse of grid-distance
      dx2=real(n2m)/ylen   ! inverse of grid-distance
      dx3=real(n3m)/zlen   ! inverse of grid-distance
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      

!   uniform grid in x
      do ic=1,n1
       x1=dble(ic-1)/dble(n1m)
       xc(ic)= xlen*x1
      end do
      do ic=1,n1m
       xm(ic)=(xc(ic)+xc(ic+1))*0.5d0
      end do

!   uniform grid in y
      do jc=1,n2
       x2=dble(jc-1)/dble(n2m)
       yc(jc)= ylen*x2
      end do
      do jc=1,n2m
       ym(jc)=(yc(jc)+yc(jc+1))*0.5d0
      end do

!   uniform grid in z
      if(istr3.eq.0) then
      do kc=1,n3
       x3=dble(kc-1)/dble(n3m)
       zc(kc)=zlen*x3
      enddo
      end if

!   clustering in z: tanh
      if (istr3.eq.4) then
        tstr3 = dtanh(dble(str3))
        zc(1)=0.0d0
        do kc=2,n3
          z2dp=dble(2*kc-n3-1)/dble(n3m)
          zc(kc)=(1+dtanh(dble(str3)*z2dp)/tstr3)*0.5*zlen
          if(zc(kc).lt.0.or.zc(kc).gt.zlen)then
            write(*,*)'Force the grid in tanh: ','zc(',kc,')=',zc(kc)
            stop
          endif
        end do
      end if

!   clustering in z: clipped chebyshev
      if(istr3.eq.6) then
        nclip = str3
        n3mo = n3+nclip+nclip
        do kc=1,n3mo
          etazm(kc) = dcos(pi*(dble(kc)-0.5d0)/dble(n3mo))
        end do
        do kc=1,n3
          etaz(kc)=etazm(kc+nclip)
        end do
        delet = etaz(1)-etaz(n3)
        do kc=1,n3
          etaz(kc)=etaz(kc)/(0.50d0*delet)
        end do
        zc(1) = 0.d0
        do kc=2,n3m
          zc(kc) = zlen*(1.d0-etaz(kc))*0.5d0
        end do
        zc(n3) = zlen
        
      endif

!       
      do kc=1,n3m
       zm(kc)=(zc(kc)+zc(kc+1))*0.5d0
      enddo
      zm(n3) = zlen

!   metric quantities
      do kc=1,n3m
        g3rm(kc)=(zc(kc+1)-zc(kc))*dx3
      enddo
      do kc=2,n3m
        g3rc(kc)=(zc(kc+1)-zc(kc-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zc(2)-zc(1))*dx3
      g3rc(n3)= (zc(n3)-zc(n3m))*dx3

      do kc=1,n3m
        udx3m(kc) = dx3/g3rm(kc)
        udx3c(kc) = dx3/g3rc(kc)
      end do
      udx3c(n3) = dx3/g3rc(n3)


! save cordinates to file 
      open(unit=78,file='outfiles/xcorr.out',status='unknown')
        do kc=1,n1
          write(78,'(i5,2(2x,es16.8))') kc,xc(kc),xm(kc)
        end do
        close(78)
        open(unit=78,file='outfiles/ycorr.out',status='unknown')
        do kc=1,n2
          write(78,'(i5,2(2x,f16.8))') kc,yc(kc),ym(kc)
        end do
        close(78)
        open(unit=78,file='outfiles/zcorr.out',status='unknown')
        do kc=1,n3
          write(78,345) kc,zc(kc),zm(kc),g3rc(kc),g3rm(kc)
        end do
      close(78)
345     format(i5,4(2x,es16.8))



      return
      end subroutine CreateGrid
