!     -----------------------------------------------------
!     -----------------------------------------------------
!
!     File     : ReadInputFile.F90
!     Contains : Subroutine for reading in input file 
!
!     -----------------------------------------------------
!     -----------------------------------------------------

      subroutine ReadInputFile
      use param
      implicit none

      character(len=4) :: dummy

      open(unit=15,file='input.in',status='old')
      read(15,301) dummy
      read(15,*)   n1m,n2m,n3m
      read(15,301) dummy
      read(15,*)   nsst,nwrit,nread
      read(15,301) dummy
      read(15,*)   ntst,tframe,tpin,tmax,walltimemax,ireset
      read(15,301) dummy
      read(15,*)   xlen,ylen,zlen
      read(15,301) dummy
      read(15,*)   istr3, str3
      read(15,301) dummy
      read(15,*)   ren
      read(15,301) dummy
      read(15,*)   ubcbot, ubctop
      read(15,301) dummy
      read(15,*)   idtv,dt,dtmax,cflmax,cfllim,resid
      read(15,301) dummy
      read(15,*)   tsta,starea
301   format(a4)
      close(15)

      n1=n1m+1                                                          
      n2=n2m+1  
      n3=n3m+1
      n1mh = n1m/2 + 1
      n2mh = n2m/2 + 1
      
      return      
      end subroutine ReadInputFile 

!     -----------------------------------------------------
!     -----------------------------------------------------


