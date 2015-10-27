*deck  @(#)trn1so.f	5.1 11/6/94
      subroutine trn1so(c,nbf,norbs,numij,t1,t2,t3,bflabl)
      implicit integer (a-z)
c***begin prologue     trn1so.f
c***date written       910101 (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           spin-orbit 
c***author             braunstein, matt
c***source             @(#)trn1so.f	5.1   11/6/94
c***purpose            transforms one-electron spin-orbit matrices. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       trn1so.f
c
      real*8 c(nbf,norbs),t1(norbs,nbf),t2(norbs,norbs)
      real*8 t3(nbf,nbf)
      character*(*) bflabl(nbf)
      logical logkey
c
      common/io/inp,iout
c
c     ---- transform one-electron spin-orbit integrals ----
      write(iout,*)' transforming one-electron spin-orbit integrals'
      call iosys('read character "basis function labels" from rwf',
     $            -1,0,0,bflabl)
c
c     ----- lx -----    
      call iosys('read real "asox integrals" from rwf',
     $            nbf*nbf,t3,0,' ')
c
c     do transformation; have to transform full matrix, but only
c     need lower triangle of the output transformed matrix since
c     m904 assumes this matrix to be antisymmetric
c
      call ebtc(t1,c,t3,norbs,nbf,nbf)
      call ebc(t2,t1,c,norbs,nbf,norbs)
      call totrso(t3,t2,norbs,numij)
c       
      call iosys('write real "msox integrals" to tints',
     $           numij,t3,0,' ')
c
c     ----- print the integrals ----
c
      if(logkey(ops,'print=int=so',.false.,' '))then
         write(iout,*)'molecular spin-orbit integrals - x'
         call wmat(t2,norbs,norbs,bflabl,bflabl)
      end if
c
c     ----- ly -----
      call iosys('read real "asoy integrals" from rwf',
     $           nbf*nbf,t3,0,' ')
      call ebtc(t1,c,t3,norbs,nbf,nbf)
      call ebc(t2,t1,c,norbs,nbf,norbs)
      call totrso(t3,t2,norbs,numij)
c    
      call iosys('write real "msoy integrals" to tints',
     $           numij,t3,0,' ')
c
      if(logkey(ops,'print=int=so',.false.,' '))then
         write(iout,*)'molecular spin-orbit integrals - y'
         call wmat(t2,norbs,norbs,bflabl,bflabl)
      end if
c
c     ----- lz -----
      call iosys('read real "asoz integrals" from rwf',
     $            nbf*nbf,t3,0,' ')
      call ebtc(t1,c,t3,norbs,nbf,nbf)
      call ebc(t2,t1,c,norbs,nbf,norbs)
      call totrso(t3,t2,norbs,numij)
c       
      call iosys('write real "msoz integrals" to tints',
     $            numij,t3,0,' ')
      if(logkey(ops,'print=int=so',.false.,' '))then
         write(iout,*)'molecular spin-orbit integrals - z'
         call wmat(t2,norbs,norbs,bflabl,bflabl)
      end if
c
c
      return
      end
