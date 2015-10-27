*deck @(#)eigint.f	5.1  11/6/94
      subroutine eigint(eigvec,fre,natoms3,nvar,vname,b,evcint,
     $                  icfx,binv,order,scr6)
c***begin prologue     eigint.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)eigint.f	5.1   11/6/94
c***purpose            transform cartesian vibrational eigenvectors 
c                      to internals. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       eigint.f
      implicit none
c     --- input variables -----
      integer nvar,natoms3,icfx
c     --- input arrays (unmodified) ---
      integer order(natoms3)
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 eigvec(natoms3,natoms3),evcint(nvar,nvar)
      real*8 fre(natoms3)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 b(natoms3,nvar),binv(nvar,nvar)
      real*8 scr6(2*nvar)
c     --- local variables ---
      integer inp,iout
      integer nvc,nrt,count,i,j,icol,jj,irt
      real*8 fac,zero,five
      real*8 det
      logical debug
c
      parameter (debug=.false.,zero=0.0d+00,five=5.0d+00)
c
      common /io/ inp,iout
c
 1000 format(5x,'cartesian normal modes ')
 1010 format(/,
     $' ** warning:  different $ of zero frequencies than expected'/)
 1020 format(5x,'c to z trans matrix - inverting')
 1030 format(5x,'internal coordinate normal modes')
c
c     --- print the cartesian normal modes.
      write(iout,1000)
      call frqprt(eigvec,natoms3,natoms3,natoms3,natoms3,0,0,
     $            vname,' ',0,fre,.true.)
c
      if(icfx.eq.1) return
c
c     --- transform the eigenvectors to internal coordinates
c         read in the b matrix generated in m732
      nvc=natoms3*nvar
      call iosys('read real "ctoz transformation matrix" from rwf',
     $            nvc,b,0,' ')
c
c     --- crunch the columns of the eigenvector matrix, getting rid of
c         the zero frequency modes
      nrt=6
      count=0
      do 10 icol=1,natoms3
         if(abs(fre(icol)).gt.five) then
            count=count+1
            fre(count)=fre(icol)
            do 5 j=1,natoms3
               eigvec(j,count)=eigvec(j,icol)
    5       continue
         endif
   10 continue
      if(count.ne.natoms3-nrt) then
         write(iout,1010)
      endif
c
c     --- b contains the derivatives of cartesian coordinates     
c         with respect to internals. we need to invert this, but 
c         first contract the rectangular matrix into a square.  
c         this removes the fixed cartesian coordinates
c
c         fill up the order array which allows us to refer to the fixed
c         coordinates (order(1)-order(6)) or to the varied coordinates
c         (order(7)-order(3n))
c
      do 27 i=1,natoms3
         order(i)=i
   27 continue
      order(6)=8
      order(7)=6
      order(8)=7
c
      irt=6
      do 88 j=1,nvar
         jj=order(j+irt)
         do 87 i=1,nvar
           binv(j,i)=b(jj,i)
   87    continue
   88 continue
      if(debug) then
         write(iout,1020)
         call matout(binv,nvar,nvar,nvar,nvar,iout)
      endif
      call minvrt(binv,nvar,nvar,det,scr6(1),scr6(nvar+1))
c
c     --- now put it back
      do 86 j=1,nvar
         jj=order(j+irt)
         do 85 i=1,nvar
            b(jj,i)=binv(i,j)
   85    continue
   86 continue
c
      call ebtc(evcint,b,eigvec,nvar,natoms3,nvar)
c
c     --- normalize the normal modes expressed in z-matrix coodinates
      do 200 j=1,nvar
         fac=zero
         do 100 i=1,nvar
            fac=fac+evcint(i,j)*evcint(i,j)
  100    continue
         do 150 i=1,nvar
            evcint(i,j)=evcint(i,j)/sqrt(fac)
  150    continue
  200 continue
c
c
      write(iout,1030) 
      call frqprt(evcint,nvar,nvar,nvar,nvar,1,0,vname,' ',
     $            0,fre,.true.)
c
c
      return
      end
