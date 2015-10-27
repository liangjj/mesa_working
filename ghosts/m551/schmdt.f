*deck %W%  %G%
      subroutine schmdt(c,n,mdim,nroot,thrn,thrs,nvec)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cmp   extended dummy c
      dimension c(2)
c
      common /io/ inp,iout
c
ccc
c     orthogonalize new vectors to the old expansion space
ccc
      ipass=0
cc
    1 continue
cc
      ipass=ipass+1
      ipack=0
      ioff=mdim*n
c
      newvec=0
      mold=mdim
c
      do 90 m=1,nvec
c
         is=ioff+1
         ie=ioff+n
c
c
         xx=0.d0
         do 10 i=is,ie
            xx=xx+c(i)*c(i)
 10      continue
c
c     write(iout,9015) xx,thrs,thrn
c9015 format(' xx thrs thrn',3(2x,d15.7))
c
         xx=sqrt(xx)
         if(xx.gt.thrn)go to 15
         ipack=ipack+1
         if(m.eq.nvec.and.xx.gt.thrs.and.newvec.eq.0)go to 15
         if(m.eq.nvec)go to 90
ccc
c     pack the new vector set
ccc
         je=(nvec-m)*n+ioff
         do 11 i=is,je
            c(i)=c(i+n)
 11      continue
         go to 90
 15      continue
c
         xnorm=1.d0/xx
c

         do 20 i=is,ie
            c(i)=c(i)*xnorm
 20      continue
c
         joff=0
         do 40 j=1,mold
            xx=0.d0
            do 30 k=1,n
               xx=xx+c(joff+k)*c(ioff+k)
 30         continue
            do 35 k=1,n
               c(ioff+k)=c(ioff+k)-xx*c(joff+k)
 35         continue
            joff=joff+n
 40      continue
c
c
         xx=0.d0
         do 50 i=is,ie
            xx=xx+c(i)*c(i)
 50      continue
c
         xx=sqrt(xx)
         if(xx.gt.thrn)go to 55
         ipack=ipack+1
         if(m.eq.nvec.and.xx.gt.thrs.and.newvec.eq.0)go to 55
         if(m.eq.nvec)go to 90
ccc
c     pack the new vector set
ccc
         je=(nvec-m)*n+ioff
         do 51 i=is,je
            c(i)=c(i+n)
 51      continue
         go to 90
 55      continue
c
         xnorm=1.d0/xx
c
         do 60 i=is,ie
            c(i)=c(i)*xnorm
 60      continue
c
         ioff=ioff+n
         newvec=newvec+1
         mold=mold+1
c
 90   continue
c
      nvec=newvec
c     write(iout,91) nvec,ipass
      if(nvec.eq.0)go to 2000
      if(ipass.eq.1.and.ipack.gt.0)go to 1
c
c  91 format('  schmdt  nvec ipass',2i8)
c
      return
 1000 continue
      write(iout,1010)
 1010 format(/,'  can not normalize new trial vector in schmdt...stop')
      call lnkerr(' ')
 2000 continue
      write(iout,2010)
 2010 format(/,'  stop .. no new vectors obtained after schmidt')
      call lnkerr(' ')
      end
