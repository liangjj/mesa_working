*deck  @(#)schmidt.f	5.1 11/6/94
      subroutine schmidt(a,temp,tempv,n,l)
c***begin prologue     schmidt.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)schmidt.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       schmidt.f
c
      implicit integer(a-z)
c
      real*8 a(l,n),temp(l,l),tempv(l),dum,fac
c
      common/io/inp,iout
c
c     schmidt orthogonalization assuming a unit metric.
c     a(l,n) matrix to be orthogonalized.
c     l .. dimension of basis
c     n .. number of vectors
c     temp,tempv scratch
c
      do 2 i=1,l
         do 1 j=1,l
            temp(j,i)=0.0d+00
    1    continue
         temp(i,i)=1.0d+00
    2 continue
c
      do 200 j=1,n
         if (j.eq.1) go to 100
         do 50 i=1,l
            dum=0.0d+00
            do 40 k=1,l
               dum=dum+temp(i,k)*a(k,j)
   40       continue
            tempv(i)=dum
   50    continue
         jpm=j-1
         do 80 jp=1,jpm
            dum=0.0d+00
            do 60 k=1,l
               dum=dum+tempv(k)*a(k,jp)
   60       continue
            do 70 k=1,l
               a(k,j)=a(k,j)-dum*a(k,jp)
   70       continue
   80    continue
c
c     ----- normalize j-th column -----
c
  100    continue
         fac=0.0d+00
         do 120 i=1,l
            dum=0.0d+00
            do 110 k=1,l
               dum=dum+temp(k,i)*a(k,j)
  110       continue
            fac=fac+dum*a(i,j)
  120    continue
         if (fac.le.0.0d+00) go to 999
         fac=1.0d+00/sqrt(fac)
         do 150 i=1,l
            a(i,j)=fac*a(i,j)
  150    continue
  200 continue
      return
c
c     ----- error handling -----
c
  999 continue
      write (iout,900)fac,j
  900 format (//,' schmidt: error in norm of vector:',g12.3,//
     #,       ' vector number',i5,//)
      call lnkerr(' ')
c
      end
