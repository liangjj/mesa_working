*deck @(#)makvec.f	1.1  11/30/90
      subroutine makvec(c,b,evec,n,mdim,nroot)
c
c***begin prologue     makvec
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)makvec.f	1.1   11/30/90
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       makvec
c
      implicit real*8(a-h,o-z)
c
      dimension c(2),b(2),evec(2)
c
      loff=0
      ioff=0
      do 30 m=1,nroot
         is=ioff+1
         ie=ioff+n
         do 10 i=is,ie
            b(i)=0.d0
 10      continue
         joff=0
         do 20 j=1,mdim
            aa=evec(j+loff)
            do 15 i=is,ie
               joff=joff+1
               b(i)=b(i)+aa*c(joff)
 15         continue
 20      continue
         ioff=ioff+n
         loff=loff+mdim
 30   continue
c
c
      return
      end
