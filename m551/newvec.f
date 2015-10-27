*deck @(#)newvec.f	1.1  11/30/90
      subroutine newvec(diag,eval,b,c,evec,n,nroot,mdim,conv)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)newvec.f	1.1   11/30/90
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
cmp   extended dummy diag,eval,b,c,evec
      dimension diag(2),eval(2),b(2),c(2),evec(2),conv(20)
c
      loff=0
      ioff=mdim*n
      do 50 m=1,nroot
         energy=eval(m)
c
         is=ioff+1
         ie=ioff+n
         do 10 i=is,ie
            c(i)=0.d0
 10      continue
c
         joff=0
         do 20 j=1,mdim
            aa=evec(j+loff)
            do 15 i=1,n
               c(ioff+i)=aa*(b(joff+i)-energy*c(joff+i))+c(ioff+i)
 15         continue
            joff=joff+n
 20      continue
c
         conx=0.d0
         do 30 i=is,ie
            conx=conx+c(i)*c(i)
 30      continue
         conv(m)=sqrt(conx)/float(n)
c
c
         do 40 i=1,n
            c(ioff+i)=c(ioff+i)/(diag(i)-energy)
 40      continue
c
         loff=loff+mdim
         ioff=ioff+n
 50   continue
c
c     call eigout(c(is),eval(nroot),n,1)
      return
      end
