      subroutine fasthv(buf,c,b,n,mdim,nvec,ipt)
c
c***begin prologue     fasthv
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
c     fast h*|c>=|b> routine making use of
c     sparse c vector on the first iteration
c
c***references
c
c***routines called    (none)
c
c***end prologue       fasthv
c
      implicit real*8(a-h,o-z)
c
      dimension buf(2),c(2),b(2),ipt(2)
c
      joff=0
      do 20 i=1,nvec
         js=joff+1
         je=joff+n
         do 10 j=js,je
            b(j)=0.d0
 10      continue
         joff=je
 20   continue
c
      joff=0
      koff=0
      do 60 m=1,nvec
         do 50 k=1,mdim
            i=ipt(k)
            ix=(i*(i-1))/2
            xx=c(koff+k)
            js=joff+1
            je=joff+i
            do 30 j=js,je
               ix=ix+1
               b(j)=b(j)+buf(ix)*xx
 30         continue
            if(i.eq.n)go to 45
            i1=i+1
            ij=ix+i
            do 40 j=i1,n
               b(joff+j)=b(joff+j)+buf(ij)*xx
               ij=ij+j
 40         continue
 45         continue
 50      continue
         koff=koff+mdim
         joff=joff+n
 60   continue
c
      return
      end
