*deck %W%  %G%
      subroutine expand(evec,c,ipt,n,nsml,nroot)
c
c***begin prologue     expand
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
c***end prologue       expand
c
      implicit real*8(a-h,o-z)
c
      dimension evec(2),c(2),ipt(2)
c
      ioff=0
      joff=0
c
      do 30 m=1,nroot
         do 10 i=1,n
            c(ioff+i)=0.d0
 10      continue
         do 20 j=1,nsml
            c(ioff+ipt(j))=evec(joff+j)
 20      continue
         ioff=ioff+n
         joff=joff+nsml
 30   continue
c
c
      return
      end
