*deck @(#)blfold.f	5.1  11/6/94
      subroutine blfold(tr,sq,n)
c
c***begin prologue     blfold
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)blfold.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       blfold
c
      implicit integer(a-z)
      real*8 tr(*),sq(n,n)
c
      ix=0
      do 1 i=1,n
         do 2 j=1,i
            ix=ix+1
            tr(ix)=sq(i,j)+sq(j,i)
 2       continue
 1    continue
c
c
      return
      end
