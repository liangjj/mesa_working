      subroutine fixbin(xbin,bin,nocc,mij)
c
c***begin prologue     fixbin
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
c***end prologue       fixbin
c
      dimension bin(2),xbin(2)
c
      ntr=(nocc*(nocc+1))/2
      noc2=nocc*nocc
c
      ix=1
      jx=1
      do 10 i=1,mij
         call totr(xbin(ix),bin(jx),nocc,ntr)
         ix=ix+ntr
         jx=jx+noc2
 10   continue
c
      return
      end
