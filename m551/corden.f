*deck @(#)corden.f	1.1  11/30/90
      subroutine corden(c,cdens,nbf,nco,nnp)
c
c***begin prologue     corden
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)corden.f	1.1   11/30/90
c
c***purpose
c                      form the frozen core density
c                       matrix (cdens).
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       corden
c
      implicit integer (a-z)
c
      real*8 c(nbf,nco),cdens(nnp)
c
      call rzero(cdens,nnp)
      ncore=0
      drtorb=0
      do 3 bf=1,nco
         ncore=ncore+1
         ij=0
         do 2 i=1,nbf
            do 1 j=1,i
               ij=ij+1
               cdens(ij)=cdens(ij)+c(i,bf)*c(j,bf)
 1          continue
 2       continue
    3 continue
c
      ntot=(nbf*(nbf+1))/2
      do 4 i=1,ntot
         cdens(i)=cdens(i)*2.0
    4 continue
c
      return
      end
