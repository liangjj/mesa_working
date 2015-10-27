*deck @(#)symord.f	5.1  11/6/94
      subroutine symord(vec,scr,nbf,nob,ksym,kksym)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)symord.f	5.1   11/6/94
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
      implicit real*8 (a-h,o-z)
      dimension vec(nbf,nob),scr(nbf,nob),kksym(nob)
c
c   reorder the orbitals in symmetry order
c
      kx=0
      do 1 i=1,ksym
         do 2 j=1,nob
            if(kksym(j).eq.i) then
               kx=kx+1
               call vmove(scr(1,kx),vec(1,j),nbf)
            endif
 2       continue
 1    continue
c
      ntot=nbf*nob
c
      call vmove(vec,scr,ntot)
c
      return
      end
