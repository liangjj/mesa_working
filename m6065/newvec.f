*deck @(#)newvec.f	1.1 9/8/91
c***begin prologue     newvec
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           transformation, kohn
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            rearrange mo transformation matrix for optical
c***                   potential. the output matrix contains only the
c***                   nl2 scattering orbitals and the non-zero ao's.
c***
c***references         none
c
c***routines called
c***end prologue       newvec
      subroutine newvec(trans,tmp,list,ncon,nmo,nkept,nmokpt)
      implicit integer (a-z)
      real *8 trans, tmp
      dimension trans(ncon,nmo), tmp(nkept,nmokpt), list(nkept)
      common /io/ inp,iout
      nsmall=nmo-nmokpt
      call rzero(tmp,nkept*nmokpt)
      do 10 i=1,nkept
         ii=list(i) 
         count=0
         do 20 j=nsmall+1,nmo
            count=count+1
            tmp(i,count)=trans(ii,j)
   20    continue         
   10 continue
      call copy(tmp,trans,nkept*nmokpt)
      return
      end



