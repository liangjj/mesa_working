*deck @(#)reordr.f	1.1  11/30/90
      subroutine reordr(evec,eval,nbf,nob,temp)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)reordr.f	1.1   11/30/90
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
      dimension evec(nbf,nob),eval(nob),temp(nbf)
c
      if(nob.eq.1) return
c
      nob1=nob-1
      do 6 i=1,nob1
         ev=eval(i)
         imin=i
         i1=i+1
ccc
c     locate smallest eigenvalue
ccc
         do 2 j=i1,nob
            if(eval(j).gt.ev)go to 2
            ev=eval(j)
            imin=j
 2       continue
c
         if(imin.eq.i)go to 6
ccc
c     reorder eigenvalues
ccc
         eold=eval(i)
         eval(i)=eval(imin)
         eval(imin)=eold
ccc
c     reorder eigenvectors
ccc
         do 3 j=1,nbf
            temp(j)=evec(j,i)
 3       continue
         do 4 j=1,nbf
            evec(j,i)=evec(j,imin)
 4       continue
         do 5 j=1,nbf
            evec(j,imin)=temp(j)
 5       continue
c
 6    continue
c
      return
      end
