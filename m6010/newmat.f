*deck @(#)newmat.f	1.1 9/7/91
c***begin prologue     newmat
c***date written       910304   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6010, link 6010, kohn data
c***author             schneider, barry (lanl)
c***source             m6010
c***purpose            re-arrange matrices
c***                   for kohn codes.
c
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       newmat
      subroutine newmat(matin,matout,scr,list,nbf,nmo,type)
      implicit integer (a-z)
      real *8 matin, matout, scr
      character *(*) type
      dimension matin(nbf,nbf), matout(nbf,nmo), scr(nmo,nmo)
      dimension list(nmo)
      if (type.eq.'ao-mo') then
          call rzero(matout,nbf*nmo)
          do 10 i=1,nmo
             ii=list(i)
             do 20 j=1,nbf
                matout(j,i)=matin(j,ii)
   20        continue
   10     continue
          call copy(matout,matin,nmo*nbf)
      elseif(type.eq.'mo-mo') then
          call rzero(scr,nmo*nmo)
          do 30 i=1,nmo
             ii=list(i)
             do 40 j=1,nmo
                jj=list(j)
                scr(i,j)=matin(ii,jj)
   40        continue
   30     continue
          call copy(scr,matin,nmo*nmo)
      else
          call lnkerr('error in newmat')  
      endif
      return
      end    
