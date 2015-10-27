*deck @(#)subvar.f	5.1  11/6/94
      subroutine subvar(bl,alpha,beta,lbl,lalpha,lbeta,values,nz,nvar)
c***begin prologue     subvar
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           z-matrix, variables
c***author             martin, richard (lanl)
c***source
c***purpose            substitutes values of variables into the z-matrix.
c***description
c     call subvar(bl,alpha,beta,lbl,lalpha,lbeta,values,nz,nvar)
c***references
c***routines called    (none)
c***end prologue       subvar
      implicit integer(a-z)
      real*8 bl(nz),alpha(nz),beta(nz),values(nvar)
      integer lbl(nz),lalpha(nz),lbeta(nz)
c
      common/io/inp,iout
c
 1000 format(1x,'variable index of ',i4,' on card',i3,' is out of ',
     $       'range.  nvar=',i3)
c
c
c     lbl,lalpha, and lbeta are arrays for mapping from the lists
c     of variables into the z-matrix.
c     for example, if lbl(i)=0, then bl(i) is all ready.
c     if lbl(i)= 4, then bl(i)= values(4).
c     if lbl(i)=-4, then bl(i)=-values(4)
c
c     if there are no variables the z-matrix is ready to go.
      if(nvar.eq.0) return
c
      do 10 i=1,nz
         if(lbl(i).eq.0) then
         else if(lbl(i).lt.0) then
            bl(i)=-values(abs(lbl(i)))
         else
            bl(i)=values(lbl(i))
         endif
c
         if(lalpha(i).eq.0) then
         else if(lalpha(i).lt.0) then
            alpha(i)=-values(abs(lalpha(i)))
         else
            alpha(i)=values(lalpha(i))
         endif
c
         if(lbeta(i).eq.0) then
         else if(lbeta(i).lt.0) then
            beta(i)=-values(abs(lbeta(i)))
         else
            beta(i)=values(lbeta(i))
         endif
c
   10 continue
c
c
      return
      end
