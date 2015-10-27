*deck @(#)slaterf.f	5.2  2/5/95
      subroutine slaterf(ngrid,mxgrd,dengrida,dengridb,alpha,derivs,
     $     fout)
c***begin prologue     slaterf.f
c***date written       930521     (yymmdd)
c***revisiondate       2/5/95
c
c***keywords           xm602, link 602, DFT, LYP, gradient corrected
c
c***author             RUSSO, thomas v.    (lanl)
c***source             @(#)slaterf.f	5.2   2/5/95
c***description        computes the value of the slater functional
c                      given the value of the density on a grid of points.
c
c***references
c        Johnson, B. G. et al., J. Chem. Phys 98(7),5612
c***routines called
c                      
c
c***end prologue       slaterf.f
c
c    ----------
c death to the bloated lackey of imperialist FORTRASH, implicit typing,
c the bane of mine existence and my mortal enemy
c
      implicit none
c
c --- input variables ---
c
c derivs>0 do both functional and derivs
c derivs=0 do only functional
c derivs<0 do only derivs
      integer ngrid,derivs,mxgrd
      real*8 alpha
c
c --- input arrays (unmolested) ---
c
      real*8 dengrida(ngrid),dengridb(ngrid)
c
c --- input arrays (scribbled all over then discarded) ---
c
c --- input arrays (modified)
c --- Rule Psix ---
c (there is No Rule Psix)
c
c --- output arrays ---
      real*8 fout(mxgrd,*)
c
c local
c
      integer i
      real*8 n4a,one,pif,fr3,un3

      parameter(one=1.0d0,fr3=4.0d0/3.0d0,un3=1.0d0/3.0d0)
      n4a=9.0d0*alpha/4.0d0
      pif=(3.0d0/(16.0d0*atan(one)))**un3

      if (derivs.ge.0) then
         do 10 i=1,ngrid
            fout(i,1)=-pif*n4a*(dengrida(i)**fr3+dengridb(i)**fr3)
 10      continue 
      endif
      
      if (derivs .ne. 0) then
         do 20 i=1,ngrid
            fout(i,2)=-3.0d0*alpha*pif*dengrida(i)**un3
            fout(i,3)=0
            fout(i,4)=-3.0d0*alpha*pif*dengridb(i)**un3
            fout(i,5)=0
 20      continue 
      endif

      return
      end
