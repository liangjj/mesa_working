      subroutine blcg1(nk,nmk,cg,
     $     nf35,r,g,bufix,lbufso,sg,rabcx,incor)
C
C***Begin prologue     blcg1
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
c
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of 1-symmetry integrals and the corresponding
c                     density matrix to c*g.
c
c --- input
c
c     nk              no of k active orbitals.
c     nmk             no of k orbitals, including core and virtual.
c     g(--)           density matrix elements.
c     nf35            fortran no for dataset containing transformed
c                     integrals.
c
c --- working storage
c
c     r(nmk,nk)
c
c --- output
c
c     cg(nmk,nk)      c*g
c
c-----------------------------------------------------------------------
c
c
C
C***References
C
C***Routines called    (none)
C
C***End prologue       blcg1
C
      implicit real*8(a-h,o-z)
c
      dimension g(2)
      dimension cg(nmk,nk)
      real*8 bufix(lbufso),sg(*),rabcx(*)
      dimension r(nmk,nk)
c
      common /pcpack/ ipkt, nhext, ipkes
      common /io/ inp,iout
c
      ijkl = 0
      ij = 1
      ix=0
      mij=0
      nr = nk * nmk
c
      if(incor.eq.0) then
         call iosys('rewind abcx on rwf',0,0,0,0)
         call iosys('read real abcx from rwf without rewinding',
     $        lbufso,bufix,0,0)
      end if
c
      do 100 i = 1, nk
         do 200 j = 1, i
            if(incor.eq.0) then
               call vmove(r,bufix(ix+1),nr)
               mij=mij+1
               ix=ix+nr
               if (ix+nr.gt.lbufso) then
                  if (mij.lt.nij) then
                     call iosys('read real abcx from rwf '//
     $                    'without rewinding',lbufso,bufix,0,0)
                     ix=0
                  end if
               end if
            else
               call vmove(r,rabcx(ix+1),nr)
               ix=ix+nr
            end if
            kl = 1
c
            do 400 k = 1, nk
               do 500 l = 1, k
                  if (kl-ij) 540, 560, 580
 560              klij = ijkl + kl + 1
                  t = g(ijkl+1) + g(ijkl+1)
                  ijkl = ijkl + 1
                  go to 590
 540              t = g(ijkl+1)
                  ijkl = ijkl + 1
                  go to 590
 580              t = g(klij)
                  klij = klij + kl
 590              do 600 m = 1, nmk
                     cg(m,l) = cg(m,l) + t * r(m,k)
                     cg(m,k) = cg(m,k) + t * r(m,l)
 600              continue
                  kl = kl + 1
 500           continue
 400        continue
            ij = ij + 1
 200     continue
 100  continue
c
c
      return
      end
