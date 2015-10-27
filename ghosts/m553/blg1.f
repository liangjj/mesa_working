      subroutine blg1(nk,nmk,lokk,lenk,mix,cm,
     $     nf35,r,g,buf,lbufso,rabcx,incor)
C
C***Begin prologue     blg1
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
c-----------------------------------------------------------------------
c
c --- description     this routine makes contributions from a block
c                     of 1-symmetry integrals to g(ijkl).  the integrals
c                     (ij kl) are stored in canonical order.
c
c --- input
c
c     nk              no of k active orbitals.
c     nmk             no of k orbitals, including core and virtual.
c     lock(nk)        lock(k) starting postion of the kth vector in
c                      the arrays mix and cm.
c     lenk(nk)        length of the kth vector.
c     mix(--)         indices of vector components.
c     cm(--)          vector components.
c     nf35            fortran no for dataset containing transformed
c                     integrals.
c
c --- working storage
c
c     r(nmk,nk)
c
c --- output
c
c     g(2)
c
c-----------------------------------------------------------------------
C
C***References
C
C***Routines called    (none)
C
C***End prologue       blg1
C
      implicit real*8(a-h,o-z)
c
      dimension lokk(2),lenk(2),mix(2),cm(nmk,nk)
      dimension g(2)
      real*8 buf(lbufso),rabcx(*)
      dimension r(nmk,nk)
c
      common /pcpack/ ipkt, nhext, ipkes
      common /io/ inp,iout
c
      if(incor.eq.0) then
         call iosys('rewind abcx on rwf',0,0,0,0)
         call iosys('read real abcx from rwf without rewinding',
     $        lbufso,buf,0,0)
      end if
c
      kijkl = 0
      ij = 1
      ix=0
      mij=0
      nnk=nk*(nk+1)/2
      nr = nk * nmk
      do 100 i = 1, nk
         do 200 j = 1, i
            if(incor.eq.0) then
               call vmove(r,buf(ix+1),nr)
               mij=mij+1
               ix=ix+nr
               if (ix+nr.gt.lbufso) then
                  if (mij.lt.nnk) then
                     call iosys('read real abcx from rwf '//
     $                    'without rewinding',lbufso,buf,0,0)
                     ix=0
                  end if
               end if
            else
               call vmove(r,rabcx(ix+1),nr)
               ix=ix+nr
            end if
            kl = 1
            do 400 k = 1, nk
               mk1 = lokk(k) + 1
               mk2 = mk1 + lenk(k) - 1
               do 500 l = 1, k
                  ml1 = lokk(l) + 1
                  ml2 = ml1 + lenk(l) - 1
                  t = 0.d0
                  if (mk2 .lt. mk1) go to 512
                  do 510 m = mk1, mk2
                     t = t + r(mix(m),l) * cm(mix(m),k)
 510              continue
 512              if (k .ne. l) go to 515
                  t = t + t
                  go to 530
 515              if (ml2 .lt. ml1) go to 530
                  do 520 m = ml1, ml2
                     t = t + r(mix(m),k) * cm(mix(m),l)
 520              continue
 530              if (kl-ij) 540, 560, 580
c------------------------------------------------c
c         make k and l contributions to g        c
c------------------------------------------------c
 560              iklij = kijkl + kl + 1
                  t = t + t
 540              g(kijkl+1) = g(kijkl+1) + t
                  kijkl = kijkl + 1
                  go to 590
c------------------------------------------------c
c         make i and j contributions to g        c
c------------------------------------------------c
 580              g(iklij) = g(iklij) + t
                  iklij = iklij + kl
 590              kl = kl + 1
 500           continue
 400        continue
            ij = ij + 1
 200     continue
 100  continue
c
c
      return
      end
