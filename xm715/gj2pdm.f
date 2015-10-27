*deck %W%  %G%
      subroutine gj2pdm(d2,dij,dkl,nprimi,nprimj,
     $                  nprimk,npriml,nfi,nfj,nfk,nfl,imgrp,jmgrp,
     $                  kmgrp,lmgrp,lend2,index,nv,len,
     $                  alpha,nshell,ndmat)
c
c***begin prologue     gj2pdm.f
c***date written       940507 
c***revision date      11/6/94
c   may 7, 1994        rlm at lanl
c      modifying ghf2dm(m712) to do just the coulomb piece.
c***keywords
c***author             martin, richard    (lanl)
c***source             %W%   %G%
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       gj2pdm.f
      implicit none
c     --- input variables -----
      integer nprimi,nprimj,nprimk,npriml
      integer nfi,nfj,nfk,nfl
      integer imgrp,jmgrp,kmgrp,lmgrp
      integer lend2,nv,len,nshell,ndmat
c     --- input arrays (unmodified) ---
      integer index(len,6)
      real*8 dij(nprimi,nprimj,nfi,nfj,ndmat)
      real*8 dkl(nprimk,npriml,nfk,nfl,ndmat)
      real*8 alpha(nshell,nshell)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 d2(lend2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer ii,jj,kk,ll,ijmg,klmg
      integer ityp,jtyp,if,jf,kf,lf
      integer i,pt
      integer lenblk,npint
      logical debug
      real*8 alp
c
      common /io/ inp,iout
c
      parameter (debug=.false.)
c
      ii=max(imgrp,jmgrp)
      jj=min(imgrp,jmgrp)
      kk=max(kmgrp,lmgrp)
      ll=min(kmgrp,lmgrp)
      ijmg=ii*(ii-1)/2+jj
      klmg=kk*(kk-1)/2+ll
c
c     --- ijmg=imgrp*(imgrp-1)/2+jmgrp
c     --- klmg=kmgrp*(kmgrp-1)/2+lmgrp
c
c     --- initialize the 2-particle density matrix.
      call rzero(d2,lend2)
c
      if (imgrp.eq.jmgrp.and.imgrp.eq.kmgrp.and.imgrp.eq.lmgrp)
     $                                                          then
c
c        ----- [ii;ii] -----
         do 116 ityp=1,ndmat
            do 115 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
c
               pt=0
               do 114 if=1,nfi
                  do 113 jf=1,nfj
                     do 112 kf=1,nfk
                        do 111 lf=1,nfl
                           do 101 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     $                          dij(index(i,1),index(i,2),if,jf,ityp)*
     $                          dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     $                                         alp
  101                      continue
  111                   continue
  112                continue
  113             continue
  114          continue
  115       continue
  116    continue
c
      else if (imgrp.eq.jmgrp.and.kmgrp.eq.lmgrp) then
c
c        ----- [ii;jj] -----
         do 216 ityp=1,ndmat
            do 215 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
c
               pt=0
               do 214 if=1,nfi
                  do 213 jf=1,nfj
                     do 212 kf=1,nfk
                        do 211 lf=1,nfl
                           do 201 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     $                         2*dij(index(i,1),index(i,2),if,jf,ityp)*
     $                           dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     $                                          alp
  201                      continue
  211                   continue
  212                continue
  213             continue
  214          continue
  215       continue
  216    continue
c
      else if (imgrp.ne.jmgrp.and.kmgrp.ne.lmgrp.and.
     $         ijmg.ne.klmg) then
c
c        ----- [ij;kl] ------
         do 316 ityp=1,ndmat
            do 315 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
c
               pt=0
               do 314 if=1,nfi
                  do 313 jf=1,nfj
                     do 312 kf=1,nfk
                        do 311 lf=1,nfl
                           do 301 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     $                            dij(index(i,1),index(i,2),if,jf,ityp)*
     $                            dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     $                                           8*alp
  301                      continue
  311                   continue
  312                continue
  313             continue
  314          continue
  315       continue
  316    continue
c
      else
c
c        ----- [ij;ll]  [ii;kl]  [ij;ij] -----
         do 416 ityp=1,ndmat
            do 415 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
c
               pt=0
               do 414 if=1,nfi
                  do 413 jf=1,nfj
                     do 412 kf=1,nfk
                        do 411 lf=1,nfl
                           do 401 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     $                           dij(index(i,1),index(i,2),if,jf,ityp)*
     $                           dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     $                                          4*alp
  401                      continue
  411                   continue
  412                continue
  413             continue
  414          continue
  415       continue
  416    continue
      end if
c
c
      if(debug) then
         lenblk=nfi*nfj*nfk*nfl
         npint=nprimi*nprimj*nprimk*npriml
         write (iout,500)
  500    format (10x,'two-particle density matrix for hf')
         call matout(d2,nv,lenblk,nv,lenblk,iout)
      endif
c
c
      return
      end
