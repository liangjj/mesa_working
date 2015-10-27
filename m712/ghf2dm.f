*deck @(#)ghf2dm.f	5.1  11/6/94
      subroutine ghf2dm(d2,dij,dkl,dik,djl,dil,djk,nprimi,nprimj,
     #                  nprimk,npriml,nfi,nfj,nfk,nfl,imgrp,jmgrp,
     #                  kmgrp,lmgrp,lend2,index,nv,len,
     #                  alpha,beta,nshell,ndmat)
c
c***begin prologue     ghf2dm
c***date written       860820   (yymmdd)
c***revision date      891201   (yymmdd)
c
c                      december 1, 1989     bhl at llnl
c                        bug fixed in the definition of ijmg and klmg
c
c***keywords
c***author             saxe, paul    (lanl)
c***source             @(#)ghf2dm.f	5.1   11/6/94
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue
c
      implicit integer (a-z)
c
      real*8 dij(nprimi,nprimj,nfi,nfj,ndmat)
      real*8 dkl(nprimk,npriml,nfk,nfl,ndmat)
      real*8 dik(nprimi,nprimk,nfi,nfk,ndmat)
      real*8 djl(nprimj,npriml,nfj,nfl,ndmat)
      real*8 dil(nprimi,npriml,nfi,nfl,ndmat)
      real*8 djk(nprimj,nprimk,nfj,nfk,ndmat)
      real*8 d2(lend2)
      real*8 alpha(nshell,nshell),beta(nshell,nshell),alp,bet
      integer index(len,6)
c
      common /io/ inp,iout
c
      ii=max(imgrp,jmgrp)
      jj=min(imgrp,jmgrp)
      kk=max(kmgrp,lmgrp)
      ll=min(kmgrp,lmgrp)
      ijmg=ii*(ii-1)/2+jj
      klmg=kk*(kk-1)/2+ll
c
c..      ijmg=imgrp*(imgrp-1)/2+jmgrp
c..      klmg=kmgrp*(kmgrp-1)/2+lmgrp
c.
c
      call rzero(d2,lend2)
c
      if (imgrp.eq.jmgrp.and.imgrp.eq.kmgrp.and.imgrp.eq.lmgrp)
     #                                                          then
c
c        ----- [ii;ii] -----
c
         do 116 ityp=1,ndmat
            do 115 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
               bet=beta(ityp,jtyp)
c
               pt=0
               do 114 if=1,nfi
                  do 113 jf=1,nfj
                     do 112 kf=1,nfk
                        do 111 lf=1,nfl
                           do 101 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     #                          dij(index(i,1),index(i,2),if,jf,ityp)*
     #                          dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     #                                         alp+
     #                          dik(index(i,1),index(i,3),if,kf,ityp)*
     #                          djl(index(i,2),index(i,4),jf,lf,jtyp)*
     #                                         bet
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
c
         do 216 ityp=1,ndmat
            do 215 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
               bet=beta(ityp,jtyp)
c
               pt=0
               do 214 if=1,nfi
                  do 213 jf=1,nfj
                     do 212 kf=1,nfk
                        do 211 lf=1,nfl
                           do 201 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     #                         2*dij(index(i,1),index(i,2),if,jf,ityp)*
     #                           dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     #                                          alp+
     #                         2*dik(index(i,1),index(i,3),if,kf,ityp)*
     #                           djl(index(i,2),index(i,4),jf,lf,jtyp)*
     #                                          bet
  201                      continue
  211                   continue
  212                continue
  213             continue
  214          continue
  215       continue
  216    continue
c
      else if (imgrp.ne.jmgrp.and.kmgrp.ne.lmgrp.and.
     #         ijmg.ne.klmg) then
c
c        ----- [ij;kl] ------
c
         do 316 ityp=1,ndmat
            do 315 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
               bet=beta(ityp,jtyp)
c
               pt=0
               do 314 if=1,nfi
                  do 313 jf=1,nfj
                     do 312 kf=1,nfk
                        do 311 lf=1,nfl
                           do 301 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     #                            dij(index(i,1),index(i,2),if,jf,ityp)*
     #                            dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     #                                           8*alp+
     #                           (dik(index(i,1),index(i,3),if,kf,ityp)*
     #                            djl(index(i,2),index(i,4),jf,lf,jtyp)+
     #                            dil(index(i,1),index(i,4),if,lf,ityp)*
     #                            djk(index(i,2),index(i,3),jf,kf,jtyp))
     #                                           *4*bet
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
c
         do 416 ityp=1,ndmat
            do 415 jtyp=1,ndmat
               alp=alpha(ityp,jtyp)
               bet=beta(ityp,jtyp)
c
               pt=0
               do 414 if=1,nfi
                  do 413 jf=1,nfj
                     do 412 kf=1,nfk
                        do 411 lf=1,nfl
                           do 401 i=1,nv
                              pt=pt+1
                              d2(pt)=d2(pt)+
     #                           dij(index(i,1),index(i,2),if,jf,ityp)*
     #                           dkl(index(i,3),index(i,4),kf,lf,jtyp)*
     #                                          4*alp+
     #                          (dik(index(i,1),index(i,3),if,kf,ityp)*
     #                           djl(index(i,2),index(i,4),jf,lf,jtyp)+
     #                           dil(index(i,1),index(i,4),if,lf,ityp)*
     #                           djk(index(i,2),index(i,3),jf,kf,jtyp))
     #                                          *2*bet
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
c$$$      lenblk=nfi*nfj*nfk*nfl
c$$$      npint=nprimi*nprimj*nprimk*npriml
c$$$      write (iout,500)
c$$$ 500  format (10x,'two-particle density matrix for hf')
c$$$      call matout(d2,nv,lenblk,nv,lenblk,iout)
c
c
      return
      end
