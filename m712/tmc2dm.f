*deck @(#)tmc2dm.f	5.1  11/6/94
      subroutine tmc2dm(d2cc,d2ac,dij,dkl,dik,djl,dil,djk,nprimi,nprimj,
     #                  nprimk,npriml,nfi,nfj,nfk,nfl,imgrp,jmgrp,
     #                  kmgrp,lmgrp,lend2,index,nv,len,
     #                  alpha,beta,nshell,ndmat)
c
c***begin prologue     tmc2dm
c***date written       871216   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield,byron    (brl)
c***source             @(#)tmc2dm.f	5.1   11/6/94
c***purpose
c    revised version of ghf2dm which separately computes
c    active-active, core-active, and core-core densities
c    for mcscf derivatives
c
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
      real*8 d2cc(lend2),d2ac(lend2),alpac,betac
      real*8 alpha(nshell,nshell),beta(nshell,nshell),alp,bet
      integer index(len,6)
c
      common /io/ inp,iout
c
      ijmg=imgrp*(imgrp-1)/2+jmgrp
      klmg=kmgrp*(kmgrp-1)/2+lmgrp
c
      call rzero(d2cc,lend2)
      call rzero(d2ac,lend2)
c
               alp=alpha(1,1)
               bet=beta(1,1)
               alpac=alpha(2,1)
               betac=beta(2,1)
c
      if (imgrp.eq.jmgrp.and.imgrp.eq.kmgrp.and.imgrp.eq.lmgrp)
     #                                                          then
c
c        ----- [ii;ii] -----
c
c
               pt=0
               do 114 if=1,nfi
                  do 113 jf=1,nfj
                     do 112 kf=1,nfk
                        do 111 lf=1,nfl
                           spt=pt
                           do 101 i=1,nv
                              pt=pt+1
                              d2cc(pt)=d2cc(pt)+
     #                          dij(index(i,1),index(i,2),if,jf,1)*
     #                          dkl(index(i,3),index(i,4),kf,lf,1)*
     #                                         alp+
     #                          dik(index(i,1),index(i,3),if,kf,1)*
     #                          djl(index(i,2),index(i,4),jf,lf,1)*
     #                                         bet
  101                      continue
c..bhl
                           pt=spt
                           do 102 i=1,nv
                              pt=pt+1
                              d2ac(pt)=d2ac(pt)+
     #                         (dij(index(i,1),index(i,2),if,jf,1)*
     #                          dkl(index(i,3),index(i,4),kf,lf,2)+
     #                          dij(index(i,1),index(i,2),if,jf,2)*
     #                          dkl(index(i,3),index(i,4),kf,lf,1))
     #                                        *alpac+
     #                         (dik(index(i,1),index(i,3),if,kf,1)*
     #                          djl(index(i,2),index(i,4),jf,lf,2)+
     #                          dik(index(i,1),index(i,3),if,kf,2)*
     #                          djl(index(i,2),index(i,4),jf,lf,1))
     #                                        *betac
  102                      continue
c..bhl
  111                   continue
  112                continue
  113             continue
  114          continue
c
      else if (imgrp.eq.jmgrp.and.kmgrp.eq.lmgrp) then
c
c        ----- [ii;jj] -----
c
               pt=0
               do 214 if=1,nfi
                  do 213 jf=1,nfj
                     do 212 kf=1,nfk
                        do 211 lf=1,nfl
                           spt=pt
                           do 201 i=1,nv
                              pt=pt+1
                              d2cc(pt)=d2cc(pt)+
     #                         2*dij(index(i,1),index(i,2),if,jf,1)*
     #                           dkl(index(i,3),index(i,4),kf,lf,1)*
     #                                          alp+
     #                         2*dik(index(i,1),index(i,3),if,kf,1)*
     #                           djl(index(i,2),index(i,4),jf,lf,1)*
     #                                          bet
  201                      continue
c..bhl
                           pt=spt
                           do 202 i=1,nv
                              pt=pt+1
                              d2ac(pt)=d2ac(pt)+
     #                        (2*dij(index(i,1),index(i,2),if,jf,1)*
     #                           dkl(index(i,3),index(i,4),kf,lf,2)+
     #                         2*dij(index(i,1),index(i,2),if,jf,2)*
     #                           dkl(index(i,3),index(i,4),kf,lf,1))
     #                                         *alpac+
     #                        (2*dik(index(i,1),index(i,3),if,kf,1)*
     #                           djl(index(i,2),index(i,4),jf,lf,2)+
     #                         2*dik(index(i,1),index(i,3),if,kf,2)*
     #                           djl(index(i,2),index(i,4),jf,lf,1))
     #                                         *betac
  202                      continue
c..bhl
  211                   continue
  212                continue
  213             continue
  214          continue
c
      else if (imgrp.ne.jmgrp.and.kmgrp.ne.lmgrp.and.
     #         ijmg.ne.klmg) then
c
c        ----- [ij;kl] ------
c

               pt=0
               do 314 if=1,nfi
                  do 313 jf=1,nfj
                     do 312 kf=1,nfk
                        do 311 lf=1,nfl
                           spt=pt
                           do 301 i=1,nv
                              pt=pt+1
                              d2cc(pt)=d2cc(pt)+
     #                            dij(index(i,1),index(i,2),if,jf,1)*
     #                            dkl(index(i,3),index(i,4),kf,lf,1)*
     #                                           8*alp+
     #                           (dik(index(i,1),index(i,3),if,kf,1)*
     #                            djl(index(i,2),index(i,4),jf,lf,1)+
     #                            dil(index(i,1),index(i,4),if,lf,1)*
     #                            djk(index(i,2),index(i,3),jf,kf,1))
     #                                           *4*bet
  301                      continue
c..bhl
                           pt=spt
                           do 302 i=1,nv
                              pt=pt+1
                              d2ac(pt)=d2ac(pt)+
     #                           (dij(index(i,1),index(i,2),if,jf,1)*
     #                            dkl(index(i,3),index(i,4),kf,lf,2)+
     #                            dij(index(i,1),index(i,2),if,jf,2)*
     #                            dkl(index(i,3),index(i,4),kf,lf,1))
     #                                          *8*alpac+
     #                           (dik(index(i,1),index(i,3),if,kf,1)*
     #                            djl(index(i,2),index(i,4),jf,lf,2)+
     #                            dil(index(i,1),index(i,4),if,lf,1)*
     #                            djk(index(i,2),index(i,3),jf,kf,2)+
     #                            dik(index(i,1),index(i,3),if,kf,2)*
     #                            djl(index(i,2),index(i,4),jf,lf,1)+
     #                            dil(index(i,1),index(i,4),if,lf,2)*
     #                            djk(index(i,2),index(i,3),jf,kf,1))
     #                                           *4*betac
  302                      continue
c..bhl
  311                   continue
  312                continue
  313             continue
  314          continue
c
      else
c
c        ----- [ij;ll]  [ii;kl]  [ij;ij] -----
c
c
               pt=0
               do 414 if=1,nfi
                  do 413 jf=1,nfj
                     do 412 kf=1,nfk
                        do 411 lf=1,nfl
                           spt=pt
                           do 401 i=1,nv
                              pt=pt+1
                              d2cc(pt)=d2cc(pt)+
     #                           dij(index(i,1),index(i,2),if,jf,1)*
     #                           dkl(index(i,3),index(i,4),kf,lf,1)*
     #                                          4*alp+
     #                          (dik(index(i,1),index(i,3),if,kf,1)*
     #                           djl(index(i,2),index(i,4),jf,lf,1)+
     #                           dil(index(i,1),index(i,4),if,lf,1)*
     #                           djk(index(i,2),index(i,3),jf,kf,1))
     #                                          *2*bet
  401                      continue
c..bhl
                           pt=spt
                           do 402 i=1,nv
                              pt=pt+1
                              d2ac(pt)=d2ac(pt)+
     #                          (dij(index(i,1),index(i,2),if,jf,1)*
     #                           dkl(index(i,3),index(i,4),kf,lf,2)+
     #                           dij(index(i,1),index(i,2),if,jf,2)*
     #                           dkl(index(i,3),index(i,4),kf,lf,1))
     #                                         *4*alpac+
     #                          (dik(index(i,1),index(i,3),if,kf,1)*
     #                           djl(index(i,2),index(i,4),jf,lf,2)+
     #                           dil(index(i,1),index(i,4),if,lf,1)*
     #                           djk(index(i,2),index(i,3),jf,kf,2)+
     #                           dik(index(i,1),index(i,3),if,kf,2)*
     #                           djl(index(i,2),index(i,4),jf,lf,1)+
     #                           dil(index(i,1),index(i,4),if,lf,2)*
     #                           djk(index(i,2),index(i,3),jf,kf,1))
     #                                          *2*betac
  402                      continue
c..bhl
  411                   continue
  412                continue
  413             continue
  414          continue
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
