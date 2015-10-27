*deck @(#)hf2pdm.f	1.1  11/30/90
      subroutine hf2pdm(d2,dij,dkl,dik,djl,dil,djk,nprimi,nprimj,
     #                  nprimk,npriml,nfi,nfj,nfk,nfl,ishell,jshell,
     #                  kshell,lshell,lend2,index,nv,len)
c
c***begin prologue     hf2pdm
c***date written       851115   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             saxe, paul    (lanl)
c***source             @(#)hf2pdm.f	1.1   11/30/90
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
      real*8 dij(nprimi,nprimj,nfi,nfj),dkl(nprimk,npriml,nfk,nfl)
      real*8 dik(nprimi,nprimk,nfi,nfk),djl(nprimj,npriml,nfj,nfl)
      real*8 dil(nprimi,npriml,nfi,nfl),djk(nprimj,nprimk,nfj,nfk)
      real*8 d2(lend2)
      integer index(len,6)
c
      ijsh=ishell*(ishell-1)/2+jshell
      klsh=kshell*(kshell-1)/2+lshell
c
      if (ishell.eq.jshell.and.ishell.eq.kshell.and.ishell.eq.lshell)
     #                                                          then
c
c        ----- [ii;ii] -----
c
         pt=0
         do 114 if=1,nfi
            do 113 jf=1,nfj
               do 112 kf=1,nfk
                  do 111 lf=1,nfl
                     do 101 i=1,nv
                        pt=pt+1
                        d2(pt)=2*dij(index(i,1),index(i,2),if,jf)*
     #                           dkl(index(i,3),index(i,4),kf,lf)-
     #                           dik(index(i,1),index(i,3),if,kf)*
     #                           djl(index(i,2),index(i,4),jf,lf)
  101                continue
  111             continue
  112          continue
  113       continue
  114    continue
c
      else if (ishell.eq.jshell.and.kshell.eq.lshell) then
c
c        ----- [ii;jj] -----
c
         pt=0
         do 214 if=1,nfi
            do 213 jf=1,nfj
               do 212 kf=1,nfk
                  do 211 lf=1,nfl
                     do 201 i=1,nv
                        pt=pt+1
                        d2(pt)=4*dij(index(i,1),index(i,2),if,jf)*
     #                           dkl(index(i,3),index(i,4),kf,lf)-
     #                         2*dik(index(i,1),index(i,3),if,kf)*
     #                           djl(index(i,2),index(i,4),jf,lf)
  201                continue
  211             continue
  212          continue
  213       continue
  214    continue
c
      else if (ishell.ne.jshell.and.kshell.ne.lshell.and.
     #         ijsh.ne.klsh) then
c
c        ----- [ij;kl] ------
c
         pt=0
         do 314 if=1,nfi
            do 313 jf=1,nfj
               do 312 kf=1,nfk
                  do 311 lf=1,nfl
                     do 301 i=1,nv
                        pt=pt+1
                        d2(pt)=16*dij(index(i,1),index(i,2),if,jf)*
     #                            dkl(index(i,3),index(i,4),kf,lf)-
     #                          4*dik(index(i,1),index(i,3),if,kf)*
     #                            djl(index(i,2),index(i,4),jf,lf)-
     #                          4*dil(index(i,1),index(i,4),if,lf)*
     #                            djk(index(i,2),index(i,3),jf,kf)
  301                continue
  311             continue
  312          continue
  313       continue
  314    continue
c
      else
c
c        ----- [ij;ll]  [ii;kl]  [ij;ij] -----
c
         pt=0
         do 414 if=1,nfi
            do 413 jf=1,nfj
               do 412 kf=1,nfk
                  do 411 lf=1,nfl
                     do 401 i=1,nv
                        pt=pt+1
                        d2(pt)=8*dij(index(i,1),index(i,2),if,jf)*
     #                           dkl(index(i,3),index(i,4),kf,lf)-
     #                         2*dik(index(i,1),index(i,3),if,kf)*
     #                           djl(index(i,2),index(i,4),jf,lf)-
     #                         2*dil(index(i,1),index(i,4),if,lf)*
     #                           djk(index(i,2),index(i,3),jf,kf)
  401                continue
  411             continue
  412          continue
  413       continue
  414    continue
      end if
c
c
      return
      end
