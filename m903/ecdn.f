*deck @(#)ecdn.f	5.1  11/6/94
      subroutine ecdn(ec,nnp,nkwks,dn,idn,wpti4,ndnlps,
     #                idnnwk,jdnnwk,jupnwk,offset,c,nwks,iwlkup,
     #                mini,maxi)
c
c***begin prologue     ecdn
c***date written       860915  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            form the product e(ij,ij)*c(j) for loops in the
c                      lower portion of the graph.
c***description
c
c***references
c***routines called
c***end prologue       ecdn
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      real*8 ec(nnp,nkwks)
c
c     ----- external arrays not modified -----
c           note that dn,idn are equivalenced in the calling routine.
c
      real*8 dn(4,ndnlps),c(nwks)
      integer idn(wpti4,ndnlps)
      integer iwlkup(jupnwk)
c
c     ----- external scalars not modified -----
c
      real*8 e
      integer nnp,nkwks,ndnlps,idnnwk,jdnnwk,offset,nwks
c
c     ----- find the range of upperwalks we can handle -----
c
      do 1 i=1,jupnwk
         if (iwlkup(i).ge.mini-1) then
            upmin=i
            go to 2
         end if
    1 continue
      return
c
    2 continue
      do 3 i=jupnwk,1,-1
         if (iwlkup(i).le.maxi-1) then
            upmax=i
            go to 4
         end if
    3 continue
      return
c
    4 continue
c
c
      do 10 bot=1,ndnlps
         ij=idn(3,bot)
         iwkmin=idn(1,bot)+1-(mini-1)*idnnwk
         jwkmin=idn(2,bot)+1-jdnnwk+offset
         e=dn(4,bot)
cdir$ ivdep
         do 5 upwk=upmin,upmax
            ec(ij,iwkmin+iwlkup(upwk)*idnnwk)=
     #               ec(ij,iwkmin+iwlkup(upwk)*idnnwk)+
     #                                 e*c(jwkmin+upwk*jdnnwk)
    5    continue
   10 continue
c
c
      return
      end
