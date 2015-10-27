*deck @(#)sup.f	5.1  11/6/94
      subroutine sup(gec,nnp,nkwks,up,iup,wpti4,nuplps,
     #               idnnwk,jdnnwk,offset,s,nwks,iwlkdn,mini,maxi)
c
c***begin prologue     sup
c***date written       860915  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            form s(j)=s(j)+gec(ij,i)*e(ij,ij) for loops in the
c                      upper portion of the graph.
c***description
c
c***references
c***routines called
c***end prologue       sup
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      real*8 s(nwks)
c
c     ----- external arrays not modified -----
c
      real*8 up(4,nuplps),gec(nnp,nkwks)
      integer iup(wpti4,nuplps)
      integer iwlkdn(jdnnwk)
c
c     ----- external scalars not modified -----
c
      real*8 e
      integer nnp,nkwks,nuplps,idnnwk,jdnnwk,offset,nwks
c
      do 10 top=1,nuplps
         if (iup(1,top).lt.mini-1.or.iup(1,top).gt.maxi-1) go to 10
         ij=iup(3,top)
         iwkmin=(iup(1,top)-mini+1)*idnnwk+1
         jwkmin=iup(2,top)*jdnnwk+offset
         e=up(4,top)
         do 5 dnwk=1,jdnnwk
            s(jwkmin+dnwk)=s(jwkmin+dnwk)+gec(ij,iwkmin+iwlkdn(dnwk))*e
    5    continue
   10 continue
c
c
      return
      end
