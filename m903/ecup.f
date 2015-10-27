*deck @(#)ecup.f	5.1  11/6/94
      subroutine ecup(ec,nnp,nkwks,up,iup,wpti4,nuplps,
     #                idnnwk,jdnnwk,offset,c,nwks,iwlkdn,mini,maxi)
c
c***begin prologue     ecup
c***date written       860915  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source
c***purpose            form the product e(ij,ij)*c(j) for loops in the
c                      upper portion of the graph.
c***description
c
c***references
c***routines called
c***end prologue       ecup
c
      implicit integer (a-z)
c
c     ----- external arrays modified -----
c
      real*8 ec(nnp,nkwks)
c
c     ----- external arrays not modified -----
c           note that up,iup are equivalenced in the calling routine.
c
      real*8 up(4,nuplps),c(nwks)
      integer iup(wpti4,nuplps)
      integer iwlkdn(jdnnwk)
c
c     ----- external scalars not modified -----
c
      integer nnp,nkwks,nuplps,idnnwk,jdnnwk,offset,nwks
c
c     ----- internal scalars not modified -----
c
      real*8 e
c
c
      do 10 top=1,nuplps
         if (iup(1,top).lt.mini-1.or.iup(1,top).gt.maxi-1) go to 10
         ij=iup(3,top)
         iwkmin=(iup(1,top)-mini+1)*idnnwk+1
         jwkmin=iup(2,top)*jdnnwk+offset
         e=up(4,top)
cdir$ ivdep
         do 5 dnwk=1,jdnnwk
            ec(ij,iwkmin+iwlkdn(dnwk))=ec(ij,iwkmin+iwlkdn(dnwk))+
     #                                 e*c(jwkmin+dnwk)
    5    continue
   10 continue
c
c
      return
      end
