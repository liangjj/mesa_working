*deck @(#)pck.f	5.1  11/6/94
      subroutine pck(idstor,ia,nbits,number)
      integer shiftl
      integer or
      dimension ia(*), idstor(*)
*mdc*if cray
*      parameter (intl=64)
*mdc*else
      parameter (intl=32)
*mdc*endif
      npw=intl/nbits
      ipwu=npw
      nw=(number-1)/npw+1
      do 20 iw=1,nw
         idstor(iw)=0
         if(iw.eq.nw) ipwu=number-(iw-1)*npw
         do 20 ipw=1,ipwu
*mdc*if cray
* 20       idstor(iw)= or(shiftl(idstor(iw),nbits),ia((iw-1)*npw+ipw))
*      if(ipwu.ne.npw) idstor(nw)=shiftl(idstor(nw),(npw-ipwu)*nbits)
*mdc*elseif sun
   20 idstor(iw)= or(shiftl(idstor(iw),nbits),ia((iw-1)*npw+ipw))
      if(ipwu.ne.npw) idstor(nw)=shiftl(idstor(nw),(npw-ipwu)*nbits)
*mdc*else
*  20 idstor(iw)=ior( ishft(idstor(iw),nbits),ia((iw-1)*npw+ipw))
*     if(ipwu.ne.npw) idstor(nw)= ishft(idstor(nw),(npw-ipwu)*nbits)
*mdc*endif
      return
      end
