*deck @(#)unpck.f	5.1  11/6/94
      subroutine unpck(idstor,ia,nbits,number)
      integer shiftr
      integer and
      dimension ia(*), idstor(*)
*mdc*if cray
*      parameter (intl=64)
*mdc*else
      parameter (intl=32)
*mdc*endif
      msk=2**nbits-1
      npw=intl/nbits
      ipwu=npw
      nw=(number-1)/npw+1
      do 20 iw=1,nw
      if(iw.eq.nw) ipwu=number-(iw-1)*npw
      do 20 ipw=1,ipwu
*mdc*if cray
*   20 ia((iw-1)*npw+ipw)= and(shiftr(idstor(iw),(npw-ipw)*nbits),msk)
*mdc*elseif sun
   20 ia((iw-1)*npw+ipw)= and(shiftr(idstor(iw),(npw-ipw)*nbits),msk)
*mdc*else
*  20 ia((iw-1)*npw+ipw)=iand( ishft(idstor(iw),(ipw-npw)*nbits),msk)
*mdc*endif
      return
      end
