*deck @(#)hmult.f	5.1  11/6/94
      subroutine hmult(b,mdim,mvec,t,hss,nrow,incore)
      implicit real*8(a-h,o-z)
      real*8 b(mdim,mvec),hss(*),t(mdim,mvec)
c
      if(incore.eq.1) then
         call ebc(t,hss,b,mdim,mdim,mvec)
      else
         ix=0
         do 1 i=1,mdim,nrow
            nr=min(nrow,mdim-i+1)
            lread=mdim*nr
            ix=ix+lread
            call mxma(hss,1,mdim,b,1,mdim,t(i,1),1,mdim,nr,mdim,mdim)
    1    continue
      endif
c
      return
      end
