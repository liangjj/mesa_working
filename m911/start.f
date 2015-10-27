*deck @(#)start.f	5.1  11/6/94
      block data start
      implicit real*8 (a-h,o-z)
      common /diag/ rep,fzcore,eguess,eci,cnverg,sqcdif,czero
     *,             refwlk,mxiter,icnvg,iter,nroot
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      data itap20,itape3,itape5,itape6,itape8,itap12/20,3,5,6,8,12/
      data itap03,itap04,itap05,itap06/103,104,105,104/
      data inxt,iter,ndvdit,iblck1/1,0,0,1/
c     nwkmx2=nwksmx*2
      end
