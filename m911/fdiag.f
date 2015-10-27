*deck @(#)fdiag.f	5.1  11/6/94
      subroutine fdiag(v,d,smaln,smalo,itrdm,ioptci,state1,state2)
      implicit real*8 (a-h,o-z)
c
      integer ijunk(100)
      integer arr,refwlk,symorb,bmax,fword,fword2,orbfrm
      integer wptoin
      real*8 rjunk(20)
      real*8 d(nwksmx),v(nwksmx)
      equivalence(ijunk(21),rjunk)
      logical pagein
      character*4 itoc
c
      integer state1,state2
c
      common /dims/ nbf,nsym,norbs,nrows,nrows4,nwks,nwks2,nlevs
     *,             nrowoc,nrow4o,nwksoc,nlevoc,norboc,levfrm
     *,             nwksmx,nlwkmx,nuwkmx,bmax,nroots,orbfrm
      common /ints/   nmax,nmax2,ngroup,nblkoc,numij,symorb,intsrt
      common /optns/  iguess,irstrt,irooti,irootf,i34x
      common /symm/ jsm,jfsym,ifsym,maxsym(8),minsym(8),ismoff(8)
      common /diag/ rep,fzcore,eguess,eci,cnverg,sqcdif,czero
     *,             refwlk,mxiter,icnvg,iter,nroot
      common /tapes/itap20,itape5,itape6,itape8,itap12,itap03,itap04
     *,             itape3,itap05,itap06
      common /all/  val1,val2,val3,arr,itr1,itr2,ia,ja,itype,isegt
     *,             lvfrm1,nlwki,nlwkj,imax,imin
      common /page/ iword3,iword4,pagein,noffi,noffj,npass,nwordl,iword6
      common /rstrt/  iblck1,inxt,lowtri,ndvdit,iblock
      common/dry2/prtciv,navail,nrused(2)
c
  205 format(' nwkmx2=',i8)
  206 format(5x,'reading ',i2,' vector(s) from tape 12.')
  207 format(1x,' root=',i3,' prin ref=',i7,' ci energy=',g17.11)
c
      sqcdif=1.0d+00
      eguess=0.0d+00
      icnvg=-1
c pws
      npass=(nwks-1)/nwksmx+1
c   prepare for unit guess on lowest diagonal matrix element
      jj=0
      nwkmx2=intowp(nwksmx+nwksmx)
      nword=nwkmx2
c      write (itape6,205) nwkmx2
      fword=1
      if(iguess.le.0) then
         call lnkerr(' illegal path for dm ')
      end if
      if(itrdm.ne.0)then
         write(itape6,206) 2
         mci=2
      else
         write(itape6,206) 1
         mci=1
      endif
      nroot=2
      ndvdit=2
c
      if(itrdm.ne.0) then
c..bhl
                  ir=state1
                  call iosys('rewind "ci root '//itoc(ir)
     #                       //'"  on rwf',0,0,0,' ')
                  call iosys('read real "ci root '//itoc(ir)
     #                     //'" from rwf',nwksmx,v,0,' ')
c..bhl
                  ir=state2
                  call iosys('rewind "ci root '//itoc(ir)
     #                       //'" on rwf',0,0,0,' ')
                  call iosys('read real "ci root '//itoc(ir)
     #                       //'" from rwf',nwksmx,d,0,' ')

         tdot=0.0d0
      else
         tdot=1.0d0
c..bhl
         ir=nrused(1)
         call iosys('rewind "ci root '//itoc(ir)//'" on rwf',
     #               0,0,0,' ')
         call iosys('read real "ci root '//itoc(ir)//'" from rwf',
     #              nwksmx,v,0,' ')
c..bhl
c
         do 19001 i=1,nwks
            d(i)=v(i)
19001    continue
      endif
c
c     check normalization
      rn1=-1.0d0
      dot=0.d0
      rn2=-1.0d0
      do 250 i=1,nwks
         rn1=rn1+v(i)**2
         rn2=rn2+d(i)**2
         dot=dot+v(i)*d(i)
  250 continue
  270 continue
c
      if(nrused(1).eq.nrused(2)) tdot=1.0d0
      if((dabs(rn1).lt.smaln).and.(dabs(rn2).lt.smaln)
     $                       .and.(dabs(dot-tdot).lt.smalo))then
      else
         write(iout,1000)rn1,rn2,dot
         call lnkerr(' vector error in fdiag ')
      end if
c
      return
 1000 format(5x,' vectors are junk:rn1 rn2 dot ',3f12.6)
 1001 format(/5x,' in dm:root=',i3,' nwks=',i7)
      end
