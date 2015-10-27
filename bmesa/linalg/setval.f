      subroutine setval (bndm,bndspn,bndpar,enbnd,ncst,cntm,cntspn
     1 ,cntpar,l0cnt,ldel,chnloc,jind,indxc,lch,ech,eng,en,rk,nsts,ncmax
     2 ,ntchn,lplsmx,neng)
      implicit integer(a-z)
      real *8enbnd, eng, en, rk, ediff, ech
      character *4 itoc
      common /io/ inp, iout
      dimension bndm(nsts), bndspn(nsts), bndpar(nsts), enbnd(nsts)
      dimension ncst(nsts), cntm(nsts), cntspn(nsts), cntpar(nsts)
      dimension l0cnt(nsts), chnloc(lplsmx,nsts), jind(ncmax,nsts)
      dimension eng(neng), en(ncmax,nsts), rk(ncmax,nsts)
      dimension indxc(ntchn), lch(ntchn), ech(ntchn)
*
*
      call iosys ('read integer "trgt m" from rwf',nsts,bndm,0,0)
      call iosys ('read integer "trgt spin" from rwf',nsts,bndspn,0,0)
      call iosys ('read integer "trgt parity" from rwf',nsts,bndpar,0,0)
      call iosys ('read real "trgt energies" from rwf',nsts,enbnd,0,0)
      call iosys ('read integer "no. cntm states" from rwf',nsts,ncst,0,
     1 0)
      call iosys ('read integer "cntm m" from rwf',nsts,cntm,0,0)
      call iosys ('read integer "cntm spin" from rwf',nsts,cntspn,0,0)
      call iosys ('read integer "cntm parity" from rwf',nsts,cntpar,0,0)
      call iosys ('read integer "cntm initial l" from rwf',nsts,l0cnt,0,
     1 0)
      call iosys ('read integer "cntm l spacing" from rwf',1,ldel,0,0)
      call iosys ('read real "scat. eng." from rwf',neng,eng,0,0)
      ic=0
      do 10 i=1,nsts
      call iosys ('read integer "chnl array'//itoc(i)//'" from rwf'
     1 ,lplsmx,chnloc(1,i),0,0)
      l=l0cnt(i)-ldel
      nn=ncst(i)
      do 10 j=1,nn
      ic=ic+1
      l=l+ldel
      lch(ic)=l
      indxc(ic)=i
      jind(j,i)=l
   10 continue
      write (iout,80)
      write (iout,90)
      do 20 is=1,nsts
   20 write (iout,100) is,bndm(is),bndspn(is),bndpar(is),enbnd(is)
      write (iout,110)
      write (iout,120)
      do 40 is=1,nsts
      nsch=ncst(is)
      write (iout,130) is,nsch,cntm(is),cntspn(is),cntpar(is)
      write (iout,140) (jind(j,is),j=1,nsch)
      do 30 i=1,nsch
   30 jind(i,is)=jind(i,is)+1
   40 continue
      call iosys ('write integer "chnl l values" to rwf',ntchn,lch,0,0)
      call iosys ('write integer "chnl state indx" to rwf',ntchn,indxc,0
     1 ,0)
*
*         now write out the energy dependent information
*
      do 70 i=1,neng
      ic=0
      do 60 j=1,nsts
      ediff=enbnd(1)-enbnd(j)
      ediff=2.d+00*ediff
      nc=ncst(j)
      do 50 k=1,nc
      ic=ic+1
      en(k,j)=eng(i)-ediff
      ech(ic)=en(k,j)
      rk(k,j)=sqrt(abs(en(k,j)))
   50 continue
   60 continue
      call iosys ('write real "chnl energy'//itoc(i)//'" to rwf',ntchn
     1 ,ech,0,0)
   70 continue
      return
c
   80 format (/,20x'target state information')
   90 format (//,5x,'state',2x,'m val',2x,'spin',2x,'parity',4x,'energy
     1 ')
  100 format (/,5x,i3,4x,i3,3x,i3,5x,i3,2x,f12.6)
  110 format (/,20x,'continuum state information')
  120 format (//,5x,'state',2x,'no. l val',2x,'m val',4x,'spin',5x,'pari
     1ty')
  130 format (/,5x,i3,6x,i3,6x,i3,5x,i3,7x,i3)
  140 format ((/,5x,'l values',2x,10i3))
      end
