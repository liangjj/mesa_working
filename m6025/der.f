*deck der
      subroutine der(vv,gsav,bsav)
      implicit real*8 ( a - h, o - z )
c
      common/kar/u1(300),nuflag(30),nbd(30)
      common/kvar/nu,nprob,ners,numax,nerspnu,phi,solves,ngrad,
     1netasr,gradsq,zsq
      common /ab/ bsav2( 300),esav( 300)
      dimension vv(3),gsav( 300,14),bsav( 300)
     x,del(30),vinc(30),vder(30)
      logical solves
      data delta/1.d-06/
      call fcn(vv,bsav,phi)
      do 1 inc=1,nu
      if (nuflag(inc) .eq. 1)   go to 1
      del(inc)=vv(inc)*delta
      if(abs(del(inc)).lt.1.d-20) del(inc)=1.d-04*delta
      vinc(inc)=vv(inc)+del(inc)
    1 continue
      nind=0
      do 10 ider=1,nu
      if(nuflag(ider).eq.1) go to 10
      nind=nind+1
      do 20 ivd=1,nu
      vder(ivd)=vv(ivd)
      if(ivd.eq.ider) vder(ivd)=vinc(ivd)
   20 continue
      call fcn(vder,bsav2,pp)
      do 30 lsav=1,ners
   30 gsav(lsav,nind)=(bsav(lsav)-bsav2(lsav))/del(ider)
   10 continue
      return
      end





