*deck phig
      subroutine phig
      implicit real *8 (a-h,o-z)
      common/kar/u(30),usav(30),grad(30),z(30),zsav(30),bl(30),bu(30),v(
     130),vsav(30),gradu(30),nuflag(30),nbd(30)
      common /kchg/ gsav( 300,14),b( 300),bsav( 300),nersdim
      common /ab/ bsav2( 300),esav( 300)
      common/kvar/nu,nprob,ners,numax,nerspnj,phi,solves,ngrad,
     1netasr,gradsq,zsq
      common/d1/ yz,timesum,nbdsum,nindep1,rjs,iter
      if(ngrad.eq.2) go to 2
      nindep1=0
      do 40 i=1,nu
      if(nuflag(i).eq.0) nindep1=nindep1+1
   40 continue
      if(ngrad.eq.0) call fcn(v,esav,phi)
      if(ngrad.eq.-1) call fcn(v,bsav,phi)
      if(ngrad.eq.1) call der(v,gsav,bsav)
      if(iabs(ngrad)-1) 2,14,2
   14 do 15 k=1,nindep1
      grad(k)=0.d0
      do 15 j=1,ners
   15 grad(k)=grad(k)-bsav(j)*gsav(j,k)
    2 continue
      return
      end
