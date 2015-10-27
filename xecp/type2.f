*deck %W%  %G%
      subroutine type2(lm,nt,np,alpp,dp,sp,exx,maxi2,nprima,
     $                  nprimb,nx,ny,nz,vca,vcb,lamax,lbmax,
     $                  igbegn,igend,jgbegn,jgend,istart,iend,
     $                  jstart,jend,qq,atheta,aphik,btheta,bphik,
     $                  xab,yab,zab,q2,anga,angb,
     $                  acoef2,len2,ptcf2,ltop)
c***begin prologue     %M%
c***date written       920615  
c***revision date      %G%      
c
c***keywords           ecp, 
c***author             martin,richard (lanl) 
c***source             %W%   %G%
c***purpose            computes type2 ecp integrals 
c***description        computes pseudopotential integrals of type 2.
c                      lm     maximum l value of the projector.
c                      nt     an array which specifies the number of terms
c                             for each l component
c                      np     the exponent of r for each l component
c                      alpp   the gaussian exponents for each l component
c                      dp     the coefficient for each l component
c                      sp     array which contains the integrals for this
c                             shell block
c                      exx    list of primitive exponents
c                      nprima the number of primitives in shell a
c                      nprimb the number of primitives in shell b
c                      maxi2 the number of integrals in this shell combination
c                               e.g. p x p =9*nprima*nprimb integrals
c    
c
c***references         l.e. mcmurchie and e.r. davidson,
c                         j.comp.phys., 44,289(1981).
c                      r.l.martin, unpublished notes.
c***routines called
c
c***end prologue       %M%
      implicit integer(a-z)
c     computes pseudopotential integrals type 2.......
c     lm is the max l value + 1 for the potential.
c     nt is the no. terms for each l.
c     np contains the value of n for each term.
c     alpp contains the value of alpha for each term.
c     dp contains the coefficient of each term.
c     sp will contain the calculated integrals.
c     qq is scratch storage for the radial integralls.
c     maxi2 is the number of integrals in the largest shell block.
c     nprima*nprimb is the number of primitives in the largest primitive block.
c
c     ----- arguments unchanged -----
      integer lm,nt(*),np(*)
      integer maxi2,nprima,nprimb,nx,ny,nz
      integer lamax,lbmax
      integer igbegn,igend,jgbegn,jgend
      integer istart,iend,jstart,jend
      integer ptcf2(0:ltop,0:ltop,0:ltop,0:lm)
      real*8 alpp(*),dp(*),exx(*)
      real*8 vca(3),vcb(3)
      complex*16 acoef2(len2)
c     ----- arguments returned -----
      real*8 sp(nprima*nprimb,maxi2)
c     ----- scratch -----
      real*8 qq(nprima*nprimb,0:ltop,0:ltop,0:ltop)
      real*8 atheta(nprima*nprimb),aphik(nprima*nprimb)
      real*8 btheta(nprima*nprimb),bphik(nprima*nprimb)
      real*8 xab(0:ltop+ltop),yab(0:ltop+ltop),zab(0:ltop+ltop)
      real*8 q2(nprima*nprimb)
      complex*16 anga(0:lm,0:lm+lamax,nprima)
      complex*16 angb(0:lm,0:lm+lbmax,nprimb)
c     ----- local variables -----
      real*8 zero,one,two
      real*8 eps1
      real*8 ca2,cb2,ca,cb,alpa,alpb,xk,yk,zk,rk
      real*8 pab,angp
      complex*16 t1,t2
c     ----- common -----
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/io/in,iout
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
      parameter (eps1=1.d-15)
c
c     ----- find distances between sites c(ecp center),a, and b.
      ca2=zero
      cb2=zero
      do 10 i=1,3
         ca2=ca2+vca(i)*vca(i)
         cb2=cb2+vcb(i)*vcb(i)
   10 continue
      ca=sqrt(ca2)
      cb=sqrt(cb2)
c
c     ----- generate the angular arguments for all primitives.
      prim=0
      do 15 igauss=igbegn,igend
         prim=prim+1
         alpa=exx(igauss)
         xk=-two*alpa*vca(1)
         yk=-two*alpa*vca(2)
         zk=-two*alpa*vca(3)
         rk=sqrt(xk*xk+yk*yk+zk*zk)
         if (rk.eq.zero) then
            xk=zero
            yk=zero
            zk=one
            atheta(prim)=zero
            aphik(prim)=zero
         else
            xk=xk/rk
            yk=yk/rk
            zk=zk/rk
            atheta(prim)=acos(zk)
            aphik(prim)=atan(yk/xk)
         endif
   15 continue
      do 20 jgauss=jgbegn,jgend
         alpb=exx(jgauss)
         prim=prim+1
         xk=-two*alpb*vcb(1)
         yk=-two*alpb*vcb(2)
         zk=-two*alpb*vcb(3)
         rk=sqrt(xk*xk+yk*yk+zk*zk)
         if (rk.eq.zero) then
            xk=zero
            yk=zero
            zk=one
            btheta(prim)=zero
            bphik(prim)=zero
         else
            xk=xk/rk
            yk=yk/rk
            zk=zk/rk
            btheta(prim)=acos(zk)
            bphik(prim)=atan(yk/xk)
         endif
   20 continue
c
c     ----- sum over angular momentum projectors -----
      do 400 l=0,lm
c
c        ----- compute basic radial integrals needed for this shell
         call qqrad (l,lamax,lbmax,np(ndone+1),alpp(ndone+1),
     $               dp(ndone+1),nt(l),igbegn,igend,jgbegn,jgend,exx,
     $               ca,cb,qq,maxi2,nprima*nprimb)
c******
c
c
c***** from here down hacked on
         intc=0
c        ----- sum over shell block components
         do 310 ii=istart,iend
            na=nx(ii)
            la=ny(ii)
            ma=nz(ii)
            do 300 iii=jstart,jend
               nb=nx(iii)
               lb=ny(iii)
               mb=nz(iii)
               call facab(na,nb,vca(1),vcb(1),xab)
               call facab(la,lb,vca(2),vcb(2),yab)
               call facab(ma,mb,vca(3),vcb(3),zab)
               intc=intc+1
c              ----- check for special cases -----
c                    test parity, and make sure the function can have
c                    a component of the projector.
               if ((ca.eq.zero).and.(cb.eq.zero)) then
c                 both function centers coincident with projector 
                  if(mod(l-lamax,2).ne.0) goto 300
                  if(mod(l-lbmax,2).ne.0) goto 300
                  if((l.gt.lamax).or.(l.gt.lbmax)) goto 300
               else if (ca.eq.zero) then
                  if(mod(l-lamax,2).ne.0) goto 300
                  if(l.gt.lamax) goto 300
               else if (cb.eq.zero) then
                  if(mod(l-lbmax,2).ne.0) goto 300
                  if(l.gt.lbmax) goto 300
               endif
c
c              ----- compute basic angular integrals -----
               do 100 prim=1,nprima*nprimb
                  q2(prim)=zero
  100          continue
               do 290 a=0,na
                  do 290 b=0,la
                     do 290 c=0,ma
                        alo=max(l-a-b-c,0)
                        if(mod(l+a+b+c,2).eq.0) then
                           if(mod(alo,2).ne.0) alo=alo+1
                        else
                           if(mod(alo,2).eq.0) alo=alo+1
                        endif
c                       ----- get the type2 angular integral for this
c                             triad of a,b,c
                        ind=ptcf2(a,b,c,l)+1
                        do 210 prim=1,nprima
                           call omega2(a,b,c,l,atheta(prim),aphik(prim),
     $                                 anga(0,0,prim),ltop,lm,
     $                                 acoef2(ind))
  210                   continue
                        do 280 d=0,nb
                           do 280 e=0,lb
                              do 280 f=0,mb
                                 pab=xab(a+d)*yab(b+e)*zab(c+f)
                                 if (pab.eq.zero) go to 280
                                 blo=max(l-d-e-f,0)
                                 if(mod(l+d+e+f,2).eq.0) then
                                    if(mod(blo,2).ne.0) blo=blo+1
                                 else
                                    if(mod(blo,2).eq.0) blo=blo+1
                                 endif
c                                ----- get the type2 angular integral
c                                      for the d,e,f triad
                        ind=ptcf2(d,e,f,l)+1
                        do 220 prim=1,nprimb
                           call omega2(d,e,f,l,btheta(prim),
     $                                 bphik(prim),angb(0,0,prim),
     $                                 ltop,lm,acoef2(ind))
 220                    continue
c
c
                        n=a+b+c+d+e+f
                        do 250 lama=alo,l+a+b+c,2
                        do 250 lamb=blo,l+d+e+f,2
                           angp=zero
                           prim=0
                           do 240 prima=1,nprima
                           do 240 primb=1,nprimb
                              prim=prim+1
c                             ----- do the m=0 term first
                              angp=real(anga(0,lama,prima))
     $                             *real(angb(0,lamb,primb))
                              do 230 m=1,l
                                 t1=anga(m,lama,prima)
                                 t2=angb(m,lamb,primb)
                                 angp=angp+t1*conjg(t2)+conjg(t1)*t2
  230                         continue
                              q2(prim)=q2(prim)+pab*angp
     $                                 *qq(prim,n,lama,lamb) 
  240                      continue
  250                   continue
c
c
  280             continue
  290          continue
c
c              ----- accumulate integrals into sp
               do 305 prim=1,nprima*nprimb
                  sp(prim,intc)=sp(prim,intc) 
     $                          +q2(prim)
  305          continue
c
c           end summation loop
  300       continue
  310    continue
c        end shell loop.
c
  400 continue
c     end projector loop.
c
      return
      end
