*deck @(#)pseud2.f	5.1  11/6/94
      subroutine pseud2(lm,ntcnt,np,alpp,dp,sp,exx,qq,maxi2,maxp2,
     $                  nx,ny,nz)
      implicit real*8(a-h,o-z)
c
c     computes pseudopotential integrals type 2.......
c     lm is the max l value + 1 for the potential.
c     nt is the no. terms for each l.
c     np contains the value of n for each term.
c     alpp contains the value of alpha for each term.
c     dp contains the coefficient of each term.
c     sp will contain the calculated integrals.
c     qq is scratch storage for the radial integralls.
c     maxi2 is the number of integrals in the largest shell block.
c     maxp2 is the number of primitives in the largest primitive block.
c
      common/dist/cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      common/center/xa,ya,za,xb,yb,zb,xc,yc,zc
      common/angmax/lamax,lbmax
      common/limit/istart,jstart,iend,jend
      common/const/zero,one,two,three,four,five,six,ten
      common/prims/igbegn,igend,jgbegn,jgend
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/io/in,iout
      dimension qq(11,9,9,maxp2)
      dimension np(1), alpp(1), dp(1), ntcnt(1)
      dimension q2sum(11,9,9), anga(18,7,9), angb(18,7,9)
      dimension exx(1),nx(1),ny(1),nz(1), sp(maxp2,maxi2)
      dimension binom(21), ind(6)
      data binom /1., 1.,1., 1.,2.,1., 1.,3.,3.,1., 1.,4.,6.,4.,1.,
     $            1.,5.,10.,10.,5.,1./
      data ind /1,2,4,7,11,16/
      data eps1 /1.d-15/
      save binom,ind,eps1
c
      lm1=lm+1
      if (lm1.eq.1) go to 270
c
      if (ca.eq.zero) then
         xka=zero
         yka=zero
         zka=one
      else 
         xka=-cax/ca
         yka=-cay/ca
         zka=-caz/ca
      endif
      if (cb.eq.zero) then
         xkb=zero
         ykb=zero
         zkb=one
      else
         xkb=-cbx/cb
         ykb=-cby/cb
         zkb=-cbz/cb
      endif
c
c
c     sum over terms in the projector.
      llo=1
      lhi=lm
      inc=1
      ndone=0
      do 240 l=llo,lhi,inc
c
c        compute basic radial integrals needed for this projector
c        and block of primitives.
         lmalo=max(l-lamax,0)+1
         lmahi=l+lamax-1
         lmblo=max(l-lbmax,0)+1
         lmbhi=l+lbmax-1
         nhi=lamax+lbmax-1
         call qqrad (l,lmalo,lmahi,lmblo,lmbhi,nhi,np(ndone+1),
     $           alpp(ndone+1),dp(ndone+1),ntcnt(l),exx,qq,maxi2,maxp2)
         ndone=ndone+ntcnt(l)
c
c        sum over shell blocks.
         mhi=l+l-1
         intc=0
         do 230 ii=istart,iend
            na1=nx(ii)+1
            la1=ny(ii)+1
            ma1=nz(ii)+1
            naind=ind(na1)
            laind=ind(la1)
            maind=ind(ma1)
            ltota=na1+la1+ma1-2
            lmalo=max(l-ltota,0)+1
            lmahi=l+ltota-1
            call ang2 (na1-1,la1-1,ma1-1,l-1,xka,yka,zka,anga)
            do 230 iii=jstart,jend
               nb1=nx(iii)+1
               lb1=ny(iii)+1
               mb1=nz(iii)+1
               nbind=ind(nb1)
               lbind=ind(lb1)
               mbind=ind(mb1)
               ltotb=nb1+lb1+mb1-2
               lmblo=max(l-ltotb,0)+1
               lmbhi=l+ltotb-1
               nhi=ltota+ltotb-1
               intc=intc+1
c
c              check for special cases (two or more centers coincident.)
               if ((ca.eq.zero).and.(cb.eq.zero)) then
                  ijlhi=min(lm,ltota)
                  ijlhi=min(ijlhi,ltotb)
                  ijllo=mod(ltota-1,2)+1
                  ijllob=mod(ltotb-1,2)+1
                  if (ijllo.ne.ijllob) go to 230
                  if (ijllo.gt.ijlhi) go to 230
                  if ((l.lt.ijllo).or.(l.gt.ijlhi)) go to 230
                  if (mod(l-ijllo,2).ne.0) go to 230
                  ijinc=2
               else if (ca.eq.zero) then
                  ijlhi=min(lm,ltota)
                  ijllo=mod(ltota-1,2)+1
                  if (ijllo.gt.ijlhi) go to 230
                  if ((l.lt.ijllo).or.(l.gt.ijlhi)) go to 230
                  if (mod(l-ijllo,2).ne.0) go to 230
                  ijinc=2
               else if (cb.eq.zero) then
                  ijlhi=min(lm,ltotb)
                  ijllo=mod(ltotb-1,2)+1
                  if (ijllo.gt.ijlhi) go to 230
                  ijinc=2
               endif
               call ang2 (nb1-1,lb1-1,mb1-1,l-1,xkb,ykb,zkb,angb)
c
c     compute basic angular integrals.
               do 110 lama=lmalo,lmahi
                  do 110 lamb=lmblo,lmbhi
                     do 110 n=1,nhi
                        q2sum(n,lama,lamb)=zero
  110          continue
               pmax=zero
               pab=one
               ijka=0
               do 210 ia=1,na1
                  pab1=pab
                  if (ia.ne.na1) then
                     pab1=pab*binom(naind+ia-1)*(cax**(na1-ia))
                  endif
                  do 210 ja=1,la1
                     pab2=pab1
                     if (ja.ne.la1) then
                        pab2=pab1*binom(laind+ja-1)*(cay**(la1-ja))
                     endif
                     do 210 ka=1,ma1
                        pab3=pab2
                        if (ka.ne.ma1) then
                           pab3=pab2*binom(maind+ka-1)*(caz**(ma1-ka))
                        endif
                        ijka=ijka+1
                        iajaka=ia+ja+ka-3
                        lamahi=iajaka+l
                        lamalo=max(l-1-iajaka,0)+1
                        if (mod(lamahi-lamalo,2).ne.0) lamalo=lamalo+1
                        ijkb=0
                        do 200 ib=1,nb1
                           pab4=pab3
                           if (ib.ne.nb1) then
                              pab4=pab3*binom(nbind+ib-1)
     $                            *(cbx**(nb1-ib))
                           endif
                           do 200 jb=1,lb1
                              pab5=pab4
                              if (jb.ne.lb1) then
                                 pab5=pab4*binom(lbind+jb-1)
     $                               *(cby**(lb1-jb))
                              endif
                              do 200 kb=1,mb1
                                 pab6=pab5
                                 if (kb.ne.mb1) then
                                    pab6=pab5*binom(mbind+kb-1)
     $                                  *(cbz**(mb1-kb))
                                 endif
                                 ijkb=ijkb+1
                                 if (pab6.eq.zero) go to 200
                                 ibjbkb=ib+jb+kb-3
                                 lambhi=ibjbkb+l
                                 lamblo=max(l-1-ibjbkb,0)+1
                                 if (mod(lambhi-lamblo,2).ne.0) 
     $                              lamblo=lamblo+1
                                 do 190 lama=lamalo,lamahi,2
                                    do 190 lamb=lamblo,lambhi,2
                                       angp=zero
                                       do 180 m=1,mhi
                                          angp=angp
     $                                        +anga(ijka,m,lama)
     $                                        *angb(ijkb,m,lamb)
  180                                  continue
                                       n=iajaka+ibjbkb+1
                                       prang=pab6*angp
                                       pmax=max(abs(prang),pmax)
                                       q2sum(n,lama,lamb)=
     $                                         q2sum(n,lama,lamb)+prang
  190                            continue
  200                   continue
  210          continue
c
c              combine radial and angular parts.
               do 220 lama=lmalo,lmahi
               do 220 lamb=lmblo,lmbhi
               do 220 n=1,nhi
                  if (q2sum(n,lama,lamb).ne.zero) then
                     do 215 iprim=1,maxp2
                        sp(iprim,intc)=sp(iprim,intc) 
     $                                +qq(n,lama,lamb,iprim)
     $                                *q2sum(n,lama,lamb)
  215                continue
                  endif
  220          continue
  230    continue
c        end shell loop.
  240 continue
c
c        end projector loop.
  270 continue
c
c
      return
      end
