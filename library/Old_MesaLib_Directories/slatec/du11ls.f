*deck du11ls
      subroutine du11ls (a, mda, m, n, ub, db, mode, np, krank, ksure,
     +   h, w, eb, ic, ir)
c***begin prologue  du11ls
c***subsidiary
c***purpose  subsidiary to dllsia
c***library   slatec
c***type      double precision (u11ls-s, du11ls-d)
c***author  (unknown)
c***description
c
c   **** double precision version of u11ls ****
c
c       this routine performs a qr factorization of a
c       using householder transformations. row and
c       column pivots are chosen to reduce the growth
c       of round-off and to help detect possible rank
c       deficiency.
c
c***see also  dllsia
c***routines called  daxpy, ddot, dnrm2, dscal, dswap, idamax, iswap,
c                    xermsg
c***revision history  (yymmdd)
c   810801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c***end prologue  du11ls
      implicit double precision (a-h,o-z)
      double precision ddot,dnrm2
      dimension a(mda,*),ub(*),db(*),h(*),w(*),eb(*)
      integer ic(*),ir(*)
c
c        initialization
c
c***first executable statement  du11ls
      j=0
      krank=n
      do 10 i=1,n
      ic(i)=i
   10 continue
      do 12 i=1,m
      ir(i)=i
   12 continue
c
c        determine rel and abs error vectors
c
c
c
c        calculate col length
c
      do 30 i=1,n
      h(i)=dnrm2(m,a(1,i),1)
      w(i)=h(i)
   30 continue
c
c         initialize error bounds
c
      do  40 i=1,n
      eb(i)=max(db(i),ub(i)*h(i))
      ub(i)=eb(i)
      db(i)=0.0d0
   40 continue
c
c          discard self dependent columns
c
      i=1
   50 if(eb(i).ge.h(i)) go to 60
      if(i.eq.krank) go to 70
      i=i+1
      go to 50
c
c          matrix reduction
c
   60 continue
      kk=krank
      krank=krank-1
      if(mode.eq.0) return
      if(i.gt.np) go to  64
      call xermsg ('slatec', 'du11ls',
     +   'first np columns are linearly dependent', 8, 0)
      krank=i-1
      return
   64 continue
      if(i.gt.krank) go to 70
      call dswap(1,eb(i),1,eb(kk),1)
      call dswap(1,ub(i),1,ub(kk),1)
      call dswap(1,w(i),1,w(kk),1)
      call dswap(1,h(i),1,h(kk),1)
      call iswap(1,ic(i),1,ic(kk),1)
      call dswap(m,a(1,i),1,a(1,kk),1)
      go to 50
c
c           test for zero rank
c
   70 if(krank.gt.0) go to 80
      krank=0
      ksure=0
      return
   80 continue
c
c        m a i n    l o o p
c
  110 continue
      j=j+1
      jp1=j+1
      jm1=j-1
      kz=krank
      if(j.le.np) kz=j
c
c        each col has mm=m-j+1 components
c
      mm=m-j+1
c
c         ub determines column pivot
c
  115 imin=j
      if(h(j).eq.0.d0) go to 170
      rmin=ub(j)/h(j)
      do 120 i=j,kz
      if(ub(i).ge.h(i)*rmin) go to 120
      rmin=ub(i)/h(i)
      imin=i
  120 continue
c
c     test for rank deficiency
c
      if(rmin.lt.1.0d0) go to 200
      tt=(eb(imin)+abs(db(imin)))/h(imin)
      if(tt.ge.1.0d0) go to 170
c     compute exact ub
      do 125 i=1,jm1
      w(i)=a(i,imin)
  125 continue
      l=jm1
  130 w(l)=w(l)/a(l,l)
      if(l.eq.1) go to 150
      lm1=l-1
      do 140 i=l,jm1
      w(lm1)=w(lm1)-a(lm1,i)*w(i)
  140 continue
      l=lm1
      go to 130
  150 tt=eb(imin)
      do 160 i=1,jm1
      tt=tt+abs(w(i))*eb(i)
  160 continue
      ub(imin)=tt
      if(ub(imin)/h(imin).ge.1.0d0) go to 170
      go to 200
c
c        matrix reduction
c
  170 continue
      kk=krank
      krank=krank-1
      kz=krank
      if(mode.eq.0) return
      if(j.gt.np) go to 172
      call xermsg ('slatec', 'du11ls',
     +   'first np columns are linearly dependent', 8, 0)
      krank=j-1
      return
  172 continue
      if(imin.gt.krank) go to 180
      call iswap(1,ic(imin),1,ic(kk),1)
      call dswap(m,a(1,imin),1,a(1,kk),1)
      call dswap(1,eb(imin),1,eb(kk),1)
      call dswap(1,ub(imin),1,ub(kk),1)
      call dswap(1,db(imin),1,db(kk),1)
      call dswap(1,w(imin),1,w(kk),1)
      call dswap(1,h(imin),1,h(kk),1)
  180 if(j.gt.krank) go to 300
      go to 115
c
c        column pivot
c
  200 if(imin.eq.j) go to 230
      call dswap(1,h(j),1,h(imin),1)
      call dswap(m,a(1,j),1,a(1,imin),1)
      call dswap(1,eb(j),1,eb(imin),1)
      call dswap(1,ub(j),1,ub(imin),1)
      call dswap(1,db(j),1,db(imin),1)
      call dswap(1,w(j),1,w(imin),1)
      call iswap(1,ic(j),1,ic(imin),1)
c
c        row pivot
c
  230 continue
      jmax=idamax(mm,a(j,j),1)
      jmax=jmax+j-1
      if(jmax.eq.j) go to 240
      call dswap(n,a(j,1),mda,a(jmax,1),mda)
      call iswap(1,ir(j),1,ir(jmax),1)
  240 continue
c
c     apply householder transformation
c
      tn=dnrm2(mm,a(j,j),1)
      if(tn.eq.0.0d0) go to 170
      if(a(j,j).ne.0.0d0) tn=sign(tn,a(j,j))
      call dscal(mm,1.0d0/tn,a(j,j),1)
      a(j,j)=a(j,j)+1.0d0
      if(j.eq.n) go to 250
      do 248 i=jp1,n
      bb=-ddot(mm,a(j,j),1,a(j,i),1)/a(j,j)
      call daxpy(mm,bb,a(j,j),1,a(j,i),1)
      if(i.le.np) go to 248
      if(h(i).eq.0.0d0) go to 248
      tt=1.0d0-(abs(a(j,i))/h(i))**2
      tt=max(tt,0.0d0)
      t=tt
      tt=1.0d0+.05d0*tt*(h(i)/w(i))**2
      if(tt.eq.1.0d0) go to 244
      h(i)=h(i)*sqrt(t)
      go to 246
  244 continue
      h(i)=dnrm2(m-j,a(j+1,i),1)
      w(i)=h(i)
  246 continue
  248 continue
  250 continue
      h(j)=a(j,j)
      a(j,j)=-tn
c
c
c          update ub, db
c
      ub(j)=ub(j)/abs(a(j,j))
      db(j)=(sign(eb(j),db(j))+db(j))/a(j,j)
      if(j.eq.krank) go to 300
      do 260 i=jp1,krank
      ub(i)=ub(i)+abs(a(j,i))*ub(j)
      db(i)=db(i)-a(j,i)*db(j)
  260 continue
      go to 110
c
c        e n d    m a i n    l o o p
c
  300 continue
c
c        compute ksure
c
      km1=krank-1
      do 318 i=1,km1
      is=0
      kmi=krank-i
      do 315 ii=1,kmi
      if(ub(ii).le.ub(ii+1)) go to 315
      is=1
      temp=ub(ii)
      ub(ii)=ub(ii+1)
      ub(ii+1)=temp
  315 continue
      if(is.eq.0) go to 320
  318 continue
  320 continue
      ksure=0
      sum=0.0d0
      do 328 i=1,krank
      r2=ub(i)*ub(i)
      if(r2+sum.ge.1.0d0) go to 330
      sum=sum+r2
      ksure=ksure+1
  328 continue
  330 continue
c
c     if system is of reduced rank and mode = 2
c     complete the decomposition for shortest least squares solution
c
      if(krank.eq.n .or. mode.lt.2) go to 360
      nmk=n-krank
      kp1=krank+1
      i=krank
  340 tn=dnrm2(nmk,a(i,kp1),mda)/a(i,i)
      tn=a(i,i)*sqrt(1.0d0+tn*tn)
      call dscal(nmk,1.0d0/tn,a(i,kp1),mda)
      w(i)=a(i,i)/tn+1.0d0
      a(i,i)=-tn
      if(i.eq.1) go to 350
      im1=i-1
      do 345 ii=1,im1
      tt=-ddot(nmk,a(ii,kp1),mda,a(i,kp1),mda)/w(i)
      tt=tt-a(ii,i)
      call daxpy(nmk,tt,a(i,kp1),mda,a(ii,kp1),mda)
      a(ii,i)=a(ii,i)+tt*w(i)
  345 continue
      i=i-1
      go to 340
  350 continue
  360 continue
      return
      end
