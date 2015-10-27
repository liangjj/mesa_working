*deck @(#)banmat.f	1.2  10/27/94
      subroutine banmat (in,l1,l2,nt,im,a,ia,y,iy,det,int)
      implicit real *8 (a-h,o-z)
      save
      dimension a(ia,ia), y(iy,iy), int(in)
      n=in
      l=im
      ml=l1
      mu=l2
      ll=ml+mu+1
      if (n.le.1.or.ia.lt.n.or.l1.le.0.or.l2.le.0) go to 190
      lu=ll+1
      lp=ll+ml
      n1=n-1
      if (nt.ne.1) go to 100
      sn=1.d0
      do 30 i=1,ml
      ii=mu+i
      k=ml+1-i
      do 10 j=1,ii
   10 a(i,j)=a(i,j+k)
      k=ii+1
      do 20 j=k,ll
   20 a(i,j)=0.d0
   30 continue
      do 90 nr=1,n1
      np=nr+1
      lr=nr+ml
      if (lr.gt.n) lr=n
      lc=nr+ll-1
      if (lc.gt.n) lc=n
      kk=lc-nr
      mx=nr
      xm=abs(a(nr,1))
      do 40 i=np,lr
      if (abs(a(i,1)).le.xm) go to 40
      mx=i
      xm=abs(a(i,1))
   40 continue
      int(nr)=mx
      if (mx.eq.nr) go to 60
      do 50 i=1,ll
      xx=a(nr,i)
      a(nr,i)=a(mx,i)
   50 a(mx,i)=xx
      sn=-sn
   60 xm=a(nr,1)
      if (xm.eq.0.d0) go to 190
      xm=-1.d0/xm
      do 80 i=np,lr
      j=ll+i-nr
      xx=a(i,1)*xm
      a(nr,j)=xx
      do 70 kt=1,kk
      kt1=(kt-1)*ia
      kti=kt1+i
      ktr=kt1+nr
      a(kti,1)=xx*a(ktr,2)+a(kti,2)
   70 continue
   80 a(i,ll)=0.d0
   90 continue
      if (a(n,1).eq.0.d0) sn=0.d0
  100 det=sn
      if (im.eq.0) return
      do 150 nr=1,n1
      np=nr+1
      lr=nr+ml
      if (lr.gt.n) lr=n
      kk=lr-nr
      if (a(nr,1).eq.0.d0) go to 190
      if (int(nr).eq.nr) go to 120
      j=int(nr)
      do 110 i=1,l
      xx=y(nr,i)
      y(nr,i)=y(j,i)
  110 y(j,i)=xx
  120 do 140 i=1,l
      do 130 kt=1,kk
      ktp=(kt-1)+np
      ktr=(kt-1)*ia+nr
      y(ktp,i)=y(nr,i)*a(ktr,lu)+y(ktp,i)
  130 continue
  140 continue
  150 continue
      xm=a(n,1)
      nr=n
      if (xm.eq.0.d0) go to 190
      xm=1.d0/xm
      do 160 i=1,l
  160 y(n,i)=y(n,i)*xm
      ns=1
      do 180 nb=ns,n1
      nr=n-nb
      np=nr+1
      lc=nr+ll-1
      if (lc.gt.n) lc=n
      kk=lc-nr
      do 170 i=1,l
  170 y(nr,i)=(y(nr,i)-sdot(kk,a(nr,2),ia,y(np,i),1))/a(nr,1)
  180 continue
      return
  190 det=0.d00
      return
      end
