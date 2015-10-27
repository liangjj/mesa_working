*deck glss
      subroutine glss(im,in,ik,il,nr,a,ia,y,iy,b,ib,x,jx,eps)
      implicit real *8 (a-h,o-z)
      dimension a(ia,ia),b(ib,ib),y(iy,iy),x(jx,jx)
      equivalence (dot,t)
      m=im
      n=in
      k=ik
      l=il
      if(m.le.0.or.n.le.0.or.k.le.0.or.l.le.0)goto98
      if(m.gt.ia.or.m.gt.iy.or.l.gt.ib.or.l.gt.jx)goto99
      goto100
   98 write(6,1515)
      call exit
   99 write(6,2525)
      call exit
  100 do1i=1,l
      do1j=1,k
    1 x(i,j)=0.d00
      sc=0.d0
      do50i=1,n
      dot=sdot(m,a(1,i),1,a(1,i),1)
      if(t.le.sc) go to 50
      sc = t
      ix = i
   50 continue
      nr=0
      if(sc.eq.0.d0)goto27
      yy=sqrt(sc)
      yx=eps*yy
      go to 29
   30 ic = nr + 1
      t = 0.0d0
      ix = 0
      do 26 i = ic,n
      yy =sdot(m,a(1,i),1,a(1,i),1)
      if(yy.le.t) go to 26
      t = yy
      ix=i
   26 continue
      if(ix.eq.0) go to 24
      yy=sqrt(t)
      if(yy.lt.yx) go to 24
   29 nr=nr+1
      if(ix.eq.nr) go to 31
      ic = ix
      do32i=1,m
      t=a(i,ic)
      a(i,ic)=a(i,nr)
   32 a(i,nr)=t
      do33i=1,l
      t=b(i,ic)
      b(i,ic)=b(i,nr)
   33 b(i,nr)=t
   31 yy=1.0d0/yy
      call sscal(m,yy,a(1,nr),1)
      call sscal(l,yy,b(1,nr),1)
      if(nr.eq.n) go to 34
      ip1=nr+1
      do10j=ip1,n
      dot=sdot(m,a(1,nr),1,a(1,j),1)
      t = -t
      call saxpy(m,t,a(1,nr),1,a(1,j),1)
      call saxpy(l,t,b(1,nr),1,b(1,j),1)
   10 continue
   34 do20j=1,k
      dot=sdot(m,a(1,nr),1,y(1,j),1)
      t = -t
      call saxpy(m,t,a(1,nr),1,y(1,j),1)
      call saxpy(l,t,b(1,nr),1,x(1,j),1)
   20 continue
      if(nr.lt.n) go to 30
   24 continue
      do25j=1,k
      do25i=1,l
   25 x(i,j)=-x(i,j)
   27 return
 1515 format(/,5x,'Glss Called With Arrays Dimensioned Wrong')
 2525 format(/,5x,'Glss Called With Non-Positive Indices')
      end
