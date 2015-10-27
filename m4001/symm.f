*deck @(#)symm.f	5.1 11/6/94
      subroutine symm
      common/csymm/ax(1000),ay(1000),az(1000),q(1000),
     1rmult(1000),lsym,ncell,name(1000)
      data mdel/4h****/
      if(lsym.lt.1)return
      kpt=ncell
      do 500 isym=1,lsym
      read(5,1001)sx1,sx2,sx3,tx
      read(5,1001)sy1,sy2,sy3,ty
      read(5,1001)sz1,sz2,sz3,tz
 1001 format(4f12.8)
      read(5,1002)nord
 1002 format(i1)
      do 400 i=1,ncell
      kpt=kpt+1
      ax(kpt)=sx1*ax(i)+sx2*ay(i)+sx3*az(i)+tx
      ay(kpt)=sy1*ax(i)+sy2*ay(i)+sy3*az(i)+ty
      az(kpt)=sz1*ax(i)+sz2*ay(i)+sz3*az(i)+tz
      name(kpt)=name(i)
      q(kpt)=q(i)
  400 rmult(kpt)=0.0
      if(nord.eq.2)go to 450
      s11=sx1
      s12=sx2
      s13=sx3
      s21=sy1
      s22=sy2
      s23=sy3
      s31=sz1
      s32=sz2
      s33=sz3
      u1=tx
      u2=ty
      u3=tz
      do 2000 ipow=2,nord-1
      a11=s11*sx1+s12*sy1+s13*sz1
      a12=s11*sx2+s12*sy2+s13*sz2
      a13=s11*sx3+s12*sy3+s13*sz3
      b1=s11*tx+s12*ty+s13*tz+u1
      a21=s21*sx1+s22*sy1+s23*sz1
      a22=s21*sx2+s22*sy2+s23*sz2
      a23=s21*sx3+s22*sy3+s23*sz3
      b2=s21*tx+s22*ty+s23*tz+u2
      a31=s31*sx1+s32*sy1+s33*sz1
      a32=s31*sx2+s32*sy2+s33*sz2
      a33=s31*sx3+s32*sy3+s33*sz3
      b3=s31*tx+s32*ty+s33*tz+u3
      sx1=a11
      sx2=a12
      sx3=a13
      sy1=a21
      sy2=a22
      sy3=a23
      sz1=a31
      sz2=a32
      sz3=a33
      tx=b1
      ty=b2
      tz=b3
      do 2400 i=1,ncell
      kpt=kpt+1
      ax(kpt)=sx1*ax(i)+sx2*ay(i)+sx3*az(i)+tx
      ay(kpt)=sy1*ax(i)+sy2*ay(i)+sy3*az(i)+ty
      az(kpt)=sz1*ax(i)+sz2*ay(i)+sz3*az(i)+tz
      name(kpt)=name(i)
      q(kpt)=q(i)
 2400 rmult(kpt)=0.0
 2000 continue
  450 ncell=ncell*nord
  500 continue
      do 600 i=1,kpt
  505 if(ax(i).ge.0.0)go to 510
      ax(i)=ax(i)+1.0
      go to 505
  510 if(ax(i).le.1.00001)go to 520
      ax(i)=ax(i)-1.0
      go to 510
  520 if(ay(i).ge.0.0)go to 530
      ay(i)=ay(i)+1.0
      go to 520
  530 if(ay(i).le.1.00001)go to 540
      ay(i)=ay(i)-1.0
      go to 530
  540 if(az(i).ge.0.0)go to 550
      az(i)=az(i)+1.0
      go to 540
  550 if(az(i).le.1.00001)go to 600
      az(i)=az(i)-1.0
      go to 550
  600 continue
      max=kpt
      do 800 i=1,kpt
      if(name(i).eq.mdel)go to 800
      ip=i+1
      do 799 j=ip,kpt
      if(abs(ax(i)-ax(j)).gt.0.0001.and.abs(abs(ax(i)-ax(j))-1.0)
     1.gt.0.0001)go to 799
      if(abs(ay(i)-ay(j)).gt.0.0001.and.abs(abs(ay(i)-ay(j))-1.0)
     1.gt.0.0001)go to 799
      if(abs(az(i)-az(j)).gt.0.0001.and.abs(abs(az(i)-az(j))-1.0)
     1.gt.0.0001)go to 799
      name(j)=mdel
      max=max-1
  799 continue
  800 continue
      inc=0
      do 900 i=1,max
  810 j=i+inc
      if(name(j).ne.mdel)go to 820
      inc=inc+1
      go to 810
  820 name(i)=name(j)
      ax(i)=ax(j)
      ay(i)=ay(j)
      az(i)=az(j)
      rmult(i)=rmult(j)
      q(i)=q(j)
  900 continue
      ncell=max
      return
      end
