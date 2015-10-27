*deck genztb
* generates the ztab array for ecp integrals for arbritary angular
* momentum. courtesy tony rappe, csu, 3/3/93.
      program genztb
      implicit real*8 (a-h,o-z)
c     Mesa commons:
      common/const/zero,one,two,three,four,five,six,ten
      common/ztabcm/lf(8),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      common/ztabwp/zlm(130)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c     realsh data
      integer*4 lx(1000),ly(1000),lz(1000)
      real*8 zz(1000)
      integer*4 lmmf(150),lmml(150),lff(20)
      call ldata
c     call MESA version:
      call ztab
      call realsh(lmax,lx,ly,lz,zz,lmmf,lmml,lff)
      call lcomp(lx,ly,lz,zz,zlm,lmmf,lmml)
      call wrtdat(lmax,lx,ly,lz,zz,lmmf,lmml,lff)
      stop
      end
      block data ltab
      implicit real*8(a-h,o-z)
c
c     this routine sets up the real spherical harmonics in the form of
c     linear combinations of cartesian products-(l,m)#l(l+1)-m+1
c     coded through i functions (l=6).
c
c     11 july 1990     rlm at lanl
c        breaking up the /ztabcm/ common block so that real*8 and integers
c        are in different blocks.  real*8 on some machines must start
c        on an even word boundary.
c
      common/ztabcm/lf(8),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      data lf /1,2,5,10,17,26,37,130/
      data lmf /1,2,3,4,5,7,8,10,11,12,14,16,18,20,22,23,25,28,30,34,36,
     1 39,41,43,45,47,50,53,57,61,64,67,70,72,76,78,81,84,87,93,97,103,1
     2 06,110,113,116,120,124,127/
      data lml /1,2,3,4,6,7,9,10,11,13,15,17,19,21,22,24,27,29,33,35,38,
     1 40,42,44,46,49,52,56,60,63,66,69,71,75,77,80,83,86,92,96,102,105,
     2 109,112,115,119,123,126,130/
      data lmx /0,1,0,0,2,0,1,0,0,0,1,3,1,2,0,1,1,0,0,0,0,1,2,0,4,2,0,3,
     1 1,2,0,0,2,1,1,0,0,0,0,0,1,1,2,0,3,1,5,3,1,0,2,4,3,1,1,3,2,0,2,0,1
     2 ,1,1,0,0,0,0,0,0,1,1,2,2,0,0,3,1,4,2,0,5,3,1,5,3,1,4,2,0,4,2,0,3,
     3 1,1,3,2,0,2,0,2,0,1,1,1,0,0,0,0,0,0,0,1,1,1,2,0,2,0,3,1,3,1,4,2,0
     4 ,6,4,2,0/
      data lmy /0,0,0,1,0,2,0,0,0,1,1,0,2,0,2,0,0,0,0,1,1,1,1,3,0,2,4,0,
     1 2,0,2,2,0,0,0,0,0,0,1,1,1,1,1,3,1,3,0,2,4,4,2,0,0,2,2,0,0,2,0,2,0
     2 ,0,0,0,0,0,1,1,1,1,1,1,1,3,3,1,3,1,3,5,1,3,5,0,2,4,0,2,4,0,2,4,0,
     3 2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,5
     4 ,0,2,4,6/
      data lmz /0,0,1,0,0,0,1,2,0,1,0,0,0,1,1,2,0,3,1,2,0,1,0,0,0,0,0,1,
     1 1,2,2,0,0,3,1,4,2,0,3,1,2,0,1,1,0,0,0,0,0,1,1,1,2,2,0,0,3,3,1,1,4
     2 ,2,0,5,3,1,4,2,0,3,1,2,0,2,0,1,1,0,0,0,0,0,0,1,1,1,2,2,2,0,0,0,3,
     3 3,1,1,4,4,2,2,0,0,5,3,1,6,4,2,0,5,3,1,4,2,0,3,3,1,1,2,2,0,0,1,1,1
     4 ,0,0,0,0/
c
      end
      subroutine lcomp(lx,ly,lz,zz,zlm,lmmf,lmml)
      integer*4 lx(1),ly(1),lz(1),lmmf(1),lmml(1)
      real*8 zz(1),zlm(1)
      common/ztabcm/lf(8),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      integer*4 lp(8)
      lp(1)=1
      lp(2)=2
      lp(3)=5
      lp(4)=12
      lp(5)=25
      lp(6)=47
      lp(7)=81
      lp(8)=131
c     lf indicates the first term for a given l value
c     lmf indicates the first term for a m value given l
c     lml indicates the last term for a m value given l
      ip=0
      do 25 l=0,6
      lstrt=lp(l+1)
      lstop=lp(l+2)-1
c     do 20 i=1,l+l+1
c     m=l-i+1
      do 15 j=lstrt,lstop
      do 10 k=lstrt,lstop
c     write(6,600) lmf(i+lf(l+1)-1),lmmf(i+lf(l+1)-1),
c    $ lml(i+lf(l+1)-1),lmml(i+lf(l+1)-1)
c 600 format(4i5)
c     do 15 j=lmf(i+lf(l+1)-1),lml(i+lf(l+1)-1)
c     do 10 k=lmf(i+lf(l+1)-1),lml(i+lf(l+1)-1)
      if(lmx(j).ne.lx(k)) go to 5
      if(lmy(j).ne.ly(k)) go to 5
      if(lmz(j).ne.lz(k)) go to 5
      ip=ip+1
      write(6,100) ip,lmx(j),lmy(j),lmz(j),lx(k),ly(k),lz(k),zlm(j),
     $ zz(k),zz(k)-zlm(j)
  100 format(i4,6i2,1p3d12.5)
    5 continue
   10 continue
   15 continue
   20 continue
   25 continue
      return
      end
      subroutine ztab
      implicit real*8(a-h,o-z)
c
c     this routine sets up the real spherical harmonics in the form of
c     linear combinations of cartesian products-(l,m)#l(l+1)-m+1
c     coded through i functions (l=6).
c
c     11 july 1990     rlm at lanl
c        breaking up the /ztabcm/ common block so that real*8 and integers
c        are in different blocks.  real*8 on some machines must start
c        on an even word boundary.
c
      common/const/zero,one,two,three,four,five,six,ten
      common/ztabwp/zlm(130)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      zlm(1)=sqrt(one/fpi)
      zlm(2)=sqrt(three/fpi)
      zlm(3)=zlm(2)
      zlm(4)=zlm(2)
      zlm(5)=sqrt(15.0d0/fpi)/two
      zlm(6)=-zlm(5)
      zlm(7)=two*zlm(5)
      zlm(8)=three*sqrt(five/fpi)/two
      zlm(9)=-zlm(8)/three
      zlm(10)=zlm(7)
      zlm(11)=zlm(7)
      zlm(12)=sqrt(35.0d0/(8.0d0*fpi))
      zlm(13)=-three*zlm(12)
      zlm(14)=sqrt(105.0d0/(four*fpi))
      zlm(15)=-zlm(14)
      zlm(16)=five*sqrt(21.0d0/(8.0d0*fpi))
      zlm(17)=-zlm(16)/five
      zlm(18)=five*sqrt(7.0d0/fpi)/two
      zlm(19)=-three*zlm(18)/five
      zlm(20)=zlm(16)
      zlm(21)=zlm(17)
      zlm(22)=two*zlm(14)
      zlm(23)=-zlm(13)
      zlm(24)=-zlm(12)
      zlm(25)=sqrt(315.0d0/(64.0d0*fpi))
      zlm(26)=-six*zlm(25)
      zlm(27)=zlm(25)
      zlm(28)=sqrt(315.0d0/(8.0d0*fpi))
      zlm(29)=-three*zlm(28)
      temp=sqrt(45.0d0/fpi)/four
      zlm(30)=7.0d0*temp
      zlm(31)=-zlm(30)
      zlm(32)=temp
      zlm(33)=-temp
      temp=sqrt(45.0d0/(8.0d0*fpi))
      zlm(34)=7.0d0*temp
      zlm(35)=-three*temp
      temp=sqrt(9.0d0/fpi)/8.0d0
      zlm(36)=35.0d0*temp
      zlm(37)=-30.0d0*temp
      zlm(38)=three*temp
      zlm(39)=zlm(34)
      zlm(40)=zlm(35)
      temp=sqrt(45.0d0/(four*fpi))
      zlm(41)=7.0d0*temp
      zlm(42)=-temp
      zlm(43)=-zlm(29)
      zlm(44)=-zlm(28)
      zlm(45)=sqrt(315.0d0/(four*fpi))
      zlm(46)=-zlm(45)
      zlm(47)=sqrt(693.0d0/(128.0d0*fpi))
      zlm(48)=-10.0d0*zlm(47)
      zlm(49)=five*zlm(47)
      zlm(50)=sqrt(3465.0d0/(64.0d0*fpi))
      zlm(51)=-six*zlm(50)
      zlm(52)=zlm(50)
      temp=sqrt(385.0d0/(128.0d0*fpi))
      zlm(53)=9.0d0*temp
      zlm(54)=-27.0d0*temp
      zlm(55)=three*temp
      zlm(56)=-temp
      temp=sqrt(1155.0d0/fpi)/four
      zlm(57)=three*temp
      zlm(58)=-zlm(57)
      zlm(59)=-temp
      zlm(60)=+temp
      temp=sqrt(165.0d0/fpi)/8.0d0
      zlm(61)=21.0d0*temp
      zlm(62)=-14.0d0*temp
      zlm(63)=temp
      temp=sqrt(11.0d0/fpi)/8.0d0
      zlm(64)=63.0d0*temp
      zlm(65)=-70.0d0*temp
      zlm(66)=15.0d0*temp
      zlm(67)=zlm(61)
      zlm(68)=zlm(62)
      zlm(69)=zlm(63)
      temp=sqrt(1155.0d0/fpi)/two
      zlm(70)=three*temp
      zlm(71)=-temp
      zlm(72)=-zlm(54)
      zlm(73)=-zlm(55)
      zlm(74)=-zlm(53)
      zlm(75)=-zlm(56)
      zlm(76)=sqrt(3465.0d0/fpi)/two
      zlm(77)=-zlm(76)
      zlm(78)=zlm(49)
      zlm(79)=zlm(48)
      zlm(80)=zlm(47)
      temp=sqrt(3003.0d0/(512.0d0*fpi))
      zlm(81)=6.0d0*temp
      zlm(82)=-20.0d0*temp
      zlm(83)=zlm(81)
      zlm(84)=sqrt(9009.0d0/(128.0d0*fpi))
      zlm(85)=-10.0d0*zlm(84)
      zlm(86)=5.0d0*zlm(84)
      temp=sqrt(819.0d0/(256.0d0*fpi))
      zlm(87)=11.0d0*temp
      zlm(88)=-66.0d0*temp
      zlm(89)=zlm(87)
      zlm(90)=-temp
      zlm(91)=6.0d0*temp
      zlm(92)=-temp
      temp=sqrt(1365.0d0/(128.0d0*fpi))
      zlm(93)=11.0d0*temp
      zlm(94)=-33.0d0*temp
      zlm(95)=9.0d0*temp
      zlm(96)=-3.0d0*temp
      temp=sqrt(1365.0d0/(512.0d0*fpi))
      zlm(97)=33.0d0*temp
      zlm(98)=-zlm(97)
      zlm(99)=-18.0d0*temp
      zlm(100)=+18.0d0*temp
      zlm(101)=temp
      zlm(102)=-temp
      temp=sqrt(273.0d0/fpi)/8.0d0
      zlm(103)=33.0d0*temp
      zlm(104)=-30.0d0*temp
      zlm(105)=5.0d0*temp
      temp=sqrt(13.0d0/fpi)/16.0d0
      zlm(106)=231.0d0*temp
      zlm(107)=-315.0d0*temp
      zlm(108)=105.0d0*temp
      zlm(109)=-5.0d0*temp
      zlm(110)=zlm(103)
      zlm(111)=zlm(104)
      zlm(112)=zlm(105)
      temp=sqrt(1365.0d0/(128.0d0*fpi))
      zlm(113)=33.0d0*temp
      zlm(114)=-18.0d0*temp
      zlm(115)=temp
      zlm(116)=-zlm(94)
      zlm(117)=-zlm(93)
      zlm(118)=-zlm(95)
      zlm(119)=-zlm(96)
      temp=sqrt(819.0d0/fpi)/4.0d0
      zlm(120)=11.0d0*temp
      zlm(121)=-zlm(120)
      zlm(122)=-temp
      zlm(123)=temp
      zlm(124)=zlm(86)
      zlm(125)=zlm(85)
      zlm(126)=zlm(84)
      zlm(127)=sqrt(3003.0d0/(512.0d0*fpi))
      zlm(128)=-15.0d0*zlm(127)
      zlm(129)=-zlm(128)
      zlm(130)=-zlm(127)
c
c
      return
      end
      subroutine ldata
      implicit real*8(a-h,o-z)
c
c     load useful constants.
c
      common/dfac/dfac(23)
      common/fact/fac(13),fprod(7,7)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/const/zero,one,two,three,four,five,six,ten
c     data zero/0.0d0/,one/1.0d0/,two/2.0d0/,three/3.0d0/,four/4.0d0/
c     data five/5.0d0/,six/6.0d0/,ten/10.0d0/
c     save zero,one,two,three,four,five,six,ten
c
      zero=0.0d+00
      one=1.0d+00
      two=2.0d+00
      three=3.0d+00
      four=4.0d+00
      five=5.0d+00
      six=6.0d+00
      ten=10.0d+00
c
      pi=four*atan(one)
      twopi=two*pi
      fpi=four*pi
      sqpi=sqrt(pi)
      sqpi2=sqrt(pi/two)
      pi3haf=pi*sqrt(pi)
      pi5hf2=twopi*pi3haf
      piquart=twopi/sqrt(pi5hf2)
      dfac(1)=one
      dfac(2)=one
      dfac(3)=one
      do 10 i=4,23
   10    dfac(i)=dfac(i-2)*float(i-2)
      fac(1)=one
      do 20 i=1,12
   20    fac(i+1)=fac(i)*float(i)
      do 30 l1=1,7
         do 30 k1=1,l1
            k=k1-1
            fprod(k1,l1)=fac(l1+k)/(fac(k1)*fac(l1-k1+1))
   30 continue
c
c
      return
      end
      subroutine realsh(lmax,lx,ly,lz,zz,lmmf,lmml,lff)
      implicit none
      integer*4 lx(1000),ly(1000),lz(1000),lmmf(49),lmml(49),lff(1)
      real*8 zz(1000)
      integer*4 i,j,k,l,m,kstart,kst,kk,kz,ip,ipp,mm,imp
      complex*16 bi,b,tmpc1,tmpc2
      real*8 fact,halfg,z,xc,xs,yc,ys,fac0,hlfg0,z0,xc0,xs0,yc0,ys0
      real*8 c0,c,s0,s
      dimension fact(40),halfg(40)
      dimension z(40)
      dimension xc(40),xs(40),yc(40),ys(40)
      common/factt/fac0,fact
      common/hlfgg/hlfg0,halfg
      common/zz/z0,z
      common/xx/xc0,xc,xs0,xs
      common/yy/yc0,yc,ys0,ys
      common/aa/c0,c,s0,s
      real*8 pi,sqrt2,sqrtpi
      real*8 stpg,step,hlfg,fac,one,anrm,coef,zzero,coefc,coefs
      integer*4 lmax,ixp,iyp
c     dimensioned for l=15,2l+1
      integer*4 il(1000),im(1000),ix(1000),iy(1000),iz(1000)
      real*8 zt(1000)
      b=(0.0d0,-1.0d0)
      bi=(0.0d0,1.0d0)
      pi=4.0d0*datan(1.0d0)
      sqrt2=1.0d0/dsqrt(2.0d0)
      sqrtpi=1.0d0/dsqrt(pi)
c     generate factorial array and half gamma array
      stpg=0.5d0
      step=1.0d0
      hlfg=0.5d0
      fac=1.0d0
      fac0=1.0d0
      hlfg0=0.5d0
      do 10 i=1,40
      fact(i)=fac
      stpg=stpg+1.0d0
      hlfg=hlfg*stpg
      halfg(i)=hlfg
      step=step+1.0d0
      fac=fac*step
   10 continue
      write(6,1)
    1 format(' enter lmax')
      read(5,2) lmax
    2 format(i2)
      ip=0
      do 110 l=0,lmax
      one=1.0d0
      if(mod(l,2).ne.0) one=-1.0d0
      do 15 k=0,l
   15 z(k)=0.0d0
c     generate m=0 terms
      anrm=sqrtpi*dsqrt(halfg(l)/(2.0d0*fact(l)*fact(l+l)))
      anrm=anrm*one
      m=0
      kstart=l/2
      do 30 k=kstart,l
      step=dfloat(k+k-l)
      coef=1.0d0
      do 20 j=1,l
      step=step+1.0d0
   20 coef=coef*step
      coef=((-1.0d0)**k)*coef*fact(l)/(fact(l-k)*fact(k))
      coef=coef*anrm
   30 z(k+k-l)=coef
      if(z(0).ne.0.0d0) then
      zzero=dabs(z(0))
      else
      zzero=dabs(z(1))
      endif
      do 40 k=0,l
   40 z(k)=z(k)/zzero
      do 45 kk=0,l
      k=l-kk
      if(z(k).ne.0.0d0) then
      ip=ip+1
      il(ip)=l
      im(ip)=0
      ix(ip)=0
      iy(ip)=0
      iz(ip)=k
      zt(ip)=z(k)*zzero
      endif
   45 continue
c     generate -m and m terms for given l
      do 100 m=1,l
      do 50 k=0,l
      xc(k)=0.0d0
      yc(k)=0.0d0
      xs(k)=0.0d0
      ys(k)=0.0d0
   50 z(k)=0.0d0
c     generate normalization term
      anrm=dsqrt(halfg(l)*fact(l-m)/(2.0d0*fact(l)*fact(l+l)*fact(l+m)))
      anrm=anrm*sqrtpi*one
c     generate z of m terms
      kstart=(l+m)/2
      do 70 k=kstart,l
      step=dfloat(k+k-l-m)
      coef=1.0d0
      kst=l+m
      do 60 j=1,kst
      step=step+1.0d0
   60 coef=coef*step
      coef=((-1.0d0)**k)*coef*fact(l)/(fact(l-k)*fact(k))
      coef=coef*anrm
   70 z(k+k-l-m)=coef
      if(z(0).ne.0.0d0) then
      zzero=dabs(z(0))
      else
      zzero=dabs(z(1))
      endif
      do 75 k=0,l
      z(k)=z(k)/zzero
   75 continue
c     generate (x-iy)**m + (x+iy)**m terms
c     referred to as +m terms
      do 76 kk=0,l
      kz=l-kk
      do 80 k=0,m
      coef=sqrt2*(fact(m)/(fact(m-k)*fact(k)))
      tmpc1=b**k
      tmpc2=bi**k
      coefc=coef*(tmpc1+tmpc2)*zzero
      ixp=m-k
      iyp=k
      if(coefc.ne.0.0d0) then
      if(z(kz).ne.0.0d0) then
      ip=ip+1
      il(ip)=l
      im(ip)=m
      ix(ip)=ixp
      iy(ip)=iyp
      iz(ip)=kz
      zt(ip)=coefc*z(kz)
c     write(6,7) m,ixp,iyp,kz,coefc*z(kz)
c   7 format(' m = ',i3,' x**',i2,' y**',i2,' z**',i2,
c    $ ' coef = ',1pd22.15)
      endif
      endif
   80 continue
   76 continue
c     generate i*((x-iy)**m - (x+iy)**m) terms
c     referred to as -m terms
      do 86 kk=0,l
      kz=l-kk
      do 90 k=0,m
      coef=sqrt2*(fact(m)/(fact(m-k)*fact(k)))
      tmpc1=b**k
      tmpc2=bi**k
      coefs=coef*(tmpc1-tmpc2)*bi*zzero
      ixp=m-k
      iyp=k
      if(coefs.ne.0.0d0) then
      if(z(kz).ne.0.0d0) then
      ip=ip+1
      il(ip)=l
      im(ip)=-m
      ix(ip)=ixp
      iy(ip)=iyp
      iz(ip)=kz
      zt(ip)=coefs*z(kz)
c     write(6,7) -m,ixp,iyp,kz,coefs*z(kz)
      endif
      endif
   90 continue
   86 continue
  100 continue
  110 continue
      do 120 i=1,ip
      write(6,105) i,il(i),im(i),ix(i),iy(i),iz(i),zt(i)
  105 format(i5,' l = ',i2, ' m = ',i2,3i3,1pd22.15)
  120 continue
c     sort to decreasing m order
      write(6,106) 1,ix(1),iy(1),iz(1),zt(1)
      ipp=1
      lx(ipp)=ix(1)
      ly(ipp)=iy(1)
      lz(ipp)=iz(1)
      zz(ipp)=zt(1)
      lmmf(1)=1
      lmml(1)=1
      imp=1
      lff(1)=1
      lff(2)=2
      do 300 l=1,lmax
      do 299 mm=-l,l
      m=-mm
      lmml(imp)=ipp
      imp=imp+1
      lmmf(imp)=ipp+1
c     write(6,691) l,m,lml(imp-1),lmf(imp)
  691 format(5i5)
c     loop to take care of multiple terms per m value
      do 295 k=1,l+3
      do 290 i=1,ip
      if(il(i).eq.l) then
      if(im(i).eq.m) then
      il(i)=0
      ipp=ipp+1
      lx(ipp)=ix(i)
      ly(ipp)=iy(i)
      lz(ipp)=iz(i)
      zz(ipp)=zt(i)
      write(6,106) ipp,ix(i),iy(i),iz(i),zt(i)
  106 format(i3,3i2,1pd22.15)
      goto 290
      endif
      endif
  290 continue
  295 continue
  299 continue
      lff(l+2)=imp+1
  300 continue
      lmml(imp)=ipp
      write(6,610) (lff(l),l=1,lmax+1)
  610 format(16i5)
      return
      end
      subroutine wrtdat(lmax,lx,ly,lz,zz,lmmf,lmml,lff)
      integer*4 lmax,lx(1),ly(1),lz(1),lmmf(1),lmml(1)
      integer*4 lff(1)
      real*8 zz(1)
      write(6,10) (lff(i),i=1,lmax+1)
   10 format('      data lf/',16i4)
      lmm=0
      do 20 l=0,lmax
      lmm=lmm+l+l+1
   20 continue
      write(6,30) (lmmf(i),i=1,lmm)
   30 format('      data lmf/',12i4,/,('     $',16i4))
      write(6,40) (lmml(i),i=1,lmm)
   40 format('      data lmf/',12i4,/,('     $',16i4))
      lxm=lmml(lmm)
      write(6,50) (lx(i),i=1,lxm)
   50 format('      data lmx/',16i3,/,('     $',18i3))
      write(6,60) (lx(i),i=1,lxm)
   60 format('      data lmy/',16i3,/,('     $',18i3))
      write(6,70) (lx(i),i=1,lxm)
   70 format('      data lmz/',16i3,/,('     $',18i3))
      write(6,80) (i,zz(i),i=1,lxm)
   80 format('      zlm(',i3,')=',1pd22.15)
      return
      end
