*deck @(#)ptprep.f	5.1  11/6/94
      subroutine ptprep(np,nhi,lmalo,lmahi,lmblo,lmbhi,alpha,rka,rkb
     1 ,argsum)
      implicit real*8(a-h,o-z)
c
c  this routine sets up the points and weights for calculating
c  radial integrals. the quadrature scheme is dependent on argsum.
c      argsum.lt.100  use hermite polynomials on the half
c                     interval (0,infinity).
c      argsum.ge.100  hermites on the interval (-inf,inf).
c
c
      common/ptwtdat/ptpow(50,11),f(50,9,9),pt(50)
      common/ptwt/npts
      common/const/zero,one,two,three,four,five,six,ten
      dimension z(125), w(125)
      dimension temp(50),ba(50,2),bb(50,2),pta(50),ptb(50)
      data f10/1.0d01/,f100/1.0d02/,f1000/1.0d03/,f100000/1.0d05/
      data (z(i),i=1,15)/
     & 0.21686942698851d-01,0.11268419695264d+00,0.27049262038900d+00,
     & 0.48690228923427d+00,0.75304357423052d+00,0.10609308720058d+01,
     & 0.14042548093706d+01,0.17786462185204d+01,0.21817079628455d+01,
     & 0.26130606725136d+01,0.30746179394975d+01,0.35714079774673d+01,
     & 0.41137359184560d+01,0.47235128949706d+01,0.54604887738649d+01/
      data (w(i),i=1,15)/
     & 0.55443354299713d-01,0.12402771569560d+00,0.17529092095372d+00,
     & 0.19148833259240d+00,0.16347380949695d+00,0.10593766037807d+00,
     & 0.50027040043584d-01,0.16442978004616d-01,0.35732067930689d-02,
     & 0.48289694247755d-03,0.37490905185694d-04,0.14936859732817d-05,
     & 0.25527085813265d-07,0.13421789241149d-09,0.95622914518487d-13/
      data (z(i),i=16,40)/
     & 0.10300087285163d-01,0.53985341518803d-01,0.13145451466029d+00,
     & 0.24094787552632d+00,0.38019082766422d+00,0.54664041574787d+00,
     & 0.73769994842358d+00,0.95088776237837d+00,0.11839535321852d+01,
     & 0.14349475698292d+01,0.17022546954938d+01,0.19846055256613d+01,
     & 0.22810767024455d+01,0.25910895759857d+01,0.29144153490462d+01,
     & 0.32511944348558d+01,0.36019794600646d+01,0.39678161994154d+01,
     & 0.43503876059168d+01,0.47522701516031d+01,0.51774085039360d+01,
     & 0.56320649416312d+01,0.61269698576629d+01,0.66832894012089d+01,
     & 0.73569549248072d+01/
      data (w(i),i=16,40)/
     & 0.26397629065285d-01,0.60636947271484d-01,0.92218749059043d-01,
     & 0.11773959347869d+00,0.13265187184088d+00,0.13289895631010d+00,
     & 0.11752399360624d+00,0.90479862760810d-01,0.59653016997591d-01,
     & 0.33091033407772d-01,0.15167308104066d-01,0.56385174832006d-02,
     & 0.16676965314824d-02,0.38448103611385d-03,0.67560321141489d-04,
     & 0.88208709846012d-05,0.83050835920950d-06,0.54370797500722d-07,
     & 0.23635062256294d-08,0.64194892468547d-10,0.10012625782071d-11,
     & 0.79150508564863d-14,0.25895258786105d-16,0.23969524295180d-19,
     & 0.24292790129830d-23/
      data (z(i),i=41,90)/
     & 0.36994166894119d-02,0.19465578559300d-01,0.47723257511078d-01,
     & 0.88299474142332d-01,0.14094424983151d+00,0.20534382473297d+00,
     & 0.28113094473115d+00,0.36789570186710d+00,0.46519663238257d+00,
     & 0.57257157583974d+00,0.68954786429340d+00,0.81565152412903d+00,
     & 0.95041529376704d+00,0.10933853708522d+01,0.12441268939405d+01,
     & 0.14022282325284d+01,0.15673042057959d+01,0.17389983773292d+01,
     & 0.19169845842537d+01,0.21009678588666d+01,0.22906848929554d+01,
     & 0.24859041828723d+01,0.26864259797693d+01,0.28920821561667d+01,
     & 0.31027360886914d+01,0.33182826484263d+01,0.35386483857011d+01,
     & 0.37637919961013d+01,0.39937051598745d+01,0.42284138590013d+01,
     & 0.44679802967817d+01,0.47125055766333d+01,0.49621333441732d+01,
     & 0.52170546662502d+01,0.54775145229933d+01,0.57438204412415d+01,
     & 0.60163540281218d+01,0.62955865198298d+01,0.65821000263865d+01,
     & 0.68766170796566d+01,0.71800426651669d+01,0.74935257051424d+01,
     & 0.78185521503957d+01,0.81570921045017d+01,0.85118452647437d+01,
     & 0.88866800592441d+01,0.92874967414165d+01,0.97241658658846d+01,
     & 0.10215886258578d+02,0.10812986072945d+02/
      data (w(i),i=41,90)/
     & 0.94907437046997d-02,0.22024804243228d-01,0.34374412277947d-01,
     & 0.46295442331579d-01,0.57427067352410d-01,0.67261674334723d-01,
     & 0.75165510542760d-01,0.80448721018930d-01,0.82486665378029d-01,
     & 0.80877284147973d-01,0.75596813121437d-01,0.67098896909613d-01,
     & 0.56303751410186d-01,0.44453364727083d-01,0.32860265552944d-01,
     & 0.22628104702244d-01,0.14442120563853d-01,0.84998891744169d-02,
     & 0.45898158294769d-02,0.22624803215574d-02,0.10129634585731d-02,
     & 0.40985505102578d-03,0.14910339655435d-03,0.48520972503506d-04,
     & 0.14049918795431d-04,0.36005744558413d-05,0.81205488566054d-06,
     & 0.16023622632303d-06,0.27491449821097d-07,0.40739387068860d-08,
     & 0.51773390940542d-09,0.55988313722907d-10,0.51083344028052d-11,
     & 0.38953118086079d-12,0.24563170831570d-13,0.12656133093179d-14,
     & 0.52558265333367d-16,0.17314937367235d-17,0.44419824036643d-19,
     & 0.86801390035972d-21,0.12580502763550d-22,0.13087212865077d-24,
     & 0.93765255356630d-27,0.43860415493907d-29,0.12469057270559d-31,
     & 0.19484701693978d-34,0.14405299157886d-37,0.39431480287447d-41,
     & 0.25121992477205d-45,0.11532495344215d-50/
      data (z(i),i= 91,110)/-.53874808898628d+01,-.46036824578120d+01,
     & -.39447640525915d+01,-.33478545826881d+01,-.27888060754991d+01,
     & -.22549740197087d+01,-.17385377287517d+01,-.12340762292472d+01,
     & -.73747373780326d+00,-.24534071156534d+00, .24534071156534d+00,
     &  .73747373780327d+00, .12340762292472d+01, .17385377287517d+01,
     &  .22549740197087d+01, .27888060754991d+01, .33478545826881d+01,
     &  .39447640525914d+01, .46036824578119d+01, .53874808898629d+01/
      data (w(i),i= 91,110)/ .22293936114860d-12, .43993406245663d-09,
     &  .10860692579369d-06, .78025556477670d-05, .22833861377064d-03,
     &  .32437730854839d-02, .24810519542637d-01, .10901720310356d+00,
     &  .28667550458475d+00, .46224367491990d+00, .46224367491990d+00,
     &  .28667550458473d+00, .10901720310356d+00, .24810519542640d-01,
     &  .32437730854853d-02, .22833861377079d-03, .78025556477728d-05,
     &  .10860692579379d-06, .43993406245699d-09, .22293936114888d-12/
      data (z(i),i=111,120)/-.34361591188396d+01,-.25327316742343d+01,
     & -.17566836493012d+01,-.10366108297905d+01,-.34290132722408d+00,
     &  .34290132722407d+00, .10366108297905d+01, .17566836493012d+01,
     &  .25327316742342d+01, .34361591188395d+01/
      data (w(i),i=111,120)/ .76404328551401d-05, .13436457467716d-02,
     &  .33874394455339d-01, .24013861108199d+00, .61086263373582d+00,
     &  .61086263373583d+00, .24013861108200d+00, .33874394455345d-01,
     &  .13436457467719d-02, .76404328551425d-05/
      data (z(i),i=121,125)/-.20201828704560d+01,-.95857246461381d+00,
     & -.86713183502837d-16, .95857246461379d+00, .20201828704561d+01/
      data (w(i),i=121,125)/ .19953242059054d-01, .39361932315223d+00,
     &  .94530872048294d+00, .39361932315225d+00, .19953242059051d-01/
      save f10,f100,f1000,f100000
      save z,w
c
      if (argsum.ge.one) go to 10
      ilo=1
      ihi=15
      npts=15
      go to 60
   10 if (argsum.ge.f10) go to 20
      ilo=16
      ihi=40
      npts=25
      go to 60
   20 if (argsum.ge.f100) go to 30
      ilo=41
      ihi=90
      npts=50
      go to 60
   30 if (argsum.ge.f1000) go to 40
      ilo=91
      ihi=110
      npts=20
      go to 60
   40 if (argsum.ge.f100000) go to 50
      ilo=111
      ihi=120
      npts=10
      go to 60
   50 ilo=121
      ihi=125
      npts=5
c
c     prepare the integrand.
   60 sqalp=sqrt(alpha)
      sqalpi=one/sqalp
      nhim1=nhi-1
      rkab=rka+rkb
      if(argsum.lt.f100) go to 100
c
c     use scalar bessel routine for short quadrature lengths.
      c=rkab/(two*alpha)
      i1=0
      do 90 i=ilo,ihi
      i1=i1+1
      pt(i1)=c+z(i)*sqalpi
      temp(i1)=w(i)
      do 70 la=lmalo,lmahi
   70 ba(i1,la-lmalo+1)=bess(rka*pt(i1),la-1)
      do 80 lb=lmblo,lmbhi
   80 bb(i1,lb-lmblo+1)=bess(rkb*pt(i1),lb-1)
   90 continue
      go to 140
c
c     use vector bessel routine for long quadrature lengths.
  100 i1=0
      do 110 i=ilo,ihi
      i1=i1+1
      pt(i1)=z(i)*sqalpi
      temp(i1)=exp(rkab*pt(i1))*w(i)
      pta(i1)=pt(i1)*rka
  110 ptb(i1)=pt(i1)*rkb
      do 120 la=lmalo,lmahi
  120 call bessv(ba(1,la-lmalo+1),pta,la-1,npts)
      do 130 lb=lmblo,lmbhi
  130 call bessv(bb(1,lb-lmblo+1),ptb,lb-1,npts)
c
  140 do 150 la=lmalo,lmahi
      do 150 lb=lmblo,lmbhi
      do 150 i=1,npts
  150 f(i,la,lb)=ba(i,la-lmalo+1)*bb(i,lb-lmblo+1)*temp(i)
c
c     prepare an array containing powers of r.
      do 160 i=1,npts
  160 ptpow(i,1)=pt(i)**(np)
      if (nhi.eq.1) go to 180
      do 170 n=1,nhim1
      do 170 i=1,npts
  170 ptpow(i,n+1)=ptpow(i,n)*pt(i)
  180 continue
c
c
c
      return
      end
