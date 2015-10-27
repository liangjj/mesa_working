*deck @(#)lebedev.f	5.1  11/28/95
      subroutine lebedev(mxleb,p,w,l,npts)
c***begin prologue     lebedev.f
c***date written       930427  
c***revision date      4/17/95      
c
c   february 15,1995   rlm at lanl
c      adding quadrature for l=35 from treutler and alrichs, jcp, 346(1995).
c   january 24, 1995   rlm at lanl
c      adding quadrature for l=41 from lebedev and skorokhodov,
c         russia acad. sci. dokl. math., 45, 587(1992).
c         this reference has l=41,47, and 53.
c   january 23, 1995   rlm at lanl
c      adding quadratures for l=29 from lebedev, sibirskii matematicheskii
c         zhurnal, 18, 132(1975).
c         this reference has l=(25,27,29).
c   january 20, 1995   rlm at lanl
c      adding quadratures for l=5,7 from abramowitz and stegun, p. 894.
c   march 13, 1994     rlm at lanl
c      fixing error in number of lebedev of points for l=13
c***keywords           spherical harmonics, numerical integration,
c                      angular quadrature 
c***author             r.l. martin, lanl
c***source             @(#)lebedev.f	5.1   11/28/95
c***purpose            returns the points and weights associated with
c                      gaussian quadrature on the unit sphere.
c***description
c                      this routine returns the points and weights associated
c                      with integration on the unit sphere. in particular,
c                      quadratures of a given order l exactly integrate
c                      all surface harmonics of degree l or less. the base
c                      points of the quadrature are invariant under the
c                      octahedral group plus inversion. the jacobian
c                      [sin(theta)dtheta dphi] is implicitly included, so
c                      the working formula for integrating a polynomial f
c                      using this routine is simply Sum f(p(i))*w(i).
c
c                  mxleb  ---  leading dimension of p.
c                      p  ---  quadrature points (x/r,y/r,z/r), output.
c                      w  ---  quadrature weights, output.
c                      l  ---  order requested, input.
c                   npts  ---  number of points in the quadrature, input.
c                              the correspondence between order and npts is:
c                              l     npts
c                              ----------
c                              3        6
c                              5       18
c                              7       26
c                              9       38
c                             11       50
c                             13       74
c                             15       86
c                             17      110
c                             19      146
c                             23      194
c                             29      302
c                             35      434
c                             41      590
c                             47      770
c                             53      974
c
c    
c
c***references
c                      v.i. lebedev, zh. vychisl. mat. mat. fiz, 15,48(1975)
c                                    ibid.,16,293(1976).
c                      see also
c                      abramowitz and stegun, p. 894.
c                      a.d. mclaren, math. comput. 17, 361(1963).
c                      a.h. stroud, approximate calculation of multiple
c                                   integrals (prentice-hall,1971)
c***routines called
c                      none
c***end prologue       lebedev.f
      implicit none
c
c     --- input variables ---
      integer mxleb,l,npts
c     --- input arrays (unmodified) ---
c     --- input arrays (modified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 p(mxleb,3),w(mxleb)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer pts,i,a,b,c,k,n1,n2,n3
      integer n1max,n2max,n3max
      parameter (n1max=12,n2max=4,n3max=12)
      real*8 a1,a2,a3,signa,signb,signc,root2i,root3i
      real*8 zero,half,one,two,three,four,pi
      data zero/0.0d+00/, half/0.5d+00/, one/1.0d+00/ 
      data two/2.0d+00/, three/3.0d+00/, four/4.0d+00/ 
      real*8 bk(n1max),lb(n1max),mb(n1max)
      real*8 ck(n2max),pc(n2max),qc(n2max)
      real*8 dk(n3max),rd(n3max),ud(n3max),wd(n3max)
c
      save zero,half,one,two,three,four,pi
c
      integer inp,iout
      common/io/inp,iout
c
c     --- set the quadrature weights depending on order ---
      if (l.eq.3) then
         a1=one/6.0d0
         a2=zero
         a3=zero
         n1=0
         n2=0
         n3=0
      else if(l.eq.5) then
c        this from abramowitz and stegun, p. 894.
c        note that the icosahedral formula (sphere.f) is more compact.
         a1=one/30.0d+00
         a2=one/15.0d+00
         a3=zero
         n1=0
         n2=0
         n3=0
      else if(l.eq.7) then
c        this from abramowitz and stegun, p. 894.
         a1=40.0d+00/840.0d+00
         a2=32.0d+00/840.0d+00
         a3=27.0d+00/840.0d+00
         n1=0
         n2=0
         n3=0
      else if(l.eq.9) then
         a1=one/105.0d+00
         a2=zero
         a3=9.0d+00/280.0d+00
         n1=0
         n2=1
         ck(1)=one/35.0d+00
         pc(1)=sqrt((one+sqrt(one-four*(1.0d0/6.0d0)))/two)
         qc(1)=sqrt(one-pc(1)*pc(1))
         n3=0
      else if(l.eq.11) then
         a1=4.0d+00/315.0d+00
         a2=64.0d+00/2835.0d+00
         a3=27.0d+00/1280.0d+00
         n1=1
         bk(1)=121.0d+00*121.0d+00/725760.0d+00 
         mb(1)=three/sqrt(11.0d+00)
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         n2=0
         n3=0
      else if(l.eq.13) then
         a1=16.0d+00/31185.0d+00
         a2=128.0d+00*43.0d+00/331485.0d+00
c        note that quadrature has negative weight.
         a3=-729.0d+00/24640.0d+00
         n1=1
         bk(1)=(13.0d+00**5)/13970880.0d+00 
         mb(1)=sqrt(7.0d+00/13.0d+00)
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         n2=1
         ck(1)=65.0d+00*65.0d+00/255717.0d+00
         pc(1)=sqrt((one+sqrt(41.0d+00/65.0d+00))/2.0d+00)
         qc(1)=sqrt(one-pc(1)*pc(1))
         n3=0
      else if(l.eq.15) then
         a1=0.0115440115441d+00
         a2=zero
         a3=0.0119439090859d+00
         n1=2
         bk(1)=0.0111105557106d+00
         mb(1)=0.852518311701d+00
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         bk(2)=0.0118765012945d+00
         mb(2)=0.189063552885d+00
         lb(2)=sqrt(half*(one-mb(2)*mb(2)))
         n2=1
         ck(1)=0.0118123037469d+00
         pc(1)=0.927330657151d+00
         qc(1)=sqrt(one-pc(1)*pc(1))
         n3=0
      else if(l.eq.17) then
c        these differ from those in original paper.
c        see notes in zh. vychisl. mat.mat.fiz. 16, 293(1976).
         a1=0.00382827049494d+00
         a2=zero
         a3=0.009793737651249d+00
         n1=3
         bk(1)=0.00821173728319d+00
         mb(1)=0.965124035087d+00
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         bk(2)=0.00959547133607d+00
         mb(2)=0.828769981253d+00
         lb(2)=sqrt(half*(one-mb(2)*mb(2)))
         bk(3)=0.00994281489118d+00
         mb(3)=0.215957291846d+00
         lb(3)=sqrt(half*(one-mb(3)*mb(3)))
         n2=1
         ck(1)=0.00969499636166d+00
         pc(1)=0.878158910604d+00
         qc(1)=sqrt(one-pc(1)*pc(1))
         n3=0
      else if(l.eq.19) then
         a1=1856.0d+00/3095235.0d+00
         a2=606208.0d+00/82219995.0d+00
         a3=6490935.0d+00/900204032.0d+00
         n1=3
         bk(1)=7.57439415905d-03
         mb(1)=0.974888643677d+00
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         bk(2)=6.75382948631d-03
         mb(2)=0.807089818360d+00
         lb(2)=sqrt(half*(one-mb(2)*mb(2)))
         bk(3)=7.11635549312d-03
         mb(3)=0.291298882210d+00
         lb(3)=sqrt(half*(one-mb(3)*mb(3)))
         n2=0
         n3=1
         dk(1)=1773593.0d+00/253693440.0d+00
         rd(1)=0.882270011260d+00 
         ud(1)=0.140355381171d+00
         wd(1)=0.449332832327d+00
      else if(l.eq.23) then
         a1=(two**7)*73.0d+00/5242545.0d+00
         a2=(two**14)*1663.0d+00/(27.0d+00*121.0d+00*1458821.0d+00)
         a3=(three**10)*1599797.0d+00
     $      /((two**9)*(173.0d+00*173.0d+00)*(13.0d+00*13.0d+00)
     $        *6545.0d+00)
         n1=4
         bk(1)=5.51877146727d-03
         mb(1)=0.777493219315d+00
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         bk(2)=5.15823771181d-03
         mb(2)=0.912509096867d+00
         lb(2)=sqrt(half*(one-mb(2)*mb(2)))
         bk(3)=5.60870408259d-03
         mb(3)=0.314196994183d+00
         lb(3)=sqrt(half*(one-mb(3)*mb(3)))
         bk(4)=4.10677702817d-03
         mb(4)=0.982972302707d+00
         lb(4)=sqrt(half*(one-mb(4)*mb(4)))
         n2=1
         ck(1)=(38.0d+00**4)/
     $         ((33.0d+00**2)*(7.0d+00**3)*(1105.0d+00))
         pc(1)=0.938319218138d+00
         qc(1)=sqrt(one-pc(1)*pc(1))
         n3=1
c        note that there is a typo in lebedev's paper. he has the
c        following expression utilizing 2**6, instead of 2**9.
c        his explicit constant is ok.
         dk(1)=(19.0d+00**2)*(23.0d+00**6)/
     $         ((7.0d+00**3)*(two**9)*(9.0d+00)*(6113965.0d+00))
         rd(1)=0.836036015482d+00
         ud(1)=0.159041710538d+00
         wd(1)=0.525118572443d+00
      else if(l.eq.29) then
c        not yet debugged
         a1=8.54591172878d-04
         a2=zero
         a3=3.59911928502d-03
         n1=6
         bk(1)=3.65004580768d-03
         mb(1)=0.129238672710d+00
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         bk(2)=3.60482260142d-03
         mb(2)=0.371034178385d+00
         lb(2)=sqrt(half*(one-mb(2)*mb(2)))
         bk(3)=3.57672966173d-03
         mb(3)=0.743452042987d+00
         lb(3)=sqrt(half*(one-mb(3)*mb(3)))
         bk(4)=3.44978842429d-03
         mb(4)=0.867643624544d+00
         lb(4)=sqrt(half*(one-mb(4)*mb(4)))
         bk(5)=3.10895312238d-03
         mb(5)=0.949454317226d+00
         lb(5)=sqrt(half*(one-mb(5)*mb(5)))
         bk(6)=2.35210141366d-03
         mb(6)=0.990705621379d+00
         lb(6)=sqrt(half*(one-mb(6)*mb(6)))
         n2=2
         ck(1)=3.60082093222d-03
         pc(1)=0.820326419828d+00
         qc(1)=sqrt(one-pc(1)*pc(1))
         ck(2)=2.98234496317d-03
         pc(2)=0.964408914879d+00
         qc(2)=sqrt(one-pc(2)*pc(2))
         n3=2
         dk(1)=3.57154055427d-03
         rd(1)=0.251003475177d+00
         ud(1)=0.800072749407d+00
         wd(1)=0.544867737258d+00
         dk(2)=3.39231220501d-03
         rd(2)=0.902442529533d+00
         ud(2)=0.412772408317d+00
         wd(2)=0.123354853258d+00
      else if(l.eq.35) then
c        not yet tested.
         a1=0.526 597 657 614 280 65d-3
         a2=0.254 821 999 094 035 21d-2
         a3=0.251 231 737 094 410 58d-2
c
         n1=7
         lb(1)=0.690 934 631 051 134 58d0       
         mb(1)=0.212 646 822 756 572 27d0
         bk(1)=0.253 040 382 240 013 23d-2
c
         lb(2)=0.645 666 470 951 949 87
         mb(2)=0.407 712 664 234 151 23
         bk(2)=0.251 326 716 847 068 78d-2
c
         lb(3)=0.491 434 265 556 395 00
         mb(3)=0.719 016 498 610 493 14        
         bk(3)=0.250 172 512 106 477 33d-2
c
         lb(4)=0.392 725 982 232 176 49       
         mb(4)=0.831 584 394 850 904 11        
         bk(4)=0.244 537 330 479 967 86d-2
c
         lb(5)=0.286 128 917 876 582 18        
         mb(5)=0.914 472 790 579 114 05        
         bk(5)=0.230 269 443 256 207 58d-2
c
         lb(6)=0.177 483 652 423 745 68       
         mb(6)=0.967 987 141 569 894 03      
         bk(6)=0.201 427 826 095 260 94d-2
c
         lb(7)=0.075 680 958 662 444 68       
         mb(7)=0.994 255 895 125 528 88      
         bk(7)=0.146 249 508 154 751 42d-2
c
         n2=2
         pc(1)=0.977 642 808 920 987 23        
         qc(1)=0.210 272 533 073 347 57        
         ck(1)=0.191 095 131 473 050 82d-2
c
         pc(2)=0.881 813 289 360 544 12        
         qc(2)=0.471 598 688 194 885 97        
         ck(2)=0.241 744 235 754 198 47d-2
c
         n3=4
         rd(1)=0.099 217 699 713 625 76        
         ud(1)=0.334 436 316 957 483 71        
         wd(1)=0.937 180 984 636 078 86        
         dk(1)=0.223 660 770 713 642 63d-2
c
         rd(2)=0.205 482 371 254 664 95        
         ud(2)=0.450 233 038 742 967 35        
         wd(2)=0.868 946 031 654 344 86       
         dk(2)=0.241 693 001 073 811 79d-2
c
         rd(3)=0.106 801 825 135 337 23        
         ud(3)=0.590 515 703 098 041 30        
         wd(3)=0.799 927 855 836 003 99       
         dk(3)=0.251 223 686 473 367 06d-2
c
         rd(4)=0.310 428 403 275 151 30        
         ud(4)=0.555 015 236 814 480 68       
         wd(4)=0.771 746 262 280 424 63      
         dk(4)=0.249 664 405 192 924 56d-2
c
      else if(l.eq.41) then
c        not yet debugged; only 10 digits in lebedev.
         a1=0.3095121295d-03
         a2=zero
         a3=0.1852379698d-02
         n1=9
         bk(1)=0.9764331164d-03
         mb(1)=0.9962781297d+00
         lb(1)=sqrt(half*(one-mb(1)*mb(1)))
         bk(2)=0.1384737234d-02
         mb(2)=0.9784805837d+00
         lb(2)=sqrt(half*(one-mb(2)*mb(2)))
         bk(3)=0.1617210647d-02
         mb(3)=0.9414141582d+00
         lb(3)=sqrt(half*(one-mb(3)*mb(3)))
         bk(4)=0.1749564657d-02
         mb(4)=0.8830787279d+00
         lb(4)=sqrt(half*(one-mb(4)*mb(4)))
         bk(5)=0.1818471778d-02
         mb(5)=0.8028368773d+00
         lb(5)=sqrt(half*(one-mb(5)*mb(5)))
         bk(6)=0.1846715956d-02
         mb(6)=0.7007685753d+00
         lb(6)=sqrt(half*(one-mb(6)*mb(6)))
         bk(7)=0.1852028828d-02
         mb(7)=0.4333738687d+00
         lb(7)=sqrt(half*(one-mb(7)*mb(7)))
         bk(8)=0.1858812585d-02
         mb(8)=0.2703560883d+00
         lb(8)=sqrt(half*(one-mb(8)*mb(8)))
         bk(9)=0.1871790639d-02
         mb(9)=0.9219040707d-01
         lb(9)=sqrt(half*(one-mb(9)*mb(9)))
         n2=3
         ck(1)=0.1300321685d-02
         pc(1)=0.9850133350d+00
         qc(1)=sqrt(one-pc(1)*pc(1))
         ck(2)=0.1705153996d-02
         pc(2)=0.9180452877d+00
         qc(2)=sqrt(one-pc(2)*pc(2))
         ck(3)=0.1857161196d-02
         pc(3)=0.7911019296d+00
         qc(3)=sqrt(one-pc(3)*pc(3))
         n3=6
         dk(1)=0.1555213603d-02
         rd(1)=0.8213021581d-01
         ud(1)=0.2778673190d+00
         wd(1)=0.9571020743d+00
         dk(2)=0.1802239128d-02
         rd(2)=0.8999205842d-01
         ud(2)=0.5033564271d+00
         wd(2)=0.8593798558d+00
         dk(3)=0.1849830560d-02
         rd(3)=0.1816640840d+00
         ud(3)=0.5984126497d+00
         wd(3)=0.7803207424d+00
         dk(4)=0.1713904507d-02
         rd(4)=0.1720795225d+00
         ud(4)=0.3791035407d+00
         wd(4)=0.9092134750d+00
         dk(5)=0.1802658934d-02
         rd(5)=0.2634716655d+00
         ud(5)=0.4742392842d+00
         wd(5)=0.8400474883d+00
         dk(6)=0.1842866472d-02
         rd(6)=0.3518280927d+00
         ud(6)=0.5610263808d+00
         wd(6)=0.7493106119d+00
      else
         call plnkerr(' order requested not yet available',3200)
      endif
c
c     --- generate the points ---
      pts=1
c     --- a1 terms ---
      if(a1.ne.zero) then
         do 10 i=1,6
            w(i)=a1
   10    continue
         p(pts,1)  = one
         p(pts,2)  = zero
         p(pts,3)  = zero
c
         p(pts+1,1)= -one
         p(pts+1,2)= zero
         p(pts+1,3)= zero
c
         p(pts+2,1)= zero
         p(pts+2,2)= one
         p(pts+2,3)= zero
c
         p(pts+3,1)= zero
         p(pts+3,2)= -one
         p(pts+3,3)= zero
c
         p(pts+4,1)= zero
         p(pts+4,2)= zero
         p(pts+4,3)= one
c
         p(pts+5,1)= zero
         p(pts+5,2)= zero
         p(pts+5,3)=-one
c
         pts=pts+6
      endif
c
c     --- a2 terms ---
      root2i=one/sqrt(two)
      if(a2.ne.zero) then
         do 30 i=1,12
            w(pts+i-1)=a2
   30    continue   
         do 50 a=1,2
            signa=(-one)**(a+1)
            do 40 b=1,2
               signb=(-one)**(b+1)
               p(pts,1)  = signa*root2i
               p(pts,2)  = signb*root2i
               p(pts,3)  = zero
c
               p(pts+1,1)= signa*root2i
               p(pts+1,2)= zero
               p(pts+1,3)= signb*root2i
c
               p(pts+2,1)= zero
               p(pts+2,2)= signa*root2i
               p(pts+2,3)= signb*root2i
c
               pts=pts+3
   40       continue
   50    continue
      endif
c
c     --- a3 terms ---
      root3i=one/sqrt(three)
      if(a3.ne.zero) then
         do 60 i=1,8
            w(pts+i-1)=a3
   60    continue
         do 90 a=1,2
            signa=(-one)**(a+1) 
            do 80 b=1,2
               signb=(-one)**(b+1)
               do 70 c=1,2
                  signc=(-one)**(c+1)
                  p(pts,1) = signa*root3i
                  p(pts,2) = signb*root3i
                  p(pts,3) = signc*root3i
                  pts=pts+1
   70          continue
   80       continue
   90    continue
      endif
c
c     --- bk terms ---
      do 140 k=1,n1
         do 100 i=1,24
            w(pts+i-1)=bk(k)
  100    continue
         do 130 a=1,2
            signa=(-one)**(a+1)
            do 120 b=1,2
               signb=(-one)**(b+1)
               do 110 c=1,2
                  signc=(-one)**(c+1)
                  p(pts,1)  =signa*lb(k)
                  p(pts,2)  =signb*lb(k)
                  p(pts,3)  =signc*mb(k)
c
                  p(pts+1,1)=signa*lb(k)
                  p(pts+1,2)=signb*mb(k)
                  p(pts+1,3)=signc*lb(k)
c
                  p(pts+2,1)=signa*mb(k)
                  p(pts+2,2)=signb*lb(k)
                  p(pts+2,3)=signc*lb(k)
c
                  pts=pts+3
  110          continue
  120       continue
c
c
  130    continue
  140 continue
c
c     --- c(k) terms ---
      do 180 k=1,n2
         do 150 i=1,24
            w(pts+i-1)=ck(k)
  150    continue
         do 170 a=1,2
            signa=(-one)**(a+1)
            do 160 b=1,2
               signb=(-one)**(b+1)
               p(pts,1)  =signa*pc(k)
               p(pts,2)  =signb*qc(k)
               p(pts,3)  = zero
c
               p(pts+1,1)=signa*pc(k)
               p(pts+1,2)=zero
               p(pts+1,3)=signb*qc(k)
c
               p(pts+2,1)=zero
               p(pts+2,2)=signa*pc(k)
               p(pts+2,3)=signb*qc(k)
c
               p(pts+3,1)=signa*qc(k)
               p(pts+3,2)=signb*pc(k)
               p(pts+3,3)=zero
c
               p(pts+4,1)=signa*qc(k)
               p(pts+4,2)=zero
               p(pts+4,3)=signb*pc(k)
c
               p(pts+5,1)=zero
               p(pts+5,2)=signa*qc(k)
               p(pts+5,3)=signb*pc(k)
c
               pts=pts+6
  160       continue
  170    continue
  180 continue
c
c     --- d(k) terms ---
      do 230 k=1,n3
         do 190 i=1,48
            w(pts+i-1)=dk(k)
  190    continue
         do 220 a=1,2
            signa=(-one)**(a+1)
            do 210 b=1,2
               signb=(-one)**(b+1)
               do 200 c=1,2
                  signc=(-one)**(c+1)
                  p(pts,1)  =signa*rd(k)
                  p(pts,2)  =signb*ud(k)
                  p(pts,3)  =signc*wd(k)
c
                  p(pts+1,1)=signa*rd(k)
                  p(pts+1,2)=signb*wd(k)
                  p(pts+1,3)=signc*ud(k)
c
                  p(pts+2,1)=signa*ud(k)
                  p(pts+2,2)=signb*rd(k)
                  p(pts+2,3)=signc*wd(k)
c
                  p(pts+3,1)=signa*ud(k)
                  p(pts+3,2)=signb*wd(k)
                  p(pts+3,3)=signc*rd(k)
c
                  p(pts+4,1)=signa*wd(k)
                  p(pts+4,2)=signb*ud(k)
                  p(pts+4,3)=signc*rd(k)
c
                  p(pts+5,1)=signa*wd(k)
                  p(pts+5,2)=signb*rd(k)
                  p(pts+5,3)=signc*ud(k)
c
                  pts=pts+6
  200          continue
  210       continue
  220    continue
  230 continue
c
c     --- check that number of points processed is same as expected ---
      pts=pts-1
      if(npts.ne.pts) then
         write(iout,1000) npts,pts
         call plnkerr('wrong number of points passed to lebedev',3201)
      endif
c
c     --- scale weights by 4pi ---
      pi=4.0d+00*atan(1.0d+00)
      do 250 i=1,npts
         w(i)=four*pi*w(i)
  250 continue
c
c
 1000 format(1x,'expected',i4,' generated',i4)
c
c
      return
      end
