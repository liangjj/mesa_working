*deck %W%  %G%
c***begin prologue     lebdev
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lebdev, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            lebedev quadrature points and weights
c
c***routines called
c***end prologue       lebdev
      subroutine lebdev(pt,angle,cthet,sthet,phpt,sphi,cphi,wt,nleb,
     1                  npts,str)
      implicit integer (a-z)
      real*8 fourpi, cthet, sthet, phpt, sphi, cphi, pt, wt, l, m
      real*8 biga, bigb, bigc, bigd, p, q, r, u, w
      real*8 sqrt2, sqrt3, sumwt, angle
      character *(*) str
      dimension cthet(*), sthet(*), phpt(*), sphi(*), cphi(*)
      dimension pt(3,npts), wt(*), l(4), m(4), biga(4)
      dimension bigb(4), bigc(4), bigd(4), p(1), q(1)
      dimension r(1), u(1), w(1), angle(npts,2)
      common /io/ inp, iout
      data sqrt2 /.70710678118654752440d0 /
      data sqrt3 /.57735026918962576451d0 /
      data fourpi /  12.566370614359172d0  /
c              the a1 nodes-these are independent of the quadrature      
       call xyz(pt(1,1),0.d0,0.d0,1.d0)
       call xyz(pt(1,2),0.d0,0.d0,-1.d0)
       call xyz(pt(1,3),0.d0,1.d0,0.d0)
       call xyz(pt(1,4),0.d0,-1.d0,0.d0)
       call xyz(pt(1,5),1.d0,0.d0,0.d0)
       call xyz(pt(1,6),-1.d0,0.d0,0.d0)
c              the a2 nodes-these are independent of the quadrature
      if (nleb.ne.17) then
          call xyz(pt(1,7),sqrt2,sqrt2,0.d0)
          call xyz(pt(1,8),sqrt2,-sqrt2,0.d0)
          call xyz(pt(1,9),-sqrt2,sqrt2,0.d0)
          call xyz(pt(1,10),-sqrt2,-sqrt2,0.d0)
          call xyz(pt(1,11),sqrt2,0.d0,sqrt2)
          call xyz(pt(1,12),sqrt2,0.d0,-sqrt2)
          call xyz(pt(1,13),-sqrt2,0.d0,sqrt2)
          call xyz(pt(1,14),-sqrt2,0.d0,-sqrt2)
          call xyz(pt(1,15),0.d0,sqrt2,sqrt2)
          call xyz(pt(1,16),0.d0,sqrt2,-sqrt2)
          call xyz(pt(1,17),0.d0,-sqrt2,sqrt2)
          call xyz(pt(1,18),0.d0,-sqrt2,-sqrt2)
          ii=18
      else
          ii=6          
      endif
c              the a3 nodes-these are independent of the quadrature
      call xyz(pt(1,ii+1),sqrt3,sqrt3,sqrt3)
      call xyz(pt(1,ii+2),sqrt3,sqrt3,-sqrt3)
      call xyz(pt(1,ii+3),sqrt3,-sqrt3,sqrt3)
      call xyz(pt(1,ii+4),sqrt3,-sqrt3,-sqrt3)
      call xyz(pt(1,ii+5),-sqrt3,sqrt3,sqrt3)
      call xyz(pt(1,ii+6),-sqrt3,sqrt3,-sqrt3)
      call xyz(pt(1,ii+7),-sqrt3,-sqrt3,sqrt3) 
      call xyz(pt(1,ii+8),-sqrt3,-sqrt3,-sqrt3)
      do 10 i=1,4
          m(i)=0.d0
          l(i)=0.d0
   10 continue
      if (nleb.eq.11) then
          m(1)=.904534033733d0
          l(1)=.301511344578d0
          bigc(1)=0.d0
          biga(1)=.0126984126985d0
          biga(2)=.0225749559083d0
          biga(3)=.0210937500000d0
          bigb(1)=.0201733355379d0
          bigb(2)=0.d0
          bigb(3)=0.d0
          bigb(4)=0.d0
          do 20 i=1,6
            wt(i)=biga(1)
            wt(i+6)=biga(2)
            wt(i+12)=biga(2)
            wt(i+18)=biga(3)
   20     continue
          wt(25)=biga(3)
          wt(26)=biga(3)
          call getbnd(pt(1,ii+9),wt(ii+9),l,m,bigb,nleb)
          npts=50
      elseif (nleb.eq.17) then
          n1=3
          bigc(1)=.00969499636166d0
          biga(1)=.00382827049494d0
          biga(2)=0.d0
          biga(3)=.00979373751249d0
          bigb(1)=.00821173728319d0
          bigb(2)=.00959547133607d0
          bigb(3)=.00994281489118d0
          bigb(4)=0.d0
          m(1)=.965124035087d0
          l(1)=.185115635345d0
          m(2)=.828769981253d0
          l(2)=.395689473056d0
          m(3)=.215957291846d0
          l(3)=.690421048382d0
          p(1)=.878158910604d0
          q(1)=.478369028812d0
          do 30 i=1,6
             wt(i)=biga(1)
   30     continue
          do 40 i=7,14    
             wt(i)=biga(3)
   40     continue     
          call getbnd(pt(1,ii+9),wt(ii+9),l,m,bigb,nleb)
          call getcnd(pt(1,ii+81),wt(ii+81),p,q,bigc)
          npts=110
       elseif (nleb.eq.23) then
          n1=4
          n2=1
          n3=1
          bigc(1)=5.05184606462d-03
          p(1)=.938319218138d0
          q(1)=.345770219761d0
          bigd(1)=5.53024891623d-03
c         --- fixing typo ala bis instructions
c         r(1)=.859059015482d0
          r(1)=.8360036015482d0
          u(1)=.159041710538d0
          w(1)=.525118572443d0
          biga(1)=1.78234044724d-03
          biga(2)=5.71690594998d-03
          biga(3)=5.57338317884d-03
          bigb(1)=5.51877146727d-03
          bigb(2)=5.15823771181d-03
          bigb(3)=5.60870408259d-03
          bigb(4)=4.10677702817d-03
          m(1)=.777493219315d0
          l(1)=.444693317871d0
          m(2)=.912509096867d0
          l(2)=.289246562758d0
          m(3)=.314196994183d0
          l(3)=.671297344270d0 
          m(4)=.982972302707d0
          l(4)=.129933544765d0
          do 50 i=1,6
             wt(i)=biga(1)
   50     continue
          do 60 i=7,18    
             wt(i)=biga(2)
   60     continue     
          do 70 i=19,26
             wt(i)=biga(3)
   70     continue          
          call getbnd(pt(1,ii+9),wt(ii+9),l,m,bigb,nleb)
          call getcnd(pt(1,ii+105),wt(ii+105),p,q,bigc)
          call getdnd(pt(1,ii+129),wt(ii+129),r,u,w,bigd)
          npts=194          
      endif
      call ang(pt,cthet,sthet,cphi,sphi,phpt,npts)
      sumwt=0.d0
      do 80 i=1,npts
         wt(i)=wt(i)*fourpi 
         sumwt=sumwt+wt(i)
   80 continue
      sumwt=sumwt/fourpi
      write(iout,1) sumwt
    1 format (/,' sum of lebedev weights = ',e15.8)
      do 90 i=1,npts
         angle(i,1)=cthet(i)
         angle(i,2)=phpt(i)
   90 continue
      call iosys ('write real "lebedev angular points '//
     1             str//'" to lamdat',2*npts,angle,0,' ')
      call iosys ('write real "lebedev angular weights '//
     1             str//'" to lamdat',npts,wt,0,' ')                          
      return
      end
