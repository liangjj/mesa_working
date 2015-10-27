*deck fbnint
c***begin prologue     fbnint
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            greens function quadrature
c***
c***description        performs the integration of the radial
c***                   greens function and an inhomogeneity using
c***                   newton-cotes quadrature               
c***references
c
c***routines called
c
c***end prologue       fbnint
c
      subroutine fbnint(flm,psilm,r,j,y,wt,scr,lmax,m,n,nshell,nr,
     1                  nints,nwts,prnt)
      implicit integer (a-z)
      real*8 flm, psilm, wt, j, y, scr, int0f, int0b, r, tol, test
      real*8 delta
      logical prnt, multip
      dimension flm(n,m:lmax), wt(nwts), psilm(nints,m:lmax)
      dimension j(n,0:lmax), y(n,0:lmax), nr(nshell), scr(nints,4)
      dimensions r(n)
      parameter ( tol=1.d-10 )
      common /io/ inp, iout
c----------------------------------------------------------------------c     
c                when we are done we have in psilm                     c
c                                                                      c
c            psilm(i) = y(i) * sum ( k=1..i ) j(k)*flm(k)*wt(k) +      c
c                       j(i) * sum ( k=i+1...n ) y(k)*flm(k)*wt(k)     c
c            which is the radial wavefunction at all of the quadrature c
c            points not including the origin.                          c
c            thus we have constructed a quadrature approximation to;   c
c                                                                      c
c            psilm(x) = y(x) * integral(0,x) [ j(x) * flm(x) ] +       c
c                       j(x) * integral(x,infinity) [ y(x)*flm(x0 ]    c
c                                                                      c
c----------------------------------------------------------------------c
      do 10 l=m,lmax
c----------------------------------------------------------------------c
c             initialize the integrals and counters at first and last  c
c                                 points                               c
c----------------------------------------------------------------------c      
         int0f=0.d0
         int0b=0.d0
         locf=1
         locb=n+1
         locwtf=1
         locwtb=nwts+1
         locntf=1
         locntb=nints+1
         multip=.false.
         nb=nshell
         first=0
         call rzero(scr,4*nints)
         call rzero(psilm(1,l),n)
         do 20 ns=1,nshell
            locb=locb-nr(nb)
            locwtb=locwtb-nr(nb)*(nr(nb)-1)
            locntb=locntb-nr(nb)+1
            if (.not.multip) then
                 test=int0f 
                 call forwrd(flm(locf,l),j(locf,l),y(locf,l),
     1                       wt(locwtf),int0f,scr(locntf,1),
     2                       scr(locntf,2),nr(ns))
                 call bkwrd(flm(locb,l),j(locb,l),y(locb,l),wt(locwtb),
     1                      int0b,scr(locntb,3),scr(locntb,4),nr(nb))
                 delta=abs(test-int0f)
                 if (delta.le.tol) then
                     multip=.true.
                     if (first.eq.0) then
                         write(iout,3) l, ns
                         first=1
                     endif        
                 endif
             else
                 call appint(y(locf,l),int0f,scr(locntf,1),
     1                       scr(locntf,2),nr(ns))
             endif        
            locf=locf+nr(ns)
            locwtf=locwtf+nr(ns)*(nr(ns)-1)
            locntf=locntf+nr(ns)-1
            nb=nb-1
  20     continue
         if (prnt) then
             write(iout,1) m, l, (scr(i,2),i=1,nints)
         endif    
         if (prnt) then
             write(iout,2) m, l, (scr(i,4),i=nints,1,-1)
         endif   
         do 40 int=1,nints-1
            psilm(int+1,l)=scr(int,1)+scr(int+1,3)
   40    continue
   10 continue
      return
    1 format(/,'forward partial integral for m = ',i3,1x,'l = ',i3,
     1          (/,5(e15.8,1x)))
    2 format(/,'backward partial integral for m = ',i3,1x,'l = ',i3,
     1          (/,5(e15.8,1x)))
    3 format(/,1x,'l = ',i3,1x,'shell = ',i3,1x,'going asymptotic')     
      end















