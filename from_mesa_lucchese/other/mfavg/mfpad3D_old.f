      program mfpad
c
c  Body-Frame photoionization angular distributions
c    modified from Shungo Myabe's version in summer of 2011 by Trevisan and McCurdy to 
c    (0)  Write out a 3D mfpad on a theta, phi grid (Trevisan 2010)
c    (1)  compute integral over polarization directions
c    (2)  compute integral over photoelectron directions
c    (3)  correct the phase inconsistency in the definition of the dipole operator
c           (see comment below)
c
c   Note:  for this code to calculate a 3D mfpad correctly, the input must contain
c   all three cartesian components of the dipole operator whether they live in 
c   different symmetries (irreps) or in the same irrep, as when doing a calculation
c   in no symmetry.
c
      implicit real*8(a-h,o-z)
      parameter (ltop=10,mutop=10,ith=500)
      parameter (nsymx=8,nlmmax=72)
      parameter(nenemax=100,nchanmx=10,nrmx=50,npmx=64*nrmx)
      dimension totc(nsymx,nenemax)
      real*8 totan(nenemax),kchan(nchanmx)
      dimension th(ith), ppe(ith)
      character*8 name(nsymx)
c
      dimension nlm(nsymx),lch(nsymx,nlmmax),mch(nsymx,nlmmax)
      dimension istart(nsymx),istop(nsymx),jstart(nsymx),jstop(nsymx)
      dimension ylm(0:ltop,-mutop:mutop)
      dimension ylmp(0:ltop,-mutop:mutop)
      dimension p(0:ltop),pp(0:ltop)
      dimension xs(ith,ith,nenemax),egy(nenemax),twoe(nenemax)
      dimension xspol(ith,ith,nenemax)
      dimension xsint(ith,ith,nenemax)
      dimension xsec(nenemax)
      dimension fac(0:60)
      dimension cspl(ith+3),dspl(ith+3)
      complex*16 ai,cmp, tmat(nlmmax,nlmmax), f(nenemax)
      complex*16  polariz_component(-1:1,nenemax) 
      complex*16  polariz_direction(nlmmax,nsymx,nenemax) 
      real*8  pol_integral(nenemax)
      real*8  photoelect_integral(nenemax) 
      ai=(0.,1.)
      data pi / 3.141592653589793d0  /
c
c  input of quadrature info, scattering angles on unit5
c  input of channel  info on unit8
c  input of "T matrices" (dipole matrix elements)
c    on unit8, this unit is rewound and reread many times.
c
c      call link("unit5=inpdcs,unit6=(outdcs,hc,create),
c     1 unit7=(plthead,hc,create)//")
c
      open(5,file='inmfpad',status='unknown')
      open(6,file='outmfpad',status='unknown')
c  unit 7 is for 3D mfpat output
      open(7,file='pltmfpad3d',status='unknown')
c  unit 97 is for mfpad integrated over polarization, differential in photoelectron direction
      open(97,file='intmfpad',status='unknown')
c  unit 98 is for mfpad integrated over photoelectron directions, differential in polarization
      open(98,file='polarization_mfpad',status='unknown')
c
      fac(0)=1.
      lmmax = 60
      do 61 i=1,lmmax
 61   fac(i)=i*fac(i-1)
c
c
      read(5,*) lmax,mumax,iprt
      write(6,158) lmax,mumax
158   format( ' maximum l and m are ',2i5)
c read symmetry and channel information
      read(5,*) nsym,nchan,ichan,jchan
c ichan must be one; it refers to the polarization vector
      if(ichan.ne.1)then
         write(6,*)"stopping because ichan is not 1"
         stop
      endif
      write(6,77)nchan,nsym,ichan,jchan
77    format(//" this is a",i3," channel calculation with",i3," symmetry
     1 components"/" cross sections for going from channel",i3,
     1 " to channel",i3)
      do 23 i=1,nsym
         read(5,750)name(i)
 750     format(a8)
         write(6,78)name(i)
 78      format(" symmetry  ",a8," has the following l,m pairs:")
c         read(5,*) nlm(i),sgni(i)
         read(5,*) nlm(i)
         read(5,*) (lch(i,j),mch(i,j),j=1,nlm(i))
         write(6,79) (lch(i,j),mch(i,j),j=1,nlm(i))
 79      format(2i5)
         read(5,*) istart(i),istop(i),jstart(i),jstop(i)
         write(6,76) istart(i),istop(i),jstart(i),jstop(i)
 76    format(" starting and stopping ranges for channels i and j are:",
     1        4i4)
23    continue
      read(5,*) nenergy,polang,phi
c modified to plot out 3D mfpad (Shungo, 8/19/09)
c Setting the polarization
c  polrad= polarization theta in radians
c  phirad= polarization phi   in radians
      polrad=polang*pi/180.
      spol=sin(polrad)
      cpol=cos(polrad)
      phirad=phi*pi/180.
      cpk=cos(phirad)
      spk=sin(phirad)
c number of points, ipts for theta, ipts1 for phi
      ipts=30
      ipts1=60
c Starting loop on phi of photo-e => perad, ppe(phi); (0~2Pi) ###################
      do 888 iphi=1,ipts1
        pephi=((dble(iphi)/dble(ipts1))+0.0001d0)*359.9999d0
        ppe(iphi)= pephi*pi/180.
        perad=pephi*pi/180.
        cpkp=cos(perad)
        spkp=sin(perad)
c Starting loop on theta of photo-e => theta, th(ipts); (0~Pi) ##################
        do 1000 ithe=1,ipts
          theta=((dble(ithe)/dble(ipts))+0.0001d0)*179.9999d0
          th(ithe)= theta*pi/180.
c  ***********fix for photelectron polar angles bigger than 180
c CWM: this mysterious little bit of code is necessary, but it isn't clear why
c  the theta angle can be ge. 180.
          if(theta.ge.180.)then
            theta=360.-theta
            phin=pephi+180.
            perad =phin*pi/180.
            cpkp=cos(perad)
            spkp=sin(perad)
          endif
c  ***********
          theta = theta*pi/180.
          s=sin(theta)
          c=cos(theta)
          do 70 izero=1,nenemax
70          xsec(izero) = 0.0
c
c Now we compute Ylm's 
c
            do 13 mu=0,mumax
              call plm(cpol,narg,mu,lmax,p)
              call plm(c,narg,mu,lmax,pp)
              if(mu.eq.0)then
               do 10 i=0,lmax
                 const=sqrt((2*i+1)/(4.*pi))
                 ylmp(i,0) =  pp(i)*const
                 ylm(i,0)=p(i)*const
 10            continue
              else
                 cmp=(cmplx(cpk,spk))**mu
                 cmphi=dble(cmp)
                 smphi=imag(cmp)
                 cmp=(cmplx(cpkp,spkp))**mu
                 cmphip=dble(cmp)
                 smphip=imag(cmp)
               do 12 i=mu,lmax
                 const=sqrt((2*i+1)/(4.*pi)*fac(i-mu)/fac(i+mu))
c  sqrt(2) factor added to normalize "real valued" ylms
                 const=const*sqrt(2.0)
                 ylm(i,-mu)=p(i)*const*cmphi
                 ylm(i,mu)=p(i)*const*smphi
                 ylmp(i,-mu)=pp(i)*const*cmphip
                 ylmp(i,mu)=pp(i)*const*smphip
 12            continue
              endif
 13         continue
c
c we have the ylm's. Now read "T matrices" = partial wave dipole amplitudes
c  including the Coulomb phases, and multiplicative constants in the original
c  papers on photoionization
c   T. N. Rescigno, B. H. Lengsfield, and A. E. Orel,J. Chem. Phys. 99, 5097 (1993). 
c
c***********note**************
c we are only going to use the tmatrix block corresponding to a
c specified channel pair.
c****************************
c Zero arrays in the amplitudes will be stored
      do 49 izero=1,nenergy
c
c  Set the accumulator for the contributions to 
c  the three components of the polarization
c  to zero for use in calculating the integral over polarization 
c  directions and over photoelectron directions
         polariz_component(-1,izero) = (0.d0,0.d0)
         polariz_component(0,izero) = (0.d0,0.d0)
         polariz_component(+1,izero) = (0.d0,0.d0)
         do jsym = 1,nsymx
          do jlm =1, nlmmax
            polariz_direction(jlm,jsym,izero) = (0.d0,0.d0)
          enddo
         enddo
c
49    f(izero) = (0.0,0.0)
      do 25 isym=1,nsym
      if(ithe.eq.1.and.iphi.eq.1)then
      write(6,"(a8)")name(isym)
      endif
      open(8,file=name(isym),form='unformatted',status='unknown')
c
      rewind(8)
      do 24 iene=1,nenergy
c
66    continue
      read(8)icw,jcw,ni,nj,aki,akj
      twoe(iene)=akj*akj
      read(8)((tmat(i,j),i=1,ni),j=1,nj)
      if(ichan.eq.icw.and.jchan.eq.jcw)go to 67
      go to 66
67    continue
222   format(4x,6e12.5)
c
c calculate total cross sections on first pass through these loops
c
      if(ithe.eq.1.and.iphi.eq.1)then
c
c  print Tmatrix if requested
         if(iprt.gt.0) then
            write(6,*)' Tmatrix ',icw,jcw
            write(6,88001)((tmat(i,j),i=1,ni),j=1,nj)
88001       format(3(4x,f12.8,2x,f12.8))
         end if 
c  write total cross section contributions for each polarization direction
c  if iprt.gt.0
         if(iprt.gt.0) then
           do i=1,ni
               sum=0.d0
            do j=1,nj
               sum = sum + abs(tmat(i,j))**2
            enddo
            write(6,'("sum of |Tij|^2 for i =",i3, " = ",e12.5)')i,sum
           enddo
          endif
c 
c..photo
c
         totc(isym,iene)=0.
         ifirst=istart(isym)
         ilast=istop(isym)
         jfirst=jstart(isym)
         jlast=jstop(isym)
         do 71 i=ifirst,ilast
            ii=i+1-ifirst
            do 71 j=jfirst,jlast
               jj=j+1-jfirst
               totc(isym,iene)=totc(isym,iene)+abs(tmat(ii,jj))**2
 71      continue
c
         if(iprt.gt.0) then
            write(6,*)' isym iene totc ',isym,iene,totc(isym,iene)
         end if
c
      end if
c
c
c accumulate scattering amplitude for a particular energy
c
      ifirst=istart(isym)
      ilast=istop(isym)
      jfirst=jstart(isym)
      jlast=jstop(isym)
c  Loops to accumulate (1) 3D MFPAD for fixed input polarization direction
c   and (2) the separate contributions to each polarization component 
c   to cross section differential in photoelectron directions
c   and (3) contributions to cross section differential only in polarization direction
        do 31 i=ifirst,ilast
         ii=i+1-ifirst
         li = lch(isym,i)
         mi = mch(isym,i)
c  Correct inconsistency in definition of phase of Ylms when being 
c  used to specify polarization direction -- CWM and CST 2011
c  Cartesian x,y,z were used in cphot, but real Ylms are used here, hence phase difference
c*tnr* this was done in the code cdipc.f, so it's being commented out here
           phase = 1.d0
c           if(mi.eq.1.or.mi.eq.-1) phase = -1.d0
         do 32 j=jfirst,jlast
          jj=j+1-jfirst
          lj = lch(isym,j)
          mj = mch(isym,j)
          f(iene)=ylm(li,mi)*phase*tmat(ii,jj)*ylmp(lj,mj)
     #           +f(iene)
c
c accumulate contributions to each component of the polarization
c mi goes from -1 to 1, and this sum of amplitudes is used 
c to calculate the integral over polarizations: |Mx|^2 + |My|^2 + |Mz|^2
c NOTE: the phase correction to the ylm(li,mi), the variable  "phase" 
c  above doesn't matter here.
c
         if(mi.lt.-1.or.mi.gt.1) then
          write(6,*) 'polarization component has |m|>1. input error?'
            stop
          else
            polariz_component(mi,iene) = tmat(ii,jj)*ylmp(lj,mj) 
     #          + polariz_component(mi,iene) 
          endif
c
c  accumulate contributions to cross section differential in polarization directions
c   NOTE ylmp(l,m) are Ylm's that are evaluated on a theta, phi grid -- originally
c   for directions of the photoelectron.  ylm(l,m)  are Ylm's that are evaluated only
c   in the direction of polarization input above.  Here we use ylmp(li,mi) to get
c   the cross section for integrated electron directions, with polarization directions
c   on the theta, phi grid. Here the phase correction to ylm(li,mi) (here ylmp(li,mi) ) matters!
c
         if(mi.lt.-1.or.mi.gt.1) then
          write(6,*) 'polarization component has |m|>1. input error?'
            stop
          else
           polariz_direction(jj,isym,iene) =  
     #        ylmp(li,mi)*phase*tmat(ii,jj)  
     #        + polariz_direction(jj,isym,iene)
          endif
c
c  end loop on j index (electronic angular momentum lm pairs)
 32     continue
c  end loop on i index (dipole operator reduced tensor components)
 31    continue
c end loop on energies 
24    continue
      close(8)
c end loop on molecular symmetries
25    continue
c
      do 50 iene=1,nenergy
c
c accumulate cross section(mfpad) and the integral of the mfpad over polarization
c directions 
c
      xsec(iene)=abs(f(iene))**2
c
      pol_integral(iene) =  abs(polariz_component(-1,iene))**2 +
     #              abs(polariz_component(0,iene))**2 +
     #             abs(polariz_component(+1,iene))**2
c accumulate the integral of the mfpad over photoelectron directions (differential in 
c polarization directions)
c the loops over energy symmetry and lm pairs are closed and we have accumulated
c accumulated the sum over mu of T(mu,l0,m0)*Y_1,mu*phase for each lm pair and energy
c  at this grid point.  Each of those sums adds incoherently now to make this
c  mfpad differential in polarization only.
      photoelect_integral(iene) = 0.d0
      Do jsym = 1,nsym
        jfirst=jstart(jsym)
        jlast=jstop(jsym)
        Do  j = jfirst, jlast
          jj=j+1-jfirst
          photoelect_integral(iene) = 
     #        abs(polariz_direction(jj,jsym,iene))**2 +
     #        photoelect_integral(iene)
         enddo
       enddo
50    continue
      write(6,157) th(ithe), ppe(iphi)
157   format(' scattering angle (theta, phi)  =',2 f15.8)
      do 80 iene=1,nenergy
c factors of 2*pi here SEEM to have to do with normalization of phi part of the real Ylm's
c as defined above. -- CWM
         xsec(iene)=xsec(iene)/(2.*pi)
         pol_integral(iene) = pol_integral(iene)/(2.*pi)
         photoelect_integral(iene) = photoelect_integral(iene)/(2.*pi)
         einc = twoe(iene)/2.
         write(6,101) einc,xsec(iene)
 101     format(' incident e =',e12.5,'  MFPAD =',e12.5)
         egy(iene)=einc
         xs(iphi,ithe,iene)=xsec(iene)
         xsint(iphi,ithe,iene)=pol_integral(iene)
         xspol(iphi,ithe,iene) = photoelect_integral(iene) 
80    continue
c ending loops on photoelectron theta and phi ###########
1000  enddo
 888  enddo
c ... and write out the results (set up for gnuplot)
      do iene=1,nenergy
        write(7,*) 'Energy:', egy(iene)
        do ithe=1,ipts
          do iphi=1,ipts1
            write(7,*) xs(iphi,ithe,iene), th(ithe), ppe(iphi)
          enddo
          write(7,*)
        enddo
      enddo
c ... now write out mfpads integrated over all polarization
c     directions. 
      do iene=1,nenergy
        write(97,*) 'Energy:', egy(iene)
        do ithe=1,ipts
          do iphi=1,ipts1
            write(97,*) xsint(iphi,ithe,iene), th(ithe), ppe(iphi)
          enddo
          write(97,*)
        enddo
      enddo
c ... now write out mfpads integrated over all photoelectron
c     directions. 
      do iene=1,nenergy
        write(98,*) 'Energy:', egy(iene)
        do ithe=1,ipts
          do iphi=1,ipts1
            write(98,*) xspol(iphi,ithe,iene), th(ithe), ppe(iphi)
          enddo
          write(98,*)
        enddo
      enddo
c  total cross sections 
      continue
      do 409 i=1,nenergy
      totan(i)=0.
      do 409 j=1,nsym
c..photo
      totan(i)=totan(i)+totc(j,i)
c..photo
409   continue
      scalbhl=.529*.529
      do 997 i=1,nenergy
      anal=totan(i)
c
      write(6,764)egy(i),anal,scalbhl*anal
764   format(" energy =",f10.5,"  total cross section=",e12.4,
     $ " total (angstroms**2) ",e12.4)
997   continue
998   format(e16.8,/,(2f20.10))
      call exit
      end
      function spline1(x,y,z,nn,c,d,isw)
      implicit real*8(a-h,o-z)
      character*8 error(2)
      dimension x(1),y(1),c(*),d(*)
    1 format("** error in spline1...",2a8/)
      error(1)="unordere"
      error(2)="d x-vals"
    2 n=nn
      np1=n+1
      np2=n+2
      go to (4,19),isw
    3 return
c     section 4 calculates the spline coefficients c(1) is the
c     constant term, c(2) is the coeff of the linear term, and
c     c(3) thru c(n+2) are the spline coefficients.
    4 c(1)=y(1)
      d(1)=1.0
      c(np1)=0.0
      d(np1)=0.0
      c(np2)=0.0
      d(np2)=0.0
      do 5 i=2,n
      c(i)=y(i)-y(1)
    5 d(i)=x(i)-x(1)
      do 14 i=3,np2
      pivot=1.0/d(i-1)
      if(i.lt.np2)then
         supd=x(i-1)-x(i-2)
         if(supd.lt.0.0)then
            write(6,1) error(1),error(2)
            spline1=0.0
            return
         else
            supd=supd**3
            go to 10
         endif
      else
         supd=1.0
      endif
   10 dfact=supd*pivot
      cfact=c(i-1)*pivot
      if(i.gt.n) go to 30
      do 11 j=i,n
      v=x(j)-x(i-2)
      c(j)=c(j)-d(j)*cfact
   11 d(j)=(v**3)-d(j)*dfact
   30 if(i.lt.np2) then
         c(np1)=c(np1)-d(np1)*cfact
         d(np1)=1.0-d(np1)*dfact
      endif
      c(np2)=c(np2)-d(np2)*cfact
 14   d(np2)=x(i-2)-d(np2)*dfact
      do 18 i=1,n
      j=np2-i
      if(j.eq.np1) then 
 15      v=1.0
         go to 17
      else
 16      v=x(j)-x(j-1)
         v=v**3
      endif
   17 c(j+1)=c(j+1)/d(j+1)
   18 c(j)=c(j)-c(j+1)*v
      c(2)=c(2)/d(2)
      go to 22
c     section 19 evaluates the spline function
   19 spline1=c(1)+c(2)*(z-x(1))
      do 20 i=1,n
      v=z-x(i)
      if(v.le.0.0) return
 20      spline1=spline1+c(i+2)*(v**3)
   22 spline1=0.0
      return
      end
      subroutine plm(x,n,mu,lmax,p)
      implicit real*8(a-h,o-z)
      dimension p(0:*)
c
c associated legengre polynomials P(l,mu) .
c for fixed mu and l=mu,,,lmax 
c
c
c mu = 0 case
c
      if(mu.eq.0)then

       p(0)=1.
2      p(1)=x
       if(lmax.lt.2)return
       do 1 l=2,lmax
3      p(l)=((2*l-1)*x*p(l-1)-(l-1)*p(l-2))/l
1      continue
       return
      endif
c
c mu = 1 case
c
      if(mu.eq.1)then
       arg=sqrt(1.-x*x)
6      p(1)=-arg
       if(lmax.lt.2)return
7      p(2)=-3.*x*arg
       if(lmax.lt.3)return
       do 8 l=3,lmax
8      p(l)=((2*l-1)*x*p(l-1)-l*p(l-2))/(l-1)
       return
      endif
c
c mu must be larger than 1
c
c
c recurr across to l=mu+1, with m=0
c
      mu1=mu+1
      p0=1.
4     p1=x
      do 5 l=2,mu1
      p2=((2*l-1)*x*p1-(l-1)*p0)/l
      p0=p1
      p1=p2
5     continue
c
c recurr across to l=mu+1, with m=1
c
      arg=sqrt(1.-x*x)
      q0=-arg
9     q1=-3.*arg*x
      do 10 l=3,mu1
      q2=((2*l-1)*x*q1-l*q0)/(l-1)
      q0=q1
10    q1=q2
c
c with l fixed at mu and mu+1, recurr down to m=mu
c
      l0=mu
      l1=mu1
      do 11 m=2,mu
      p2=-2.*(m-1)*x*q0/arg-(l0-m+2)*(l0+m-1)*p0
      q2=-2.*(m-1)*x*q1/arg-(l1-m+2)*(l1+m-1)*p1
      p0=q0
      p1=q1
      q0=p2
11    q1=q2
c
c compute the desired vector and quit
c
      p(mu)=q0
12    continue
      if(lmax.eq.mu)return
13    p(mu1)=q1
      if(lmax.eq.mu1)return
      mu2=mu1+1
      do 14 l=mu2,lmax
14    p(l)=((2*l-1)*x*p(l-1)-(l+mu-1)*p(l-2))/(l-mu)
      return
      end
