      program mfavg
      implicit real*8 (a-h,o-z)
      parameter (ltop=10,mutop=10,ith=500)
      parameter (nsymx=8,nlmmax=72)
      parameter(nenemax=100,nchanmx=10,nrmx=50,npmx=64*nrmx)
      dimension totc(nsymx,nenemax)
      real*8 totan(nenemax),kchan(nchanmx)
      dimension th(ith), ppe(ith)
      character*8 name(nsymx)
c
      dimension sgni(nchanmx)
      dimension nlm(nsymx),lch(nsymx,nlmmax),mch(nsymx,nlmmax)
      dimension istart(nsymx),istop(nsymx),jstart(nsymx),jstop(nsymx)
      dimension ylm(0:ltop,-mutop:mutop)
      dimension ylmp(0:ltop,-mutop:mutop)
      dimension p(0:ltop),pp(0:ltop)
      dimension xs(ith,ith,nenemax),egy(nenemax),twoe(nenemax)
      dimension xsec(nenemax)
      dimension fac(0:60)
      dimension cspl(ith+3),dspl(ith+3)
      complex*16 ai,cmp, tmat(nlmmax,nlmmax), f(nenemax)
      ai=(0.,1.)
      data pi / 3.141592653589793d0  /
      twopi=2.d0*pi
      open(5,file='inmfpad',status='unknown')
      open(6,file='outmfpad',status='unknown')
c  unit 7 is for 3D mfpat output
      open(7,file='pltmfpad3d',status='unknown')
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
         read(5,*) nlm(i),sgni(i)
c         read(5,*) nlm(i)
         read(5,*) (lch(i,j),mch(i,j),j=1,nlm(i))
         write(6,79) (lch(i,j),mch(i,j),j=1,nlm(i))
 79      format(2i5)
         read(5,*) istart(i),istop(i),jstart(i),jstop(i)
         write(6,76) istart(i),istop(i),jstart(i),jstop(i)
 76    format(" starting and stopping ranges for channels i and j are:",
     1        4i4)
23    continue
      read(5,*) nenergy,polang,phi
      write(6,*)'polang,phi',polang,phi
c polarization angles are here defined in the body-(ie MESA) frame, 
c not necessarily the laboatory z-axis. Here, it is assumed to be the x-axis.
c  polrad= polarization theta in radians
c  phirad= polarization phi   in radians
      polrad=polang*pi/180.
      phirad=phi*pi/180.
c Get the polarization vector cartesian coordinates
      xpol=cos(polrad)
      zpol=sin(polrad)*cos(phirad)            
      ypol=sin(polrad)*sin(phirad)
c now get the polarization angles with respect to the true z-axis
      cpol=zpol
      spol=sqrt(1.d0-zpol**2)
      cpk=xpol/spol
      spk=ypol/spol
c number of points, ipts for theta, ipts1 for phi
      ipts=31
c       ipts=61
      ipts1=61
c      ipts1=121
c Starting loop on phi of photo-e => perad, ppe(phi); (0~2Pi) ###################
c The molecular axis (which here is assumed to be the x-axis) and the polarization 
c vector define a plane; when outputing the grid, the photoelectron azimuthal angle 
c is defined with respect to that plane.
c
      do 888 iphi=1,ipts1
         pephi=(iphi-1)*6+.0001d0
c         pephi=float((iphi-1))*3+.001d0
         if(iphi.eq.ipts1)pephi=359.9999
         ppe(iphi)= pephi*pi/180.
c Take phirad into account in defining the photoelectron coordinates in the body-frame.
         perad=pephi*pi/180. + phirad
         if(perad.gt.twopi)perad=perad-twopi
         cp=cos(perad)
         sp=sin(perad)
c Starting loop on theta of photo-e => theta, th(ipts); (0~Pi) ##################
         do 1000 ithe=1,ipts
            theta=(ithe-1)*6+.0001d0
c            theta=float((ithe-1)*3)+.001d0
            theta = theta*pi/180.
            th(ithe)=theta
            ss=sin(theta)
            cc=cos(theta)
c Get the photoelectron cartesian coordinates
            xel=cc
            zel=ss*cp            
            yel=ss*sp
c now get the angles with respect to the true z-axis
            c=zel
            thtrue=acos(zel)
            s=sqrt(1.d0-zel**2)
            cpkp=xel/s
            spkp=yel/s
            phtrue=acos(cpkp)
c            write(6,157)th(ithe),ppe(iphi),c,s,cpkp,
c     $           spkp,thtrue,phtrue
c 157        format(8f10.5)            
            do 70 izero=1,nenemax
 70         xsec(izero) = 0.0
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
 10               continue
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
 12               continue
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
 49         f(izero) = (0.0,0.0)
            do 25 isym=1,nsym
            if(ithe.eq.1.and.iphi.eq.1)then
               write(6,"(a8)")name(isym)
            endif
            open(8,file=name(isym),form='unformatted',status='unknown')
c
            rewind(8)
            do 24 iene=1,nenergy
c
 66            continue
               read(8)icw,jcw,ni,nj,aki,akj
               twoe(iene)=akj*akj
               read(8)((tmat(i,j),i=1,ni),j=1,nj)
               if(ichan.eq.icw.and.jchan.eq.jcw)go to 67
               go to 66
 67            continue
 222           format(4x,6e12.5)
               degi=sgni(isym)
               do 304 i=1,ni
                  do 304 j=1,nj
 304           tmat(i,j)=tmat(i,j)*degi
c
c calculate total cross sections on first pass through these loops
c
               if(ithe.eq.1.and.iphi.eq.1)then
c
c  print Tmatrix if requested
                  if(iprt.gt.0) then
                     write(6,*)' Tmatrix ',icw,jcw
                     write(6,88001)((tmat(i,j),i=1,ni),j=1,nj)
88001                format(3(4x,f12.8,2x,f12.8))
                  end if 
c  write total cross section contributions for each polarization direction
c  if iprt.gt.0
                  if(iprt.gt.0) then
                     do i=1,ni
                        sum=0.d0
                        do j=1,nj
                           sum = sum + abs(tmat(i,j))**2
                        enddo
                        write(6,
     $                '("sum of |Tij|^2 for i =",i3, " = ",e12.5)')i,sum
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
                        totc(isym,iene)=totc(isym,iene)
     $                       +abs(tmat(ii,jj))**2
 71               continue
c
                  if(iprt.gt.0) then
                     write(6,*)' isym iene totc '
     $                    ,isym,iene,totc(isym,iene)
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
               do 31 i=ifirst,ilast
                  ii=i+1-ifirst
                  li = lch(isym,i)
                  mi = mch(isym,i)
c  Correct inconsistency in definition of phase of Ylms when being 
c  used to specify polarization direction -- CWM and CST 2011
c  Cartesian x,y,z were used in cphot, but real Ylms are used here, hence phase difference
c*tnr* this was done in the code cdipc.f, so it's being commented out here
                  phase = 1.d0
c                     if(mi.eq.1.or.mi.eq.-1) phase = -1.d0
                  do 32 j=jfirst,jlast
                     jj=j+1-jfirst
                     lj = lch(isym,j)
                     mj = mch(isym,j)
                     f(iene)=ylm(li,mi)*phase*tmat(ii,jj)*ylmp(lj,mj)
     $                    +f(iene)
c
c  end loop on j index (electronic angular momentum lm pairs)
 32               continue
c  end loop on i index (dipole operator reduced tensor components)
 31            continue
c end loop on energies 
 24         continue
            close(8)
c end loop on molecular symmetries
 25         continue
            do 50 iene=1,nenergy
c
c accumulate cross section(mfpad) 
c
               xsec(iene)=abs(f(iene))**2
 50         continue
            write(6,157) th(ithe), ppe(iphi)
 157        format(' scattering angle (theta, phi)  =',2 f15.8)
            do 80 iene=1,nenergy
c factors of 2*pi here SEEM to have to do with normalization of phi part of the real Ylm's
c as defined above. -- CWM
               xsec(iene)=xsec(iene)/(2.*pi)
               einc = twoe(iene)/2.
               write(6,101) einc,xsec(iene)
 101           format(' incident e =',e12.5,'  MFPAD =',e12.5)
               egy(iene)=einc
               xs(iphi,ithe,iene)=xsec(iene)
 80         continue
c ending loops on photoelectron theta and phi ###########
 1000    enddo
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
      stop
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
