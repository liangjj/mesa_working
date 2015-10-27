*deck @(#)vwnf.f	5.3  4/17/95
      subroutine vwnf(ngrid,mxgrd,dengrida,dengridb,calc,derivs,
     $               fout,rho,x,zeta,g,h,t1,t2,t3,t4,t5,t6,t7,t8,
     $               t9,t10,t11)
c***begin prologue     vwnf.f
c***date written       930604  
c***revision date      4/17/95      
c   march 4, 1995      rlm at lanl
c      fixing a bug in open-shell implementation.  
c
c***keywords           vosko, wilk, nusair, electron-gas, dft, correlation 
c***author             martin, richard(lanl) 
c***source             @(#)vwnf.f	5.3   4/17/95
c***purpose            evaluates the vosko,wilk,nusair parameterization of 
c                      the uniform electron gas correlation functional from
c                      ceperly and alder. also does derivatives of the
c                      functional if requested. 
c***description
c   computes the vosko,wilk,nusair correlation functional(Version V) and
c   (optionally) its derivatives.
c   
c   derivs             >0 do both functional and derivative
c                      =0 do only functional
c                      <0 do only derivatives
c     
c   fout(1,1)          functional
c   fout(1,2)          derivative wrt rho(alpha) 
c   fout(1,3)
c   fout(1,4)          derivative wrt rho(beta)
c   fout(1,5)       
c
c***references         
c   S.H. Vosko, L. Wilk, and M.Nusair, Can. J. Phys. 58,1200(1980). Version V.
c   D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45,566(1980).
c   U. v. Barth and L. Hedin, J. Phys. C 5,1629 (1972)
c   O. Gunnarsson and B. I. Lundqvist, Phys. Rev. B. 13,4274, (1976)
c
c***routines called
c
c***end prologue       vwnf.f
      implicit none
c     --- input variables -----
      integer ngrid,derivs,mxgrd
      character*(*) calc
c     --- input arrays (unmodified) ---
      real*8 dengrida(ngrid),dengridb(ngrid)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 fout(mxgrd,*)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 t1(ngrid),t2(ngrid),x(ngrid),zeta(ngrid),rho(ngrid)
      real*8 g(ngrid),h(ngrid),t3(ngrid),t4(ngrid)
      real*8 t5(ngrid),t6(ngrid),t7(ngrid),t8(ngrid)
      real*8 t9(ngrid),t10(ngrid),t11(ngrid)
c     --- local variables ---
      integer i,inp,iout
      real*8 ap,bp,cp,x0p
      real*8 af,bf,cf,x0f
      real*8 aa,ba,ca,x0a
      real*8 p,s
      real*8 one,two,three,four,nine,pi,sixth
      real*8 onethrd,frthrd,vwncut
      data one/1.0d+00/,two/2.0d+00/,three/3.0d+00/,four/4.0d+00/
      data nine/9.0d+00/,vwncut/1.0d-20/
      parameter (onethrd=1.0d+00/3.0d+00,frthrd=4.0d+00/3.0d+00)
c
c     see vosko,wilk,nusair; pp.1207,1208,1209.  note that A has been
c     divided by 2 in order to convert it to atomic units. these
c     are paramagnetic, ferromagnetic, and rpa(Johnson et al. call it 
c     'antiferromagnetic') values of the constants.
      data ap,bp,cp,x0p/3.10907d-02,3.72744d+00,1.29352d+01,
     $                 -1.04980d-01/
      data af,bf,cf,x0f/1.554535d-02,7.06042d+00,1.80578d+01,
     $                 -3.25000d-01/         
      data aa,ba,ca,x0a/-1.688686d-02,1.13107d+00,1.30045d+01,
     $                 -4.75840d-03/
c
      save one,two,three,four,nine,vwncut
      save ap,bp,cp,x0p
      save af,bf,cf,x0f
      save aa,ba,ca,x0a
c
      common/io/inp,iout
c
c
      pi=four*atan(one)
c     assign the A value for the rpa component full machine precision.
      aa=-one/(6.0d+00*pi*pi)
      sixth=one/6.0d+00
c
c     --- form the total density, and the spin-independent variable x
      call vadd(rho,dengrida,dengridb,ngrid)
      p=three/(four*pi)
      do 10 i=1,ngrid
         x(i)=(p/rho(i))**sixth
   10 continue
c
      if(calc.eq.'closed') then
         if(derivs.ne.0) then
            call epsc(ap,bp,cp,x0p,x,t3,ngrid,t1,t2)
            call depsc(ap,bp,cp,x0p,x,t4,ngrid,t1,t2)
            call smul(t1,x,sixth,ngrid)
            call vneg(t1,t1,ngrid)
            call vmul(t4,t4,t1,ngrid)
            call vadd(t1,t3,t4,ngrid)
            call vadd(fout(1,2),fout(1,2),t1,ngrid)
            call vadd(fout(1,4),fout(1,4),t1,ngrid)
            if(derivs.gt.0) then
c              --- do the functional ---
               call vmul(t3,rho,t3,ngrid)
               call vadd(fout(1,1),fout(1,1),t3,ngrid)
            endif
         else if (derivs.eq.0) then
c           --- functional only ---
            call epsc(ap,bp,cp,x0p,x,t3,ngrid,t1,t2)
            call vmul(t3,rho,t3,ngrid)
            call vadd(fout(1,1),fout(1,1),t3,ngrid)
         endif 
      else
c        --- open shell ---
c            form the spin-dependent variable zeta
         call vsub(zeta,dengrida,dengridb,ngrid)
         call vdiv(zeta,zeta,rho,ngrid) 
c        --- form g(zeta) ---
         s=9.0d+00/8.0d+00
         do 100 i=1,ngrid
            g(i)=s*((one+zeta(i))**frthrd + (one-zeta(i))**frthrd -2)
  100    continue
c        --- need ec(x,zeta) for both functional and derivative ---
         call epsc(ap,bp,cp,x0p,x,t5,ngrid,t1,t2)
c        --- form h ---
         call epsc(af,bf,cf,x0f,x,t3,ngrid,t1,t2)
         call vsub(t3,t3,t5,ngrid)
         call epsc(aa,ba,ca,x0a,x,t4,ngrid,t1,t2)
         call vdiv(t3,t3,t4,ngrid)
         s=four/(nine*((two**onethrd)-1))
         call smul(t3,t3,s,ngrid)
         call ssub(h,t3,one,ngrid)
c        --- form zeta**4 and zeta**3---
         call vmul(t3,zeta,zeta,ngrid)
         call vmul(t10,t3,zeta,ngrid)
         call vmul(t3,t3,t3,ngrid)
c        --- at this point, h and g are formed, 
c            t4 contains ec(a) and t3 contains zeta**4
c            now form (hzeta**4+1)*g*e_c^a into t1
         call vmul(t1,h,t3,ngrid)
         call sadd(t1,t1,one,ngrid)
         call vmul(t1,g,t1,ngrid)
         call vmul(t1,t4,t1,ngrid)
c        --- accumulate ec(x,zeta) =e_c^p + (1+h*zeta**4)*g*e_c^A in t9 ---
         call vadd(t9,t1,t5,ngrid)
c
         if(derivs.ne.0) then
c           --- at this point, h and g are formed, 
c               t4 contains ec(a) and t3 contains zeta**4
c               form g' 
            s=three/two
            do 120 i=1,ngrid
               t1(i)=s*((one+zeta(i))**onethrd - (one-zeta(i))**onethrd)
  120       continue
c           --- (1+h*zeta**4)*g' into t2
            call vmul(t2,h,t3,ngrid) 
            call sadd(t2,t2,one,ngrid)
            call vmul(t2,t2,t1,ngrid)
c           --- 4*g*h*zeta**3 into t5
            call vmul(t5,g,h,ngrid)
            call smul(t5,t5,four,ngrid)
            call vmul(t5,t5,t10,ngrid)
c           Ok.  One more time, this is your brain, this is VWN, 
c           this is your brain on VWN. We now have:
c           t1=g', t2=(h*zeta**4+1)*g', t3=zeta**4, t4=e_c^A, t5=g*h*4*zeta**3
c           t6=empty, t7=empty, t8=empty, t9=e_c(x,zeta), t10=zeta**3, t11=empty
c           --- sum all the d/dzeta terms, multiply by ec(a), and store in t2
            call vadd(t2,t5,t2,ngrid)
            call vmul(t2,t2,t4,ngrid)
c           --- so t2 has e_c^A(g'*(1+h*zeta**4)+4*g*h*zeta**3), Ok?
c           --- form h'into t5---
            call epsc(ap,bp,cp,x0p,x,t5,ngrid,t10,t11)
            call epsc(af,bf,cf,x0f,x,t6,ngrid,t10,t11)
            call vsub(t5,t5,t6,ngrid)  
            call vdiv(t5,t5,t4,ngrid)
            call depsc(aa,ba,ca,x0a,x,t7,ngrid,t10,t11)
            call vmul(t5,t5,t7,ngrid)
            call depsc(ap,bp,cp,x0p,x,t6,ngrid,t10,t11)
            call vsub(t5,t5,t6,ngrid)
            call depsc(af,bf,cf,x0f,x,t8,ngrid,t10,t11)
            call vadd(t5,t5,t8,ngrid)
            s=four/(nine*((two**onethrd)-1))
            call smul(t5,t5,s,ngrid)
            call vdiv(t5,t5,t4,ngrid)
c           --- now t6 has e_c^P', t7 has e_c^A', t8 has e_c^F, t5 has h'
c           --- ec(a)*g*h'*zeta**4 into t5
            call vmul(t5,t5,t3,ngrid)
            call vmul(t5,g,t5,ngrid)
            call vmul(t5,t5,t4,ngrid)
c           --- ec(a)'*g*(1+h*zeta**4) into t1
            call vmul(t1,h,t3,ngrid) 
            call sadd(t1,t1,one,ngrid)
            call vmul(t1,t1,g,ngrid)
            call vmul(t1,t1,t7,ngrid)
c           --- add in ec(p)'
            call vadd(t1,t1,t6,ngrid)
c           --- add in the ec(a) term 
            call vadd(t1,t1,t5,ngrid)
c           --- now times -x/6*rho (only don't do /rho, coz we do a mul later
            call vmul(t1,t1,x,ngrid)
            call smul(t1,t1,-sixth,ngrid)
c
c           --- finally, we are nearing home ---
c               derivative with respect to alpha spin.
c           ---put d(zeta)/drhoa into t5
            call ssub(t5,zeta,one,ngrid)
            call smul(t5,t5,-one,ngrid)
c           --- form the d(zeta)/drhoa term
            call vmul(t5,t2,t5,ngrid)
c           --- add in the stuff that doesn't have d(zeta)/drhoa
            call vadd(t5,t5,t1,ngrid)
c           --- add in the functional to form the derivative in t10
            call vadd(t10,t9,t5,ngrid)
            call vadd(fout(1,2),fout(1,2),t10,ngrid)
c
c           --- now the derivative with respect to beta.
c           ---put d(zeta)/drhob into t5
            call sadd(t5,zeta,one,ngrid)
            call smul(t5,t5,-one,ngrid)
c           --- form the d(zeta)/drhob term
            call vmul(t5,t2,t5,ngrid)
c           --- add in the stuff that doesn't have d(zeta)/drhob
            call vadd(t5,t5,t1,ngrid)
            call vadd(t11,t9,t5,ngrid)
            call vadd(fout(1,4),fout(1,4),t11,ngrid)
            if(derivs.gt.0) then
c              --- want functional as well
               call vmul(t9,t9,rho,ngrid)
               call vadd(fout(1,1),fout(1,1),t9,ngrid)
            endif
         else 
c           functional only.
            call vmul(t9,t9,rho,ngrid)
            call vadd(fout(1,1),fout(1,1),t9,ngrid)
         endif
      endif
c
c
      return
      end
      subroutine epsc(a,b,c,x0,x,ec,ngrid,t1,capxi)
c***begin prologue     vwnf.f
c***date written       930604  
c***revision date      4/17/95      
c
c***keywords           dft,electron-gas,vwn,correlation
c***author             martin, richard (lanl) 
c***source             @(#)vwnf.f	5.3   4/17/95
c***purpose            
c   evaluates the spin-independent portion of the correlation energy
c   of vosko, wilk, and nusair, which is derived from the electron gas
c   calculations of ceperly and alder.   
c***description
c   a,b,c,x0           input constants which depend on the portion of
c                      the correlation potential being evaluated,
c                      i.e., para-magnetic or ferro-magnetic. 
c   x                  (3/(4*pi*rho))**1/6 provided on a grid of
c                      ngrid points on input.
c   ec                 the correlation function, returned on output.
c   ngrid              the number of grid points.
c   t1,capxi           scratch arrays(ngrid) 
c***references
c   B.G. Johnson, P.M.W. Gill and J.A. Pople, J.Chem.Phys. 98,5612(1993).
c   S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58,1200(1980)
c
c***routines called
c
c***end prologue       vwnf.f
      implicit none
c     --- input variables -----
      integer ngrid
      real*8 a,b,c,x0
c     --- input arrays (unmodified) ---
      real*8 x(ngrid)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 ec(ngrid)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 t1(ngrid),capxi(ngrid)
c     --- local variables ---
      integer i
      integer inp,iout
      real*8 p,q,r,s
      real*8 one,two,four
      real*8 atan
c
      parameter (one=1.0d+00,two=2.0d+00,four=4.0d+00)
      common/io/inp,iout
c
      q=sqrt(four*c-b*b)
      p=two*(two*x0+b)/q
      do 10 i=1,ngrid
         t1(i)=atan(q/(two*x(i)+b))
         capxi(i)=one/(x(i)*x(i)+b*x(i)+c)
         ec(i)=p*t1(i)
   10 continue
c
c
      r=-b*x0/(x0*x0+b*x0+c)
      do 20 i=1,ngrid
         ec(i)=r*(log((x(i)-x0)*(x(i)-x0)*capxi(i))+ec(i))
   20 continue
c
      s=two*b/q 
      do 30 i=1,ngrid
         ec(i)=a*(ec(i)+log(x(i)*x(i)*capxi(i))+s*t1(i)) 
   30 continue
c
c
      return
      end
      subroutine depsc(a,b,c,x0,x,dec,ngrid,t1,t2,t3,t4)
c***begin prologue     vwnf.f
c***date written       930608  
c***revision date      4/17/95      
c
c***keywords           dft,electron-gas,vwn,correlation,derivatives
c***author             martin, richard (lanl) 
c***source             @(#)vwnf.f	5.3   4/17/95
c***purpose            
c   evaluates the derivative with respect to x of the
c   spin-independent portion of the correlation energy of 
c   vosko,wilk and nusair.
c***description
c   a,b,c,x0           input constants which depend on the portion of
c                      the correlation function being evaluated,
c                      i.e., para-magnetic or ferro-magnetic. 
c   x                  (3/(4*pi*rho))**1/6 provided on a grid of
c                      ngrid points on input.
c   dec                the derivative of the correlation function, 
c                      returned on output.
c   ngrid              the number of grid points.
c   t1,t2,t3,t4        scratch arrays(ngrid) 
c***references
c   B.G. Johnson, P.M.W. Gill and J.A. Pople, J.Chem.Phys. 98,5612(1993).
c   S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58,1200(1980)
c
c***routines called
c
c***end prologue       vwnf.f
      implicit none
c     --- input variables -----
      integer ngrid
      real*8 a,b,c,x0
c     --- input arrays (unmodified) ---
      real*8 x(ngrid)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 dec(ngrid)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 t1(ngrid),t2(ngrid)
      real*8 t3(ngrid),t4(ngrid)
c     --- local variables ---
      integer i
      integer inp,iout
      real*8 p,qsq,r
      real*8 one,two,four
c
      parameter (one=1.0d+00,two=2.0d+00,four=4.0d+00)
      common/io/inp,iout
c
      call rzero(dec,ngrid)
      qsq=four*c-b*b
      do 10 i=1,ngrid
         t1(i)=one/((two*x(i)+b)*(two*x(i)+b)+qsq)
         t2(i)=-(two*x(i)+b)/(x(i)*x(i)+b*x(i)+c)
   10 continue
c
      p=-four*(two*x0+b)
      r=-b*x0/(x0*x0+b*x0+c)
      call vwxs(dec,t2,t1,p,+1,ngrid)
      do 20 i=1,ngrid
         dec(i)=r*(dec(i)+two/(x(i)-x0))
   20 continue
c
      r=-four*b
      call vwxs(dec,dec,t1,r,+1,ngrid)
      call vadd(dec,dec,t2,ngrid)
      do 30 i=1,ngrid
         dec(i)=a*(dec(i)+two/x(i))
   30 continue
c
c
      return
      end
