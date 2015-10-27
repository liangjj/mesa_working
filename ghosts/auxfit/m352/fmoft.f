      subroutine fmoft(f,mmax,t)
c***begin prologue     fmoft.f
c***date written       920427  (yymmdd) 
c***revision date      yymmdd  (yymmdd)
c
c***keywords           
c***author             martin, richard(lanl)
c***source             @(#)fmoft.f	2.1   4/29/92
c***purpose            
c                      form the integral f(m,t) for m=0,mmax
c                            1
c                        int   [(u**2m) exp(-tu**2)  du]
c                            0

c***description
c
c        fmoft computes f (t)  for m=0,mmax
c                        m
c
c
c***references
c      mcmurchie and davidson, j.comp.phys., 26, 218(1978).
c      this subroutine was adapted from MELD(McMurchie,Elbert,Langhoff,Davidson)
c***routines called
c
c***end prologue       fmoft.f
c
      implicit integer(a-z)
c
      real*8 t
c
c     output arrays
      real*8 f(0:mmax)
c
      real*8 cc,fq,fqi
      real*8 zero,f12,f15,f18,f24,f30,f36,f2mp36
      real*8 fiveh,tenth,twoten,sixth,quartr,third,half,one
      real*8 two,ten
      real*8 const1,const2,const3,const4,const5,const6,const7
      real*8 const8,const9,cons10,cons11
      parameter (zero=0.0d0,f12=12.0d0,f15=15.0d0,f18=18.0d0)
      parameter (f24=24.0d0,f30=30.0d0,f36=36.0d0)
      real*8 d,y,ti,twoti
c
      common /over/cc(2760),fq(0:16),fqi(0:16)
      data const1 /0.88622692545276d0/,const2 /-0.38115593460d0/,
     *     const3 /0.3211809090d0/,    const4 /0.24736316860d0/,
     *     const5 /0.4999489092d0/,    const6 /0.2464284500d0/,
     *     const7 /0.2424943800d0/,    const8 /0.4998436875d0/,
     *     const9 /0.4990931620d0/,    cons10 /0.2152832d0/,
     *     cons11 /0.490d0/
      data fiveh/0.05d0/,tenth/0.1d0/,twoten/0.2d0/,
     *     sixth/0.16666666666666d0/,quartr/0.25d0/,
     *     third/0.33333333333333d0/,half/0.5d0/,one/1.0d0/,
     *     two/2.0d0/,ten/10.0d0/
      save const1,const2,const3,const4,const5,const6,const7,const8
      save const9,cons10,cons11
      save fiveh,tenth,twoten,sixth,quartr,third,half,one,two,ten
c
c
      f2mp36=float(2*mmax)+f36
      if(t.lt.zero) then
         call lnkerr('bad argument --negative t-- to fmoft')
      else if(t.eq.zero) then
c
c     ----- t=0 -----
         do 10 m=0,mmax
            f(m)=fqi(m)
   10    continue
      else if(t.gt.zero.and.t.lt.f12) then
c
c     ----- 0<t<12 -----
c        use series expansion in (t*-t) to get f(mmax)
c        find t*-t
         j=ifix(t*ten)
         d=tenth*float(j)+fiveh-t
c        find pointer into series expansion
         k=23*j+17-mmax
c        evaluate seven term taylor expansion
         f(mmax)=(d*(half*d*(third*d*(quartr*d*(twoten*d*(
     *              sixth*d*cc(k)+cc(k+1))+cc(k+2))+cc(k+3))
     *              +cc(k+4))+cc(k+5))+cc(k+6))
c
c        downward recursion for f(m,t)   j=mmax-1,0
         y=exp(-t)
         do 20 m=mmax-1,0,-1
            f(m)=(two*t*f(m+1)+y)*fqi(m)
   20    continue
      else if(t.ge.f12.and.t.lt.f30) then
c
c     ----- 12 <= t < 30 -----
c        form f(0,t)
         if(t.lt.f15) then
            ti=one/t
            f(0)=((const2*ti + const3)*ti - const4)*ti + const5
         else if(t.lt.f18) then
            ti=one/t
            f(0)=(const6*ti - const7)*ti + const8
         else if(t.lt.f24) then
            f(0)=const9 -  cons10/t
         else
            f(0)=cons11
         endif
         y=exp(-t)
         twoti=one/(two*t)
         f(0)=(const1/sqrt(t) -y*f(0)/t)
c        generate remaining terms with upward recursion
         do 30 m=0,mmax-1
            f(m+1)= twoti*(fq(m)*f(m) -y)
   30    continue
      else if(t.ge.f30.and.t.lt.f2mp36) then
c
c     ----- 30 <= t < 2*mmax+36
c        form f(0)
         y=exp(-t)
         twoti=one/(two*t)
         f(0)=const1/sqrt(t)        
c        generate remaining terms with exact upward recursion
         do 40 m=0,mmax-1
            f(m+1)=twoti*(fq(m)*f(m) -y)
   40    continue 
      else if(t.ge.f2mp36) then
c
c     ----- t >= 2*mmax+36
c        form f(0)
         twoti=one/(two*t)
         f(0)=const1/sqrt(t)
c        generate remaining terms with approximate upward recursion.
         do 50 m=0,mmax-1
            f(m+1)=twoti*fq(m)*f(m)
   50    continue
      endif
c
c
      return
      end
