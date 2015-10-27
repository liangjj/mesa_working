*deck abcoef
      subroutine abcoef(a,b,a0,b0,c0,d0,eta,l,n,energy,prnt)
c***begin prologue     abcoef
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            expansion coefficients for coulomb functions
c***                   at small rho.
c***description
c***references         nbs mathematical handbook
c
c***routines called
c
c***end prologue       abcoef
c
      implicit integer (a-z)
      real*8 a, b, eta, twoeta, plfun, plf, zero, one, two, sqrt2
      real*8 gamma, psi, a0, b0, c0, d0, argdum, zpoint
      dimension a(l+1:l+1+n), b(-l:-l+n)
      dimension a0(0:n), b0(0:n), c0(0:n), d0(0:n)
      character*(*) energy
      logical prnt
      data zero, one, two, sqrt2 /0.d0,1.d0,2.d0,
     1                            1.41421356237309505d+00/
      common /io/ inp, iout
      if (energy.eq.'positive') then
          twoeta=two*eta
c**********************************************************************c
c             generate the a coefficients needed for the regular       c
c                            series solution                           c
c**********************************************************************c
          a(l+1)=one
          a(l+2)=eta/(l+1)
          do 10 i=l+3,l+1+n
             a(i)=twoeta*a(i-1)-a(i-2)
             a(i)=a(i)/((i+l)*(i-l-1))
   10     continue
c**********************************************************************c
c              now for the b coefficients needed for the irregular     c
c              solution, things are more complicated since             c
c              they depend on the previous generation of the a's       c
c**********************************************************************c
          plf=plfun(eta,l)
          b(-l)=one
          do 20 i=-l+1,-l+n
             if (i.le.l) then
                 zpoint=zero
             else
                 zpoint=a(i)
             endif
             if (i.eq.l+1) then
                 b(i)=zero
             else
                 b(i)=twoeta*b(i-1)-b(i-2)-(i+i-1)*plf*zpoint
                 b(i)=b(i)/((i-l-1)*(i+l))
             endif
   20     continue     
          if (prnt) then
              write(iout,*)
              write(iout,*) '          the a and b expansion '//
     1                      'coefficients'
              write(iout,*)
              write(iout,*) '      ia          a           ib'//
     1                      '          b'
              do 30 i=0,n
                 ia=l+1+i
                 ib=-l+i
                 write(iout,1) ia,a(ia),ib,b(ib)
   30         continue
          endif
      elseif (energy.eq.'negative') then
c**********************************************************************c
c             coefficients for series solution for negative energy     c
c             regular and irregular coulomb functions at small         c
c                               distances.                             c
c          * the indexing is less fancy here and follows the paper     c
c                     of Henry and Roundtree in CPC                    c
c          * there are errors in that paper which will be noted        c
c          * the coefficients defined here are identical to the        c
c            paper cited and the errors have been fixed in the         c
c            definition of the small rho and large rho behavior        c
c                           of the functions                           c
c**********************************************************************c
          argdum=eta+l+1
          a0(0)=gamma(argdum)
          argdum=l+l+2
          a0(0)=a0(0)/gamma(argdum)
          do 40 i=1,n
             a0(i)=(i+eta+l)*a0(i-1)/(i*(i+l+l+1))
   40     continue
          argdum=1.d0
          d0(0)=psi(argdum)
          argdum=l+l+2
          d0(0)=d0(0)+psi(argdum)
          argdum=eta+l+1
          d0(0)=d0(0)-psi(argdum)              
          do 50 i=1,n
             d0(i)=d0(i-1)+1.d0/i+1.d0/(i+l+l+1.d0)-1.d0/(i+eta+l)
   50     continue
          do 60 i=0,n
             b0(i)=a0(i)*d0(i)
   60     continue
          argdum=l+l+1.d0
          c0(0)=gamma(eta-l)*gamma(argdum)
          if (l.gt.0) then
              do 70 i=1,l+l
                 c0(i)=(i+eta-l-1)*c0(i-1)/(i*(l+l+1-i))
   70         continue
          endif
          if (prnt) then
              write(iout,*)
              write(iout,*) '          the a,b,c and d expansion '//
     1                      'coefficients'
              write(iout,*)
              write(iout,*) '      i          a                b'
              do 80 i=0,n
                 write(iout,2) i,a0(i),b0(i)
   80         continue
              write(iout,*)
              write(iout,*) '      i          c                d'
              twoel=l+l
              do 90 i=0,n
                 if (i.le.twoel) then
                     write(iout,2) i,c0(i),d0(i)
                 else
                     write(iout,3) i,d0(i)
                 endif
   90         continue
          endif
      endif     
    1 format(5x,i3,3x,e15.8,3x,i3,3x,e15.8)   
    2 format(5x,i3,3x,e15.8,3x,e15.8)   
    3 format(5x,i3,21x,e15.8)   
      return
      end















