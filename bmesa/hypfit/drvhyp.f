*deck drvhyp.f
c***begin prologue     drvhyp
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for hyperspherical function.
c***                   
c***references         
c
c***routines called    
c***end prologue       drvhyp
      subroutine drvhyp(pgrid,energy,mmax,n,nen,type,dertyp,
     1                  fitprd,prn)
      implicit integer (a-z)
      integer*8 pgrid(2)
      real*8 grid1, grid2, energy, fmk, scr, k
      logical prn, fitprd
      character*(*) type, dertyp
      character*80 title
      character*2 mstr, itoc
      character*24 str
      dimension energy(nen), n(2), ngot(2), prn(*)
      common/io/inp, iout
      pointer (pg1,grid1(1))
      pointer (pg2,grid2(1))
      pointer (pmk,fmk(1))
      pointer (pscr,scr(1))
c
c     define locations of points, weights and functions.
c
      pg1=pgrid(1)
      pg2=pgrid(2)
      r1=1
      wtr1=r1+n(1)
      fr1=wtr1+n(1)
      dfr1=fr1+n(1)*n(1)
      r2=1
      wtr2=r2+n(2)
      fr2=wtr2+n(2)
      dfr2=fr2+n(2)*n(2)
c
c     bessel parameters for recursion and memory allocation
c
      prd=n(1)*n(2)
      rho=1
      jf=rho+prd
      nf=jf+prd
      djf=nf+prd
      dnf=djf+prd
      ang=dnf+prd
      y=ang+prd
      dy=y+prd
      c=dy+prd
      dfdr1=c+prd
      dfdr2=dfdr1+prd
      df=dfdr2+prd
      df1=df+prd
      need=df1+prd
      if(fitprd) then
         ys=need
         dys1=ys+prd
         dys2=dys1+prd
         need=dys2+prd
      endif                  
      need=wpadti(need)
      call memory(need,pmk,ngot(1),'main',0)
c
c     get some scratch space for the angular and bessel functions
c     at one point and the factorial function.
c
      lmax = m + m +2
      jb=1
      djb=jb+lmax+1
      yb=djb+lmax+1
      dyb=yb+lmax+1
      fact=dyb+lmax+1
      need=wpadti(fact+101)
      call memory(need,pscr,ngot(2),'scr',0) 
      write(iout,1) n(1), n(2), mmax, lmax
c
c     calculate the arrays giving the radial and angular variables
c
      call factl(scr(fact),100)
      call radang(fmk(rho),fmk(ang),grid1(r1),grid2(r2),n(1),n(2))
      if(prn(1)) then
         title='rho(r2,r1)'
         call prntfm(title,fmk(rho),n(2),n(1),n(2),n(1),iout)
         title='alpha(r2,r1)'
         call prntfm(title,fmk(ang),n(2),n(1),n(2),n(1),iout)
      endif 
c
      do 10 ene=1,nen
         k=sqrt(2.d0*energy(ene))
         do 20 m=0,mmax
            mstr=itoc(m)
            l1=length(mstr)
c
c           compute the functions and their derivatives wrt to both
c           angular, radial and cartesian variables.
c
            call hypang(fmk(ang),fmk(y),fmk(dy),m,n(1),n(2))
            call hyprad(fmk(rho),fmk(jf),fmk(djf),
     1                  fmk(nf),fmk(dnf),k,scr(jb),scr(djb),
     3                  scr(yb),scr(dyb),scr(fact),m,n(1),n(2),lmax,
     4                  type)
            if(prn(2)) then
               str='(alpha:'//mstr(1:l1)//')'
               l2=length(str)
               title='Exact Angular Function'//str(1:l2)
               l3=length(title)
               call prnfn(fmk(y),fmk(ang),title(1:l3),'Angle',
     1                    n(1),n(2))
               title='Exact Derivative of Angular Function'//str(1:l2)//
     1               ' wrt angle'
               l3=length(title)
               call prnfn(fmk(dy),fmk(ang),title(1:l3),
     1                    'Angle',n(1),n(2))
               str='(kr:'//mstr(1:l1)//')'
               l2=length(str)
               title='Regular Function:'//str(1:l2)
               l3=length(title)
               call prnfn(fmk(jf),fmk(rho),title(1:l3),
     1                    'Radius',n(1),n(2))
               title='Derivative of Regular Function'//str(1:l2)//
     1               ' wrt r'
               l3=length(title)
               call prnfn(fmk(djf),fmk(rho),title(1:l3),
     1                    'Radius',n(1),n(2))
               title='Irregular Function'//str(1:l2)
               l3=length(title)
               call prnfn(fmk(nf),fmk(rho),title(1:l3),
     1                    'Radius',n(1),n(2))
               title='Derivative of Irregular Function'//str(1:l2)
     1               //' wrt r'
               l3=length(title)
               call prnfn(fmk(dnf),fmk(rho),title(1:l3),
     1                    'Radius',n(1),n(2))
            endif
            if(fitprd) then
               str='General'
c
c              make product hyperspherical functions and derivatives
c              from components.
c
               call prodfn(fmk(ys),fmk(dys1),fmk(dys2),
     1                     fmk(y),fmk(dy),fmk(jf),fmk(djf),
     2                     n(1),n(2))
c
c              get dvr coefficients to fit function
c
               call coefs(fmk(c),fmk(ys),grid1(wtr1),grid1(fr1),
     1                    grid2(wtr2),grid2(fr2),n(1),n(2))
               if(dertyp.eq.'cartesian') then
c
c                 get cartesian derivatives of exact and approximate
c                 functions and compare.
c
                  call cartder(fmk(dfdr1),fmk(dfdr2),
     1                         fmk(dys1),fmk(dys2),
     2                         fmk(df),fmk(df1),
     3                         fmk(rho),fmk(ang),fmk(c),
     4                         grid1(r1),grid2(r2),
     5                         grid1(fr1),grid2(fr2),
     6                         grid1(dfr1),grid2(dfr2),
     7                         str,n(1),n(2))
               elseif(dertyp.eq.'hyperspherical') then
                  call hypder(fmk(dfdr1),fmk(dfdr2),
     1                        fmk(dys1),fmk(dys2),
     2                        fmk(df),fmk(df1),
     2                        fmk(rho),fmk(ang),fmk(c),
     3                        grid1(fr1),grid2(fr2),
     4                        grid1(dfr1),grid2(dfr2),
     5                        str,n(1),n(2))
               else
                  call lnkerr('error in derivative type')
               endif
               call prodfn(fmk(ys),fmk(dys1),fmk(dys2),
     1                     fmk(y),fmk(dy),fmk(nf),fmk(dnf),
     2                     n(1),n(2))
c
c              get dvr coefficients to fit function
c
               call coefs(fmk(c),fmk(ys),grid1(wtr1),grid1(fr1),
     1                    grid2(wtr2),grid2(fr2),n(1),n(2))
               if(dertyp.eq.'cartesian') then
c
c                 get cartesian derivatives of exact and approximate
c                 functions and compare.
c
                  call cartder(fmk(dfdr1),fmk(dfdr2),
     1                         fmk(dys1),fmk(dys2),
     2                         fmk(df),fmk(df1),
     3                         fmk(rho),fmk(ang),fmk(c),
     4                         grid1(r1),grid2(r2),
     5                         grid1(fr1),grid2(fr2),
     6                         grid1(dfr1),grid2(dfr2),
     7                         str,n(1),n(2))
               elseif(dertyp.eq.'hyperspherical') then
                  call hypder(fmk(dfdr1),fmk(dfdr2),
     1                        fmk(dys1),fmk(dys2),
     2                        fmk(df),fmk(df1),
     2                        fmk(rho),fmk(ang),fmk(c),
     3                        grid1(fr1),grid2(fr2),
     4                        grid1(dfr1),grid2(dfr2),
     5                        str,n(1),n(2))
               else
                  call lnkerr('error in derivative type')
               endif
            else
               call coefs(fmk(c),fmk(y),grid1(wtr1),grid1(fr1),
     1                    grid2(wtr2),grid2(fr2),n(1),n(2))
               str='Angular'
               if(dertyp.eq.'cartesian') then
c
c                 get cartesian derivatives of exact and approximate
c                 functions and compare.
c

                  call cartder(fmk(dfdr1),fmk(dfdr2),
     1                         fmk(dy),fmk(dy),
     2                         fmk(df),fmk(df1),
     3                         fmk(rho),fmk(ang),fmk(c),
     4                         grid1(r1),grid2(r2),
     5                         grid1(fr1),grid2(fr2),
     6                         grid1(dfr1),grid2(dfr2),
     7                         str,n(1),n(2))
               elseif(dertyp.eq.'hyperspherical') then
                  call hypder(fmk(dfdr1),fmk(dfdr2),
     1                        fmk(dy),fmk(dy),
     2                        fmk(df),fmk(df1),
     2                        fmk(rho),fmk(ang),fmk(c),
     3                        grid1(fr1),grid2(fr2),
     4                        grid1(dfr1),grid2(dfr2),
     5                        str,n(1),n(2))
               else
                  call lnkerr('error in derivative type')
               endif
               call coefs(fmk(c),fmk(jf),grid1(wtr1),grid1(fr1),
     1                    grid2(wtr2),grid2(fr2),n(1),n(2))
               str='Radial'
               if(dertyp.eq.'cartesian') then
c
c                 get cartesian derivatives of exact and approximate
c                 functions and compare.
c
                  call cartder(fmk(dfdr1),fmk(dfdr2),
     1                         fmk(djf),fmk(djf),
     2                         fmk(df),fmk(df1),
     3                         fmk(rho),fmk(ang),fmk(c),
     4                         grid1(r1),grid2(r2),
     5                         grid1(fr1),grid2(fr2),
     6                         grid1(dfr1),grid2(dfr2),
     7                         str,n(1),n(2))
               elseif(dertyp.eq.'hyperspherical') then
                  call hypder(fmk(dfdr1),fmk(dfdr2),
     1                        fmk(djf),fmk(djf),
     2                        fmk(df),fmk(df1),
     2                        fmk(rho),fmk(ang),fmk(c),
     3                        grid1(fr1),grid2(fr2),
     4                        grid1(dfr1),grid2(dfr2),
     5                        str,n(1),n(2))
               else
                  call lnkerr('error in derivative type')
               endif
               call coefs(fmk(c),fmk(nf),grid1(wtr1),grid1(fr1),
     1                    grid2(wtr2),grid2(fr2),n(1),n(2))
               str='Radial'
               if(dertyp.eq.'cartesian') then
c
c                 get cartesian derivatives of exact and approximate
c                 functions and compare.
c
                  call cartder(fmk(dfdr1),fmk(dfdr2),
     1                         fmk(dnf),fmk(dnf),
     2                         fmk(df),fmk(df1),
     3                         fmk(rho),fmk(ang),fmk(c),
     4                         grid1(r1),grid2(r2),
     5                         grid1(fr1),grid2(fr2),
     6                         grid1(dfr1),grid2(dfr2),
     7                         str,n(1),n(2))
               elseif(dertyp.eq.'hyperspherical') then
                  call hypder(fmk(dfdr1),fmk(dfdr2),
     1                        fmk(dnf),fmk(dnf),
     2                        fmk(df),fmk(df1),
     2                        fmk(rho),fmk(ang),fmk(c),
     3                        grid1(fr1),grid2(fr2),
     4                        grid1(dfr1),grid2(dfr2),
     5                        str,n(1),n(2))
               else
                  call lnkerr('error in derivative type')
               endif
            endif
 20      continue   
 10   continue   
      call memory(-ngot(1),pmk,idum,'hyperfn',idum)
      call memory(-ngot(2),pscr,idum,'scr',idum) 
      return
 1    format(/,5x,'number points in r1              = ',i4,
     1       /,5x,'number points in r2              = ',i4,
     2       /,5x,'upper value of m                 = ',i3,
     3       /,5x,'upper value of l                 = ',i4)
      end       

















