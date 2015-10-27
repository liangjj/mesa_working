*deck rcbssl
c***begin prologue     rcbssl
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate ricatti-bessel functions
c***description        ricatti-bessel functions calculated using forward or
c***                   backward recursion. if greater accuracy needed
c***                   change parameter statement.
c***references       
c
c***routines called
c***end prologue       rcbssl
      subroutine rcbssl(x,j,jp,y,yp,np,lmax,ltop,der,type,prnt)
      implicit integer (a-z)
      parameter (acc=30)
      common /io/ inp, iout
      real*8 one, two, x, j, jp, y, yp, tmp
      logical prnt
      character*(*) der, type
      dimension x(np), j(np,0:ltop), jp(np,0:ltop)
      dimension y(np,0:ltop), yp(np,0:ltop)
      data one, two/ 1.d+00, 2.d+00 /
c----------------------------------------------------------------------c
c            estimate starting l                                       c
c----------------------------------------------------------------------c
      tmp=lmax*acc
      strtl=lmax+sqrt(tmp)
      strtl=max(strtl,lmax)
      if (strtl.gt.ltop) then
          call lnkerr('starting l bigger than ltop')
      endif
c----------------------------------------------------------------------c
c               make first two bessel functions                        c
c----------------------------------------------------------------------c
      if ( type.eq.'regular'.or.type.eq.'both') then
           do 10 i=1,np
              j(i,0)=sin(x(i))
   10      continue
           if (der.eq.'derivatives') then
               do 5 i=1,np
                  jp(i,0)=cos(x(i))
    5          continue
           endif
      endif   
      if ( type.eq.'irregular'.or.type.eq.'both') then
           do 20 i=1,np
              y(i,0)=-cos(x(i))
   20      continue
           if (der.eq.'derivatives') then   
               do 25 i=1,np
                  yp(i,0)=sin(x(i))
   25          continue
           endif    
      endif   
      if (lmax.ge.1) then
          if ( type.eq.'regular'.or.type.eq.'both') then       
               do 30 i=1,np
                  j(i,1)=j(i,0)/x(i)-cos(x(i))
   30          continue
               if (der.eq.'derivatives') then
                   do 35 i=1,np 
                      jp(i,1)=j(i,0)-j(i,1)/x(i)
   35              continue
               endif                      
          endif
          if ( type.eq.'irregular'.or.type.eq.'both') then
               do 40 i=1,np
                  y(i,1)=y(i,0)/x(i)-sin(x(i))
   40          continue
               if (der.eq.'derivatives') then
                   do 45 i=1,np
                      yp(i,1)=y(i,0)-y(i,1)/x(i) 
   45              continue
               endif
          endif                          
      endif
      if (lmax.ge.2) then
c----------------------------------------------------------------------c
c                   lmax is greater than one                           c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c              calculate y by upward recursion                         c
c                    its always stable                                 c
c----------------------------------------------------------------------c
          if (type.eq.'irregular'.or.type.eq.'both') then
              do 50 l=1,strtl-1
                 ll1=l+l+1
                 lm1=l-1
                 lp1=l+1
                 do 60 i=1,np
                    y(i,lp1)=ll1*y(i,l)/x(i)-y(i,lm1)
   60            continue
   50         continue
              if (der.eq.'derivatives') then
                  do 65 l=1,lmax
                     lm=l-1
                     do 75 i=1,np
                        yp(i,l)=y(i,lm)-y(i,l)/x(i) 
   75                continue
   65             continue
              endif
          endif                 
c----------------------------------------------------------------------c
c             calculate j by upward or downward recursion              c
c             depending on the value of x                              c
c----------------------------------------------------------------------c
          if ( type.eq.'regular'.or.type.eq.'both') then  
               call rzero(j(1,strtl),np)
               onelss=strtl-1
               call vfill(j(1,onelss),1.d-100,np)             
               do 70 i=1,np
                  if (x(i).gt.dfloat(lmax)) then
c----------------------------------------------------------------------c
c               upward recursion                                       c
c----------------------------------------------------------------------c
                      do 80 l=1,strtl-1
                         ll1=l+l+1
                         lm1=l-1
                         lp1=l+1
                         j(i,lp1)=ll1*j(i,l)/x(i)-j(i,lm1)
   80                 continue
                  else
c----------------------------------------------------------------------c
c               downward recursion                                     c
c----------------------------------------------------------------------c
                      do 90 l=onelss-1,0,-1
                         ll3=l+l+3
                         lp1=l+1
                         lp2=l+2
                         j(i,l)=ll3*j(i,lp1)/x(i)-j(i,lp2)
   90                continue
c----------------------------------------------------------------------c
c                normalize the j                                       c
c----------------------------------------------------------------------c
                     do 100 l=strtl,0,-1
                        j(i,l)=j(i,l)*sin(x(i))/j(i,0)
  100                continue
                  endif
   70         continue
              if (der.eq.'derivatives') then
                  do 200 l=1,lmax
                     lm=l-1
                     do 300 i=1,np
                        jp(i,l)=j(i,lm)-j(i,l)/x(i)
  300                continue
  200             continue   
              endif
          endif   
      endif
      if (prnt) then
          if ( type.eq.'regular'.or.type.eq.'both') then      
               do 800 l=0,lmax
                  write (iout,900) l
                  write (iout,1000) (j(i,l),i=1,np)
  800          continue
               if (der.eq.'derivatives') then
                   do 810 l=0,lmax
                      write (iout,1010) l
                      write (iout,1000) (jp(i,l),i=1,np)
  810              continue
               endif
          endif               
          if ( type.eq.'irregular'.or.type.eq.'both') then      
               do 820 l=0,lmax
                  write (iout,1020) l
                  write (iout,1000) (y(i,l),i=1,np)
  820          continue
               if (der.eq.'derivatives') then
                   do 830 l=0,lmax
                      write (iout,1030) l
                      write (iout,1000) (yp(i,l),i=1,np)
  830              continue
               endif                                         
          endif
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
 1040 format(/,5x,'wronskian for l = ',1x,i3,1x,e15.8)
      return
      end
