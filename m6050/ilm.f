*deck @(#)ilm.f	
c***begin prologue     ilm
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           ilm, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            special function for optical potential
c***                   construction for two electron systems
c***references       
c
c***routines called    
c***end prologue       ilm
      function ilm (l,m,a,b,fact,prnt)
      implicit real *8 (a-h,o-z)
      common /io/ inp, iout
      logical prnt
      complex *16 ilm, a, b, c, sum, pre
      character *80 title
      character *2 itoc
      dimension fact(0:100)
c
      c=a+b
      l1=l+1
      m1=m+1 
      l2=l+2
      m2=m+2
      l3=l+3
      m3=m+3
      ilm = fact(l2)*fact(m1)/(a**l3)
      ilm = ilm * ( 1.d0/b**m2 - 1.d0/c**m2 )
      pre = 1.d0 /( ( a*a )*( c**(l+m+3) ) )
      sum = pre*fact(l1+m1)/fact(l1)
      if (l.gt.0) then
          do 10 isum=2,l1
              pre = pre * c / a
              sum = sum + isum * pre * fact(l+m+3-isum) / fact(l2-isum)
   10     continue
      endif
      ilm =2.d0*( ( ilm - fact(l1) * sum ) )
      if (prnt) then
          title='ilm-'//itoc(l)//'-'//itoc(m)
      endif
   20 format(2x,e15.8)   
      return
      end
