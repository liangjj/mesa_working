*deck vexp123.f
c***begin prologue     vexp123
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            additive pair potential
c***                   
c***references         
c
c***routines called    
c***end prologue       vexp123
      subroutine vexp123(v,q1,q2,q3,a12,a13,a23,
     1                              b12,b13,b23,n1,n2,n3,prn)
      implicit integer (a-z)
      real*8 v, q1, q2, q3, a12, a13, a23, b12, b13, b23, add
      logical prn
      character*80 title
      character*16 fptoc
      dimension v(n3,n2,n1), q1(n1), q2(n2), q3(n3)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            add = a12 * exp ( - b12 * abs( q2(j) - q1(i) ) )
            do 30 k=1,n3
               v(k,j,i) = v(k,j,i) + add
 30         continue
 20      continue   
 10   continue   
      do 40 i=1,n1
         do 50 k=1,n3
            add = a13 * exp ( - b13 * abs( q3(k) - q1(i) ) )
            do 60 j=1,n2
               v(k,j,i) = v(k,j,i) + add
 60         continue
 50      continue   
 40   continue   
      do 70 j=1,n2
         do 80 k=1,n3
            add = a23 * exp ( - b23 * abs( q3(k) - q2(j) ) )
            do 90 i=1,n1
               v(k,j,i) = v(k,j,i) + add
 90         continue
 80      continue   
 70   continue   
      if(prn) then
         do 100 i=1,n1
            title='sum of two dimensional exponentials x= '
     1                                         //fptoc(q1(i))
            call prntrm(title,v(1,1,i),n3,n2,n3,n2,iout)
 100     continue   
      endif
      return
      end       



