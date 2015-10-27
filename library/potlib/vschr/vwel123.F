*deck vwel123.f
c***begin prologue     vwel123
c***date written       000619   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            two dimensional potential well
c***                   
c***references         
c
c***routines called    
c***end prologue       vwel123
      subroutine vwel123(v,q1,q2,q3,a12,a13,a23,
     1                              d12,d13,d23,n1,n2,n3,prn)
      implicit integer (a-z)
      real*8 v, q1, q2, q3, a12, a13, a23, d12, d13, d23
      logical prn
      character*80 title
      character*16 fptoc
      dimension v(n3,n2,n1), q1(n1), q2(n2), q3(n3)
      common/io/inp, iout
      do 10 i=1,n1
         do 20 j=1,n2
            if( abs(q1(i)-q2(j)).le.a12) then
                do 30 k=1,n3
                   v(k,j,i) = v(k,j,i) + d12
 30             continue   
            endif
 20      continue   
 10   continue   
      do 40 i=1,n1
         do 50 k=1,n3
            if( abs(q1(i)-q3(k)).le.a13) then
                do 60 j=1,n2
                   v(k,j,i) = v(k,j,i) + d13
 60             continue   
            endif
 50      continue   
 40   continue   
      do 70 j=1,n2
         do 80 k=1,n3
            if( abs(q2(j)-q3(k)).le.a23) then
                do 90 i=1,n1
                   v(k,j,i) = v(k,j,i) + d23
 90             continue   
            endif
 80      continue   
 70   continue   
      if(prn) then
         do 100 i=1,n1
            title='sum of two dimensional wells x= '//fptoc(q1(i))
            call prntrm(title,v(1,1,i),n3,n2,n3,n2,iout)
 100     continue   
      endif
      return
 1    format(/,9x,'coordinate',8x,'potential')
 2    format(5x,e15.8,3x,e15.8)
      end       



