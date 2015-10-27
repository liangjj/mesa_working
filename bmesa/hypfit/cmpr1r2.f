*deck cmpr1r2.f
c***begin prologue     cmpr1r2
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            compare and print derivatives wrt (r1,r2)
c***references         
c
c***routines called    
c***end prologue       cmpr1r2
      subroutine cmpr1r2(exact,approx,r1,r2,type,coord,n1,n2)
      implicit integer (a-z)
      real*8 exact, approx, r1, r2 
      character*80 title
      character*2 itoc, mstr
      character*(*) type, coord
      dimension exact(n2,n1), approx(n2,n1), r1(n1), r2(n2)
      common/io/inp, iout
      mstr=itoc(m)
      l1=length(type)
      title='Comparison of Exact and Approximate '//coord(1:2)//
     1      ' Derivative of '//type(1:l1)//' Function'
      l2=length(title)
      write(iout,*) title(1:l2)
      write(iout,1)
      do 10 i=1,n1
         do 20 j=1,n2
            write(iout,2) r1(i), r2(j), exact(j,i), approx(j,i)
 20      continue   
 10   continue   
      return
 1    format(/,1x,'      R1      ',1x,'      R2      ',
     1         1x,'     Exact     ',1x,'  Approximate  ')
 2    format(1x,d15.8,1x,d15.8,1x,d15.8,1x,d15.8)
      end       















