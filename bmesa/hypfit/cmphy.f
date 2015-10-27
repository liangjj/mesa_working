*deck cmphy.f
c***begin prologue     cmphy
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            compare derivatives in (r,alpha)
c***references         
c
c***routines called    
c***end prologue       cmphy
      subroutine cmphy(exact,approx,var,type,n1,n2)
      implicit integer (a-z)
      real*8 exact, approx, var 
      character*2 itoc
      character*80 title
      character*(*) type
      dimension exact(n2,n1), approx(n2,n1), var(n2,n1)
      common/io/inp, iout
      l1=length(type)
      title='Comparison of Exact and Approximate Derivative of '
     1            //type(1:l1)//' Function' 
      l2=length(title)
      write(iout,*) title(1:l2)
      write(iout,1)
      do 10 i=1,n1
         do 20 j=1,n2
            write(iout,2) var(j,i), exact(j,i), approx(j,i)
 20      continue   
 10   continue   
      return
 1    format(/,1x,'   grid point  ',1x,'     Exact     ',1x,
     1            '  Approximate  ')
 2    format(1x,d15.8,1x,d15.8,1x,d15.8)
      end       















