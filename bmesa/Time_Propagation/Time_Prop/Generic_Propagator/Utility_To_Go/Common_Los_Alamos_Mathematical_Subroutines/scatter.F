*deck @(#)scatter.f	5.1 11/6/94
      subroutine scatter(n,a,index,b)
c***begin prologue     scatter
c***date written       yymmdd  
c***revision date      yymmdd
c
c***keywords           
c***author             unknown 
c***source             @(#)scatter.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       scatter
c
      implicit integer (a-z)
c
      real*8 a(*),b(n)
      integer index(n)
c
      do 1 i=1,n
         a(index(i))=b(i)
    1 continue
c
      return
      end
