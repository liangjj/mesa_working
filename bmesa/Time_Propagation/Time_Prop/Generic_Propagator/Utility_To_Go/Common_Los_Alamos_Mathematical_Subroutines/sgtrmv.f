*deck  @(#)sgtrmv.f	1.1 8/6/91
      subroutine sgtrmv(d,subd,superd,x,y,n)                          
c***begin prologue  sgtrmv                                               
c***date written            (yymmdd)                                    
c***revision date  831207   (yymmdd)                                    
c***category no.  d1b6                                                  
c***keywords  matrix multiplication,tridiagonal, vector
c***author  barry schneider (nsf)                 
c***purpose  performs multiplication of a tridiagonal matrix times 
c***         a vector x to get a new vector y        
c            provided.                                                  
c***description                                                         

c     -arguments-                                                       
c                                                                       
c      on input                                                         
c
c       n = the matrix dimension.                     
c                                                                       
c       d = the diagonal elements of the tridiagonal matrix.
c
c       subd = the subdiagonal elements of the tridiagonal matrix.
c
c       superd = the superdiagonal elements of the tridiagonal matrix.
c
c       x = the input vector.                                           
c
c       y = the output vector.
c
c***references  (none)                                                  
c***routines called  saxpy                                              
c***end prologue  sgmv                                                 
      implicit real*8 (a-h,o-z)
      dimension d(n), subd(n), superd(n), x(n), y(n)
c***first executable statement  sgmv
      if(n.eq.1) then
         y(1) = d(1)*x(1)
      elseif(n.eq.2) then
         y(1) = d(1)*x(1) + superd(1)*x(2)
         y(2) = subd(2)*x(1) + d(2)*x(2)
      elseif(n.eq.3) then
         y(1) = d(1)*x(1) + superd(1)*x(2)
         y(2) = subd(2)*x(1) + d(2)*x(2) + superd(2)*x(3)
         y(3) = subd(3)*x(1) + d(3)*x(2)
      else
         y(1) = d(1)*x(1) + superd(1)*x(2)
         do 10 i=2,n-1
            y(i) = subd(i)*x(i-1) + d(i)*x(i) + superd(i)*x(i+1)
 10      continue
         y(n)=subd(n)*x(n-1) +d(n)*x(n)   
      endif
      return                                                            
      end                                                               

