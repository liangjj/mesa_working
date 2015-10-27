*deck  @(#)sgmv.f	1.1 8/6/91
      subroutine sgmv(m,n,a,ia,x,ix,y,iy,job)                          
c***begin prologue  sgmv                                               
c***date written            (yymmdd)                                    
c***revision date  831207   (yymmdd)                                    
c***category no.  d1b6                                                  
c***keywords  matrix multiplication,vector                              
c***author  jordan, tom, los alamos national laboratory                 
c***purpose  performs multiplication of a matrix a times a vector x or a
c            transposed times x.  four forms of accumulation are        
c            provided.                                                  
c***description                                                         
c                                                                       
c     given storage array a, this routine can be                        
c                      t                                                
c     used to compute a x=y or ax=y.  in addition, four forms           
c     of accumulation are provided for each form of transposition.      
c     these are ax=y, y+ax=y, -ax=y and y-ax=y.  this routine is        
c     designed to provide optimal efficiency on the cray-1 and          
c     near optimal efficiency on the cdc 7600.                          
c                                                                       
c     -arguments-                                                       
c                                                                       
c      on input                                                         
c       m = the number of rows of array a, if matrix a is to be used.   
c         = the number of columns of array a, if a(transpose)is used.   
c                                                                       
c       n = the number of columns of array a, if matrix a is used.      
c         = the number of rows of array a, if a(transpose) is used.     
c         = the number of elements of the vector x.                     
c                                                                       
c       a = the input matrix a, dimensioned a(ia,n).                    
c       ia = the leading dimension of a.                                
c       x = the input vector.                                           
c       ix = the increment between elements of x.                       
c       y = the product vector.  (output)                               
c       iy = the spacing between elements of the result vector.         
c       job = +1  y =   + a*x                                           
c             +2  y = y + a*x                                           
c             +3  y =   + a(transpose)*x                                
c             +4  y = y + a(transpose)*x                                
c                                                                       
c     if the job is negative, then subtraction is performed instead     
c     of addition.                                                      
c                                                                       
c      on output                                                        
c       y = the result vector as defined by job.  y may be a vector     
c           with arbitrary but uniform spacing, iy, of elements.        
c                                                                       
c***references  (none)                                                  
c***routines called  saxpy                                              
c***end prologue  sgmv                                                 
      implicit real*8 (a-h,o-z)
      dimension a(ia,n),x(ix,n),y(iy,n)                                 
c***first executable statement  sgmv                                   
      if (n .le. 0)                                 return              
      ii = 1                                                            
      ij = ia                                                           
      if (((iabs(job)-1)/2) .eq. 0)                      go to 1        
          ij = 1                                                        
          ii = ia                                                       
    1 continue                                                          
      if (mod(iabs(job)-1,2) .ne. 0)                    go to 4         
c                                                                       
ccc set y = 0.                                                          
c                                                                       
      do  2    j=1,m                                                    
          y(1,j) = 0.                                                   
    2 continue                                                          
    4 continue                                                          
c                                                                       
ccc accumulate columns of a to y                                        
c                                                                       
      if (job .lt. 0)                            go to 5                
      do  3    j=1,n                                                    
          call saxpy (m,x(1,j),a(1+(j-1)*ij,1),ii,y(1,1),iy)            
    3 continue                                                          
      return                                                            
    5 continue                                                          
      do  6  j=1,n                                                      
          call saxpy (m,-x(1,j),a(1+(j-1)*ij,1),ii,y(1,1),iy)           
    6 continue                                                          
      return                                                            
      end                                                               
