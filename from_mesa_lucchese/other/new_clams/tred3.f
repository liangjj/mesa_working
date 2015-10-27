*deck %W%  %G%
      subroutine tred3(n,nv,a,d,e,e2)                                   
c***begin prologue  tred3                                               
c***date written   760101   (yymmdd)                                    
c***revision date  830518   (yymmdd)                                    
c***category no.  d4c1b1                                                
c***keywords  eigenvalues,eigenvectors,eispack                          
c***author  smith, b. t., et al.                                        
c***purpose  reduce real symmetric matrix stored in packed form to      
c            symmetric tridiagonal matrix using orthogonal              
c            transformations.                                           
c***description                                                         
c                                                                       
c     this subroutine is a translation of the algol procedure tred3,    
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
c                                                                       
c     this subroutine reduces a real symmetric matrix, stored as        
c     a one-dimensional array, to a symmetric tridiagonal matrix        
c     using orthogonal similarity transformations.                      
c                                                                       
c     on input                                                          
c                                                                       
c        n is the order of the matrix.                                  
c                                                                       
c        nv must be set to the dimension of the array parameter a       
c          as declared in the calling program dimension statement.      
c                                                                       
c        a contains the lower triangle of the real symmetric            
c          input matrix, stored row-wise as a one-dimensional           
c          array, in its first n*(n+1)/2 positions.                     
c                                                                       
c     on output                                                         
c                                                                       
c        a contains information about the orthogonal                    
c          transformations used in the reduction.                       
c                                                                       
c        d contains the diagonal elements of the tridiagonal matrix.    
c                                                                       
c        e contains the subdiagonal elements of the tridiagonal         
c          matrix in its last n-1 positions.  e(1) is set to zero.      
c                                                                       
c        e2 contains the squares of the corresponding elements of e.    
c          e2 may coincide with e if the squares are not needed.        
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c     ------------------------------------------------------------------
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow, 
c                 y. ikebe, v. c. klema, c. b. moler, *matrix eigen-    
c                 system routines - eispack guide*, springer-verlag,    
c                 1976.                                                 
c***routines called  (none)                                             
c***end prologue  tred3                                                 
c                                                                       
      implicit real*8 (a-h,o-z)
      integer i,j,k,l,n,ii,iz,jk,nv                                     
      real*8 a(nv),d(n),e(n),e2(n)                                      
      real*8 f,g,h,hh,scale                                             
c                                                                       
c     .......... for i=n step -1 until 1 do -- ..........               
c***first executable statement  tred3                                   
      do  300 ii = 1, n                                                 
         i = n + 1 - ii                                                 
         l = i - 1                                                      
         iz = (i * l) / 2                                               
         h = 0.0d0                                                      
         scale = 0.0d0                                                  
         if (l .lt. 1) go to 130                                        
c     .......... scale row (algol tol then not needed) ..........       
         do 120 k = 1, l                                                
            iz = iz + 1                                                 
            d(k) = a(iz)                                                
            scale = scale + abs(d(k))                                   
  120    continue                                                       
c                                                                       
         if (scale .ne. 0.0d0) go to 140                                
  130    e(i) = 0.0d0                                                   
         e2(i) = 0.0d0                                                  
         go to 290                                                      
c                                                                       
  140    do 150 k = 1, l                                                
            d(k) = d(k) / scale                                         
            h = h + d(k) * d(k)                                         
  150    continue                                                       
c                                                                       
         e2(i) = scale * scale * h                                      
         f = d(l)                                                       
         g = -sign(sqrt(h),f)                                           
         e(i) = scale * g                                               
         h = h - f * g                                                  
         d(l) = f - g                                                   
         a(iz) = scale * d(l)                                           
         if (l .eq. 1) go to 290                                        
         f = 0.0d0                                                      
c                                                                       
         do 240 j = 1, l                                                
            g = 0.0d0                                                   
            jk = (j * (j-1)) / 2                                        
c     .......... form element of a*u ..........                         
            do 180 k = 1, l                                             
               jk = jk + 1                                              
               if (k .gt. j) jk = jk + k - 2                            
               g = g + a(jk) * d(k)                                     
  180       continue                                                    
c     .......... form element of p ..........                           
            e(j) = g / h                                                
            f = f + e(j) * d(j)                                         
  240    continue                                                       
c                                                                       
         hh = f / (h + h)                                               
         jk = 0                                                         
c     .......... form reduced a ..........                              
         do 260 j = 1, l                                                
            f = d(j)                                                    
            g = e(j) - hh * f                                           
            e(j) = g                                                    
c                                                                       
            do 260 k = 1, j                                             
               jk = jk + 1                                              
               a(jk) = a(jk) - f * e(k) - g * d(k)                      
  260    continue                                                       
c                                                                       
  290    d(i) = a(iz+1)                                                 
         a(iz+1) = scale * sqrt(h)                                      
  300 continue                                                          
c                                                                       
      return                                                            
      end                                                               
