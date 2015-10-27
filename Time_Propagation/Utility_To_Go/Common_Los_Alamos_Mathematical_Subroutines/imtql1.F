*deck @@(#)imtql2.f	1.1  11/30/90
      subroutine imtql2(nm,n,d,e,z,ierr)                                
c***begin prologue          imtql2                                              
c***date written            760101   (yymmdd) 
c***revision date           830518   (yymmdd) 
c***category no.            d4a5,d4c2a  
c***keywords                eigenvalues,eigenvectors,eispack
c***author                  smith, b. t., et al. 
c***purpose                 computes eigenvalues and eigenvectors of symmetric
c                           tridiagonal matrix using implicit ql method. 
c***description                                                         
c                                                                       
c     this subroutine is a translation of the algol procedure imtql2,   
c     num. math. 12, 377-383(1968) by martin and wilkinson,             
c     as modified in num. math. 15, 450(1970) by dubrulle.              
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).   
c                                                                       
c     this subroutine finds the eigenvalues and eigenvectors            
c     of a symmetric tridiagonal matrix by the implicit ql method.      
c     the eigenvectors of a full symmetric matrix can also              
c     be found if  tred2  has been used to reduce this                  
c     full matrix to tridiagonal form.                                  
c                                                                       
c     on input                                                          
c                                                                       
c        nm must be set to the row dimension of two-dimensional         
c          array parameters as declared in the calling program          
c          dimension statement.                                         
c                                                                       
c        n is the order of the matrix.                                  
c                                                                       
c        d contains the diagonal elements of the input matrix.          
c                                                                       
c        e contains the subdiagonal elements of the input matrix        
c          in its last n-1 positions.  e(1) is arbitrary.               
c                                                                       
c        z contains the transformation matrix produced in the           
c          reduction by  tred2, if performed.  if the eigenvectors      
c          of the tridiagonal matrix are desired, z must contain        
c          the identity matrix.                                         
c                                                                       
c      on output                                                        
c                                                                       
c        d contains the eigenvalues in ascending order.  if an          
c          error exit is made, the eigenvalues are correct but          
c          unordered for indices 1,2,...,ierr-1.                        
c                                                                       
c        e has been destroyed.                                          
c                                                                       
c        z contains orthonormal eigenvectors of the symmetric           
c          tridiagonal (or full) matrix.  if an error exit is made,     
c          z contains the eigenvectors associated with the stored       
c          eigenvalues.                                                 
c                                                                       
c        ierr is set to                                                 
c          zero       for normal return,                                
c          j          if the j-th eigenvalue has not been               
c                     determined after 30 iterations.                   
c                                                                       
c     calls pythag(a,b) for sqrt(a**2 + b**2).                          
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c     ------------------------------------------------------------------
c***references              b. t. smith, j. m. boyle, j. j. dongarra,
c                           b. s. garbow, y. ikebe, v. c. klema, c. b. moler, 
c                           matrix eigen-system routines - 
c                           eispack guide*, springer-verlag, 1976.
c***routines called         pythag
c***end prologue            imtql2 
c                                                                       
      implicit real*8 (a-h,o-z)
      integer i,j,k,l,m,n,ii,nm,mml,ierr                                
      real*8 d(n),e(n),z(nm,n)                                          
      real*8 b,c,f,g,p,r,s,s1,s2                                        
      real*8 pythag                                                     
c                                                                       
c***first executable statement  imtql2                                  
      ierr = 0                                                          
      if (n .eq. 1) go to 1001                                          
c                                                                       
      do 100 i = 2, n                                                   
  100 e(i-1) = e(i)                                                     
c                                                                       
      e(n) = 0.0d0                                                      
c                                                                       
      do 240 l = 1, n                                                   
         j = 0                                                          
c     .......... look for small sub-diagonal element ..........         
  105    do 110 m = l, n                                                
            if (m .eq. n) go to 120                                     
            s1 = abs(d(m)) + abs(d(m+1))                                
            s2 = s1 + abs(e(m))                                         
            if (s2 .eq. s1) go to 120                                   
  110    continue                                                       
c                                                                       
  120    p = d(l)                                                       
         if (m .eq. l) go to 240                                        
         if (j .eq. 30) go to 1000                                      
         j = j + 1                                                      
c     .......... form shift ..........                                  
         g = (d(l+1) - p) / (2.0d0 * e(l))                              
         r = pythag(g,1.0d0)                                            
         g = d(m) - p + e(l) / (g + sign(r,g))                          
         s = 1.0d0                                                      
         c = 1.0d0                                                      
         p = 0.0d0                                                      
         mml = m - l                                                    
c     .......... for i=m-1 step -1 until l do -- ..........             
         do 200 ii = 1, mml                                             
            i = m - ii                                                  
            f = s * e(i)                                                
            b = c * e(i)                                                
            if (abs(f) .lt. abs(g)) go to 150                           
            c = g / f                                                   
            r = sqrt(c*c+1.0d0)                                         
            e(i+1) = f * r                                              
            s = 1.0d0 / r                                               
            c = c * s                                                   
            go to 160                                                   
  150       s = f / g                                                   
            r = sqrt(s*s+1.0d0)                                         
            e(i+1) = g * r                                              
            c = 1.0d0 / r                                               
            s = s * c                                                   
  160       g = d(i+1) - p                                              
            r = (d(i) - g) * s + 2.0d0 * c * b                          
            p = s * r                                                   
            d(i+1) = g + p                                              
            g = c * r - b                                               
c     .......... form vector ..........                                 
            do 180 k = 1, n                                             
               f = z(k,i+1)                                             
               z(k,i+1) = s * z(k,i) + c * f                            
               z(k,i) = c * z(k,i) - s * f                              
  180       continue                                                    
c                                                                       
  200    continue                                                       
c                                                                       
         d(l) = d(l) - p                                                
         e(l) = g                                                       
         e(m) = 0.0d0                                                   
         go to 105                                                      
  240 continue                                                          
c     .......... order eigenvalues and eigenvectors ..........          
      do 300 ii = 2, n                                                  
         i = ii - 1                                                     
         k = i                                                          
         p = d(i)                                                       
c                                                                       
         do 260 j = ii, n                                               
            if (d(j) .ge. p) go to 260                                  
            k = j                                                       
            p = d(j)                                                    
  260    continue                                                       
c                                                                       
         if (k .eq. i) go to 300                                        
         d(k) = d(i)                                                    
         d(i) = p                                                       
c                                                                       
         do 280 j = 1, n                                                
            p = z(j,i)                                                  
            z(j,i) = z(j,k)                                             
            z(j,k) = p                                                  
  280    continue                                                       
c                                                                       
  300 continue                                                          
c                                                                       
      go to 1001                                                        
c     .......... set error -- no convergence to an                      
c                eigenvalue after 30 iterations ..........              
 1000 ierr = l                                                          
 1001 return                                                            
      end                                                               
