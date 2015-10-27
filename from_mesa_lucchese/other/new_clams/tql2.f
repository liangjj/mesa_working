*deck %W%  %G%
      subroutine tql2(nm,n,d,e,z,ierr)                                  
c***begin prologue  tql2                                                
c***date written   760101   (yymmdd)                                    
c***revision date  830518   (yymmdd)                                    
c***category no.  d4a5,d4c2a                                            
c***keywords  eigenvalues,eigenvectors,eispack                          
c***author  smith, b. t., et al.                                        
c***purpose  compute eigenvalues and eigenvectors of symmetric          
c            tridiagonal matrix.                                        
c***description                                                         
c                                                                       
c     this subroutine is a translation of the algol procedure tql2,     
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and     
c     wilkinson.                                                        
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).   
c                                                                       
c     this subroutine finds the eigenvalues and eigenvectors            
c     of a symmetric tridiagonal matrix by the ql method.               
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
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow, 
c                 y. ikebe, v. c. klema, c. b. moler, *matrix eigen-    
c                 system routines - eispack guide*, springer-verlag,    
c                 1976.                                                 
c***routines called  pythag                                             
c***end prologue  tql2                                                  
c                                                                       
      implicit real*8 (a-h,o-z)
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr                          
      real*8 d(n),e(n),z(nm,n)                                          
      real*8 b,c,c2,c3,dl1,el1,f,g,h,p,r,s,s2                           
      real*8 pythag                                                     
c                                                                       
c***first executable statement  tql2                                    
      ierr = 0                                                          
      if (n .eq. 1) go to 1001                                          
c                                                                       
      do 100 i = 2, n                                                   
  100 e(i-1) = e(i)                                                     
c                                                                       
      f = 0.0d0                                                         
      b = 0.0d0                                                         
      e(n) = 0.0d0                                                      
c                                                                       
      do 240 l = 1, n                                                   
         j = 0                                                          
         h = abs(d(l)) + abs(e(l))                                      
         if (b .lt. h) b = h                                            
c     .......... look for small sub-diagonal element ..........         
         do 110 m = l, n                                                
            if (b + abs(e(m)) .eq. b) go to 120                         
c     .......... e(n) is always zero, so there is no exit               
c                through the bottom of the loop ..........              
  110    continue                                                       
c                                                                       
  120    if (m .eq. l) go to 220                                        
  130    if (j .eq. 30) go to 1000                                      
         j = j + 1                                                      
c     .......... form shift ..........                                  
         l1 = l + 1                                                     
         l2 = l1 + 1                                                    
         g = d(l)                                                       
         p = (d(l1) - g) / (2.0d0 * e(l))                               
         r = pythag(p,1.0d0)                                            
         d(l) = e(l) / (p + sign(r,p))                                  
         d(l1) = e(l) * (p + sign(r,p))                                 
         dl1 = d(l1)                                                    
         h = g - d(l)                                                   
         if (l2 .gt. n) go to 145                                       
c                                                                       
         do 140 i = l2, n                                               
  140    d(i) = d(i) - h                                                
c                                                                       
  145    f = f + h                                                      
c     .......... ql transformation ..........                           
         p = d(m)                                                       
         c = 1.0d0                                                      
         c2 = c                                                         
         el1 = e(l1)                                                    
         s = 0.0d0                                                      
         mml = m - l                                                    
c     .......... for i=m-1 step -1 until l do -- ..........             
         do 200 ii = 1, mml                                             
            c3 = c2                                                     
            c2 = c                                                      
            s2 = s                                                      
            i = m - ii                                                  
            g = c * e(i)                                                
            h = c * p                                                   
            if (abs(p) .lt. abs(e(i))) go to 150                        
            c = e(i) / p                                                
            r = sqrt(c*c+1.0d0)                                         
            e(i+1) = s * p * r                                          
            s = c / r                                                   
            c = 1.0d0 / r                                               
            go to 160                                                   
  150       c = p / e(i)                                                
            r = sqrt(c*c+1.0d0)                                         
            e(i+1) = s * e(i) * r                                       
            s = 1.0d0 / r                                               
            c = c * s                                                   
  160       p = c * d(i) - s * g                                        
            d(i+1) = h + s * (c * g + s * d(i))                         
c     .......... form vector ..........                                 
            do 180 k = 1, n                                             
               h = z(k,i+1)                                             
               z(k,i+1) = s * z(k,i) + c * h                            
               z(k,i) = c * z(k,i) - s * h                              
  180       continue                                                    
c                                                                       
  200    continue                                                       
c                                                                       
         p = -s * s2 * c3 * el1 * e(l) / dl1                            
         e(l) = s * p                                                   
         d(l) = c * p                                                   
         if (b + abs(e(l)) .gt. b) go to 130                            
  220    d(l) = d(l) + f                                                
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
