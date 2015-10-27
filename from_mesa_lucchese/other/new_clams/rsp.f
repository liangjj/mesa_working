*deck %W%  %G%
      subroutine rsp(nm,n,nv,a,w,matz,z,fv1,fv2,ierr)                   
c***begin prologue  rsp                                                 
c***date written   760101   (yymmdd)                                    
c***revision date  830518   (yymmdd)                                    
c***category no.  d4a1                                                  
c***keywords  eigenvalues,eigenvectors,eispack                          
c***author  smith, b. t., et al.                                        
c***purpose  compute eigenvalues and, optionally, eigenvectors of       
c            real symmetric matrix packed into a one dimensional        
c***description                                                         
c                                                                       
c     this subroutine calls the recommended sequence of                 
c     subroutines from the eigensystem subroutine package (eispack)     
c     to find the eigenvalues and eigenvectors (if desired)             
c     of a real symmetric packed matrix.                                
c                                                                       
c     on input                                                          
c                                                                       
c        nm  must be set to the row dimension of the two-dimensional    
c        array parameters as declared in the calling program            
c        dimension statement.                                           
c                                                                       
c        n  is the order of the matrix  a.                              
c                                                                       
c        nv  is an integer variable set equal to the                    
c        dimension of the array  a  as specified for                    
c        a  in the calling program.  nv  must not be                    
c        less than  n*(n+1)/2.                                          
c                                                                       
c        a  contains the lower triangle of the real symmetric           
c        packed matrix stored row-wise.                                 
c                                                                       
c        matz  is an integer variable set equal to zero if              
c        only eigenvalues are desired.  otherwise it is set to          
c        any non-zero integer for both eigenvalues and eigenvectors.    
c                                                                       
c     on output                                                         
c                                                                       
c        w  contains the eigenvalues in ascending order.                
c                                                                       
c        z  contains the eigenvectors if matz is not zero.              
c                                                                       
c        ierr  is an integer output variable set equal to an            
c        error completion code described in section 2b of the           
c        documentation.  the normal completion code is zero.            
c                                                                       
c        fv1  and  fv2  are temporary storage arrays.                   
c                                                                       
c     questions and comments should be directed to b. s. garbow,        
c     applied mathematics division, argonne national laboratory         
c     ------------------------------------------------------------------
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow, 
c                 y. ikebe, v. c. klema, c. b. moler, *matrix eigen-    
c                 system routines - eispack guide*, springer-verlag,    
c                 1976.                                                 
c***routines called  tql2,tqlrat,trbak3,tred3                           
c***end prologue  rsp                                                   
c                                                                       
      implicit real*8 (a-h,o-z)
      integer i,j,n,nm,nv,ierr,matz                                     
      real*8 a(nv),w(n),z(nm,n),fv1(n),fv2(n)                           
c                                                                       
c***first executable statement  rsp                                     
      if (n .le. nm) go to 5                                            
      ierr = 10 * n                                                     
      go to 50                                                          
    5 if (nv .ge. (n * (n + 1)) / 2) go to 10                           
      ierr = 20 * n                                                     
      go to 50                                                          
c                                                                       
   10 call  tred3(n,nv,a,w,fv1,fv2)                                     
      if (matz .ne. 0) go to 20                                         
c     .......... find eigenvalues only ..........                       
      call  tqlrat(n,w,fv2,ierr)                                        
      go to 50                                                          
c     .......... find both eigenvalues and eigenvectors ..........      
   20 do 40 i = 1, n                                                    
c                                                                       
         do 30 j = 1, n                                                 
            z(j,i) = 0.0d0                                              
   30    continue                                                       
c                                                                       
         z(i,i) = 1.0d0                                                 
   40 continue                                                          
c                                                                       
      call  tql2(nm,n,w,fv1,z,ierr)                                     
      if (ierr .ne. 0) go to 50                                         
      call  trbak3(nm,n,nv,a,n,z)                                       
   50 return                                                            
      end                                                               
