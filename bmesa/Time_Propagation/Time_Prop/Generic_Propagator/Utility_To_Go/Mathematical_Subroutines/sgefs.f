*deck  @(#)sgefs.f	1.2 8/8/91
      subroutine sgefs(a,lda,n,v,itask,ind,work,iwork)                  
c***begin prologue        sgefs 
c***date written          800317   (yymmdd)  
c***revision date         901711   (yymmdd) 
c    november 17, 1990  rlm at lanl
c    changed the xerror calls to lnkerr in order to get rid of
c    all the slatec xerror stuff that must be loaded.
c***category no.          d2a1  
c***keywords              general system of linear equations,linear equations
c***author                voorhees, e., (lanl)
c***purpose               sgefs solves a general single precision real
c                         nxn system of linear equations.
c***description        
c                                                                       
c    subroutine sgefs solves a general nxn system of single             
c    precision linear equations using linpack subroutines sgeco         
c    and sgesl.  that is, if a is an nxn real matrix and if x           
c    and b are real n-vectors, then sgefs solves the equation           
c                                                                       
c                          a*x=b.                                       
c                                                                       
c    the matrix a is first factored into upper and lower tri-           
c    angular matrices u and l using partial pivoting.  these            
c    factors and the pivoting information are used to find the          
c    solution vector x.  an approximate condition number is             
c    calculated to provide a rough estimate of the number of            
c    digits of accuracy in the computed solution.                       
c                                                                       
c    if the equation a*x=b is to be solved for more than one vector     
c    b, the factoring of a does not need to be performed again and      
c    the option to only solve (itask .gt. 1) will be faster for         
c    the succeeding solutions.  in this case, the contents of a,        
c    lda, n and iwork must not have been altered by the user follow-    
c    ing factorization (itask=1).  ind will not be changed by sgefs     
c    in this case.                                                      
c                                                                       
c  argument description ***                                             
c                                                                       
c    a      real(lda,n)                                                 
c             on entry, the doubly subscripted array with dimension     
c               (lda,n) which contains the coefficient matrix.          
c             on return, an upper triangular matrix u and the           
c               multipliers necessary to construct a matrix l           
c               so that a=l*u.                                          
c    lda    integer                                                     
c             the leading dimension of the array a.  lda must be great- 
c             er than or equal to n.  (terminal error message ind=-1)   
c    n      integer                                                     
c             the order of the matrix a.  the first n elements of       
c             the array a are the elements of the first column of       
c             the  matrix a.  n must be greater than or equal to 1.     
c             (terminal error message ind=-2)                           
c    v      real(n)                                                     
c             on entry, the singly subscripted array(vector) of di-     
c               mension n which contains the right hand side b of a     
c               system of simultaneous linear equations a*x=b.          
c             on return, v contains the solution vector, x .            
c    itask  integer                                                     
c             if itask=1, the matrix a is factored and then the         
c               linear equation is solved.                              
c             if itask .gt. 1, the equation is solved using the existing
c               factored matrix a and iwork.                            
c             if itask .lt. 1, then terminal error message ind=-3 is    
c               printed.                                                
c    ind    integer                                                     
c             gt. 0  ind is a rough estimate of the number of digits    
c                     of accuracy in the solution, x.                   
c             lt. 0  see error message corresponding to ind below.      
c    work   real(n)                                                     
c             a singly subscripted array of dimension at least n.       
c    iwork  integer(n)                                                  
c             a singly subscripted array of dimension at least n.       
c                                                                       
c  error messages printed ***                                           
c                                                                       
c    ind=-1  terminal   n is greater than lda.                          
c    ind=-2  terminal   n is less than 1.                               
c    ind=-3  terminal   itask is less than 1.                           
c    ind=-4  terminal   the matrix a is computationally singular.       
c                         a solution has not been computed.             
c    ind=-10 warning    the solution has no apparent significance.      
c                         the solution may be inaccurate or the matrix  
c                         a may be poorly scaled.                       
c                                                                       
c               note-  the above terminal(*fatal*) error messages are   
c                      designed to be handled by lnkerr.
c***references            subroutine sgefs was developed by group c-3, los alamos 
c                         scientific laboratory, los alamos, nm 87545.          
c                         the linpack subroutines used by sgefs are described in
c                         detail in the *linpack users guide* published by      
c                         the society for industrial and applied mathematics    
c                         (siam) dated 1979.                                    
c***routines called       sgeco,sgesl
c***end prologue          sgefs                                                 
c                                                                       
      integer lda,n,itask,ind,iwork(n)                                  
      real*8 a(lda,n),v(n),work(n)
      real*8 rcond,r1mach
c
      common/io/inp,iout
c***first executable statement  sgefs                                   
      if (lda.lt.n) then
         ind=-1
         call lnkerr('sgefs error (ind=-1) -- lda=i1 is less than n=i2')
         return 
      endif
      if (n.le.0) then
         ind=-2
         call lnkerr( 'sgefs error (ind=-2) -- n=i1 is less than 1')
         return
      endif
      if (itask.lt.1) then
         ind=-3
         call lnkerr( 'sgefs error (ind=-3) -- itask=i1 is less than 1')
         return
      endif
c
c
      if (itask.gt.1) then
c        matrix already factored
         call sgesl(a,lda,n,iwork,v,0)
      else 
c        factor matrix a into lu                                           
         call sgeco(a,lda,n,iwork,rcond,work)
c                                                                       
c        check for computationally singular matrix                         
         if (rcond.eq.0.0d+00) then
            ind=-4  
            call lnkerr( 'sgefs error (ind=-4) -- singular matrix a')
         endif
c                                                                       
c        compute ind (estimate of no. of significant digits)               
c           note that this is inherently machine-dependent.  r1mach is
c           the largest relative spacing, r1mach(4)=b**(1-t). this routine
c           is therefore somewhat machine dependent, and r1mach should be made
c           a function in mdutil.  however, it is just used for an estimate
c           of the precision, and so we will use the value appropriate for
c           real*8 on the sun. this should cover most of the representations
c           used on 32-bit word length machines. 
c           it should also be roughly appropriate for single precision on
c           the cray.
         r1mach=1.0d-14
         ind=-int(log10(r1mach/rcond))      
c                                                                       
c        check for ind greater than zero                                   
         if (ind.le.0) then 
c           warning... solution may not be significant 
            ind=-10
            write(iout,*) 
     $            'sgefs--ind=-10;solution may not be significant'
         endif
c        solve after factoring
         call sgesl(a,lda,n,iwork,v,0)   
      endif
c
c
      return                                                            
      end                                                               
