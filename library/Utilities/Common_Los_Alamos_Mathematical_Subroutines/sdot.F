*deck @(#)sdot.F	5.1  11/29/95
      function sdot(n,sx,incx,sy,incy)                           
c***begin prologue  sdot                                                
c***date written   791001   (yymmdd)                                    
c***revision date  820801   (yymmdd)                                    
c***category no.  d1a4                                                  
c***keywords  blas,inner product,linear algebra,vector                  
c***author  lawson, c. l., (jpl)                                        
c           hanson, r. j., (snla)                                       
c           kincaid, d. r., (u. of texas)                               
c           krogh, f. t., (jpl)                                         
c***purpose  s.p. inner product of s.p. vectors                         
c***description                                                         
c                                                                       
c                b l a s  subprogram                                    
c    description of parameters                                          
c                                                                       
c     --input--                                                         
c        n  number of elements in input vector(s)                       
c       sx  single precision vector with n elements                     
c     incx  storage spacing between elements of sx                      
c       sy  single precision vector with n elements                     
c     incy  storage spacing between elements of sy                      
c                                                                       
c     --output--                                                        
c     sdot  single precision dot product (zero if n .le. 0)             
c                                                                       
c     returns the dot product of single precision sx and sy.            
c     sdot = sum for i = 0 to n-1 of  sx(lx+i*incx) * sy(ly+i*incy),    
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is       
c     defined in a similar way using incy.                              
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  sdot                                                  
c                                                                       
      implicit real*8 (a-h,o-z)
      real*8 sdot
      real*8 sx(1),sy(1)                                                
c***first executable statement  sdot                                    
#ifdef LOCALBLAS
      real*8 temp
      sdot = 0.0d0
      temp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        temp = temp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = temp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        temp = temp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        temp = temp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = temp
      return
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      real*8 ddot
      sdot=ddot(n,sx,incx,sy,incy)
      return                                                            
      end                         
#endif MACHINENBLAS                                      
