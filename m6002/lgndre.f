*deck @(#)lgndre.f	1.1 9/7/91
c***begin prologue     lgndre
c***date written                (yymmdd)
c***revision date      901003   (yymmdd)
c***keywords           lgndre, link 6002
c***authors            unknown
c***                   
c***source             m6002
c***purpose            cartesian spherical harmonics
c***references       
c
c***routines called    
c***end prologue      lgndre
      subroutine lgndre(ylm,grid,npt,l,m)
      implicit integer (a-z)
      real*8 facd, rval, ylm, grid
      dimension grid(4,npt), ylm(npt), facd(3)
      common /io/ inp, iout
      data lmax, mmax /2,2/
      facd(1)=.5d+00
      facd(2)=3.d+00
      facd(3)=6.d+00
c----------------------------------------------------------------------c
c           calculate cartesian spherical harmonics                    c
c                           for                                        c
c                       small l and m                                  c
c                    * not normalized                                  c
c----------------------------------------------------------------------c
 
      if (l.gt.lmax) then
         call lnkerr('l greater than lmax:quit')
      endif
      if (l.eq.0) then
c----------------------------------------------------------------------c
c               s wave                                                 c
c----------------------------------------------------------------------c
          do 10 i=1,npt
             ylm(i)=1.d+00
   10     continue
      endif   
      if(l.eq.1) then
c----------------------------------------------------------------------c
c               p wave                                                 c
c----------------------------------------------------------------------c
         if (m.eq.0) then 
            do 20 i=1,npt
               rval=sqrt(grid(1,i)*grid(1,i)
     1              +grid(2,i)*grid(2,i) 
     2              +grid(3,i)*grid(3,i))
               ylm(i)=grid(3,i)/rval
   20       continue
         endif
         if (m.eq.-1) then
            do 30 i=1,npt
               rval=sqrt(grid(1,i)*grid(1,i)
     1              +grid(2,i)*grid(2,i) 
     2              +grid(3,i)*grid(3,i))
               ylm(i)=grid(1,i)/rval
   30       continue 
         endif
         if (m.eq.1) then                
            do 40 i=1,npt
               rval=sqrt(grid(1,i)*grid(1,i)
     1              +grid(2,i)*grid(2,i) 
     2              +grid(3,i)*grid(3,i))
               ylm(i)=grid(2,i)/rval
   40       continue  
         endif 
      endif
      if(l.eq.2) then 
c----------------------------------------------------------------------c 
c            d wave                                                    c
c----------------------------------------------------------------------c
        if (m.eq.0) then
           do 50 i=1,npt
              rval=grid(1,i)*grid(1,i)
     1             +grid(2,i)*grid(2,i) 
     2             +grid(3,i)*grid(3,i)        
              ylm(i)=(2.d+00*grid(3,i)*grid(3,i)
     1                -grid(1,i)*grid(1,i)
     2                -grid(2,i)*grid(2,i))*facd(1)/rval
   50      continue
        endif
        if (m.eq.-1) then                     
           do 60 i=1,npt
              rval=grid(1,i)*grid(1,i)
     1             +grid(2,i)*grid(2,i) 
     2             +grid(3,i)*grid(3,i)        
              ylm(i)=grid(1,i)*grid(3,i)*facd(2)/rval
   60      continue
        endif
        if (m.eq.1) then                    
           do 70 i=1,npt
              rval=grid(1,i)*grid(1,i)
     1             +grid(2,i)*grid(2,i) 
     2             +grid(3,i)*grid(3,i)        
              ylm(i)=grid(2,i)*grid(3,i)*facd(2)/rval
   70      continue
        endif
        if (m.eq.-2) then
           do 80 i=1,npt
              rval=grid(1,i)*grid(1,i)
     1             +grid(2,i)*grid(2,i) 
     2             +grid(3,i)*grid(3,i)        
              ylm(i)=(grid(1,i)*grid(1,i)-grid(2,i)*grid(2,i))
     1                                    *facd(2)/rval
   80      continue
        endif
        if (m.eq.2) then
           do 90 i=1,npt
              rval=grid(1,i)*grid(1,i)
     1             +grid(2,i)*grid(2,i) 
     2             +grid(3,i)*grid(3,i)        
              ylm(i)=grid(1,i)*grid(2,i)*facd(3)/rval
   90      continue
        endif
      endif
      return
      end
