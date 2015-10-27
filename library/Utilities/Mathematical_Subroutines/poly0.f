*deck poly0.f
c***begin prologue     poly0
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            compute values or first or second derivatives of
c***                   (f(x)-a)**n * (f(x)-b)**m
c***description        
c***                     
c***                                                          
c***references         
c
c***routines called    
c***end prologue       poly0
      subroutine poly0(p,x,left,right,nl,nr,der,npts)
      implicit integer (a-z)
      real*8 p, x, left, right, t1, t2, t3, t4, t5, t6
      dimension p(npts), x(npts)
      common/io/inp, iout
      if(der.eq.0) then
         if(nl.eq.0.and.nr.eq.0) then
            call vfill(p,1.d0,npts)
         elseif(nl.ne.0.and.nr.eq.0) then
            do 10 i=1,npts
               p(i) = ( x(i) - left )**nl
 10         continue
         elseif(nl.eq.0.and.nr.ne.0) then
            do 20 i=1,npts
               p(i) = ( x(i) - right )**nr
 20         continue
         elseif(nl.ne.0.and.nr.ne.0) then
            do 30 i=1,npts
               p(i) = ( ( x(i) - right )**nr )*( ( x(i) - left )**nl )
 30         continue
         endif
      elseif(der.eq.1) then
         if(nl.eq.0.and.nr.eq.0) then
            call rzero(p,npts)
         elseif(nl.eq.0.and.nr.ne.0) then
            do 40 i=1,npts
               p(i) = nr*( x(i) - right )**(nr-1)
 40         continue   
         elseif(nl.ne.0.and.nr.eq.0) then
            do 50 i=1,npts
               p(i) = nl*( x(i) - left )**(nl-1)
 50         continue   
         elseif(nl.ne.0.and.nr.ne.0) then
            do 60 i=1,npts
               t1=( x(i) - left )**(nl-1)
               t2=( x(i) - right )**(nr-1)
               t3=t1*( x(i) - left )
               t4=t2*( x(i) - right )
               p(i) = nl*t1*t4+nr*t2*t3
 60         continue
         endif 
      elseif(der.eq.2) then
         if(nl.eq.0.and.nr.eq.0) then
            call rzero(p,npts)
         elseif(nl.ne.0.and.nr.eq.0) then
            if(nl.eq.1) then
               call rzero(p,npts)
            else
               do 70 i=1,npts
                  p(i) = nl*(nl-1)* ( (x(i)-left) )**(nl-2)
 70            continue
            endif
         elseif(nl.eq.0.and.nr.ne.0) then
            if(nr.eq.1) then
               call rzero(p,npts)
            else
               do 80 i=1,npts
                  p(i) = nr*(nr-1)* ( (x(i)-right) )**(nr-2)
 80            continue
            endif
         elseif(nl.ne.0.and.nr.ne.0) then
               if(nl.eq.1.and.nr.eq.1) then
                  call vfill(p,2.d0,npts)
               elseif(nl.eq.1.and.nr.gt.1) then
                  do 90 i=1,npts
                     p(i) = 2.d0*nr*( x(i)-right )**(nr-1)
     1                             +
     2                      nr*(nr-1)*( x(i)-left ) *
     3                                ( x(i)-right )**(nr-2)
 90               continue
               elseif(nl.gt.1.and.nr.eq.1) then
                  do 100 i=1,npts
                     p(i) = 2.d0*nl*( x(i)-left )**(nl-1)
     1                             +
     2                      nl*(nl-1)*( x(i)-right )*
     3                                ( x(i)-left )**(nl-2)
 100              continue
               elseif(nr.gt.1.and.nl.gt.1) then   
                  do 200 i=1,npts 
                     t1 = ( x(i)-left )**(nl-2)
                     t2 = ( x(i)-right )**(nr-2)
                     t3 = t1*( x(i)-left )
                     t4 = t2*( x(i)-right )
                     t5 = t3*( x(i)-left )
                     t6 = t4*( x(i)-right )
                     p(i) = nl*(nl-1)*t1*t6
     1                            +
     2                      2.d0*nl*nr*t3*t4
     3                            +
     4                      nr*(nr-1)*t2*t5
 200              continue
               endif
         endif
      endif
      return
      end       

