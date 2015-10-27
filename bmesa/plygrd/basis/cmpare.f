*deck cmpare.f
c***begin prologue     cmpare
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       cmpare
      subroutine cmpare(vin,vout,p1,p2,p3,p4,q1,q2,q3,q4,n1,n2,n3,n4,
     1                  nvc,n,dim,type,prnt)
      implicit integer (a-z)
      real*8 vin, vout, p1, p2, p3, p4, q1, q2, q3, q4, fac, fe
      character*80 title
      character*(*) type
      logical prnt
      dimension vin(n,nvc), vout(n,nvc)
      dimension q1(n1), q2(n2), q3(n4), q4(n4)
      dimension p1(n1,n1), p2(n2,n2), p3(n3,n3), p4(n4,n4)
      common/io/inp, iout
      if(prnt) then
         title='interpolated coefficients'
         call prntrm(title,vin,n,nvc,n,nvc,iout)
      endif
      if(dim.eq.1) then
         call v1di2o(vout,vin,p1,n1,nvc)
      elseif(dim.eq.2) then
         call v2di2o(vout,vin,p1,p2,n1,n2,nvc)
      elseif(dim.eq.3) then
         call v3di2o(vout,vin,p1,p2,p3,n1,n2,n3,nvc)
      elseif(dim.eq.4) then
         call v4di2o(vout,vin,p1,p2,p3,p4,n1,n2,n3,n4,nvc)
      endif
      if(type.eq.'sine') then
         do 10 i=1,n1
            fe = sin(q1(i))
            write(iout,1) q1(i), fe, vin(i,1)
 10      continue   
      elseif(type.eq.'cosine') then
         do 20 i=1,n1
            fe = cos(q1(i))
            write(iout,1) q1(i), fe, vin(i,1)
 20      continue   
      elseif(type.eq.'quadratic') then
         do 30 i=1,n1
            fe = q1(i)*q1(i)
            write(iout,1) q1(i), fe, vin(i,1)
 30      continue   
      elseif(type.eq.'exponential') then
         do 40 i=1,n1
            fe = exp(-q1(i))
            write(iout,1) q1(i), fe, vin(i,1)
 40      continue   
      elseif(type.eq.'gaussian') then
         do 50 i=1,n1
            fe = exp(-q1(i)*q1(i))
            write(iout,1) q1(i), fe, vin(i,1)
 50      continue   
      else
         call lnkerr('error in function type')
      endif
      return
 1    format(/,1x,'x = ',e15.8,1x,'exact = ',e15.8,1x,
     1                            'approximate = ',e15.8)
      end
