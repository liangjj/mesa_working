*deck ab3d.f
c***begin prologue     ab3d
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           transformation
c***author             schneider, barry (nsf)
c***source             pdelib              
c***purpose            transformation of representations 
c***description                          
c***references         
c***routines called    
c***end prologue       ab2d
      subroutine ab3d(u3,u2,u1,vin,vout,t3,t2,n3,n2,n1,dim)
      implicit integer (a-z)
      real*8 u3, u2, u1, vin, vout, t3, t2
      dimension u3(n3,n3), u2(n2,n2), u1(n1,n1)
      dimension vin(n3,n2,n1), vout(n3,n2,n1)
      dimension t3(n3,*), t2(n3,n2,n1)
      common/io/inp, iout
      if(dim.eq.1) then
         call ebtc(vout,u3,vin,n3,n3,1)
      elseif(dim.eq.2) then
          call ebtc(t3,u3,vin,n3,n3,n2)
          call ebc(vout,t3,u2,n3,n2,n2)
      elseif(dim.eq.3) then             
          do 10 i=1,n1
             call ebtc(t3,u3,vin(1,1,i),n3,n3,n2)
             call ebc(t2(1,1,i),t3,u2,n3,n2,n2)
 10      continue
         call ebc(vout,t2,u1,n3*n2,n1,n1)          
      else
         call lnkerr('error in dimension')
      endif            
      return
      end       



