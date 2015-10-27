*deck honv.f
c***begin prologue     honv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           matrix vector
c***author             schneider, barry (nsf)
c***source             tproplib              
c***purpose            matrix vector multiply for dvr hamiltonian 
c***description        the vector vecin is transformed to the vector
c***                   vecout by operating with the hamiltonian.              
c
c                      n = n3*n2*n1
c                      h03(n3,n3) = matrix representation of T(3) + V0(3)
c                      h02(n2,n2) = matrix representation of T(2) + V0(2)
c                      h01(n1,n1) = matrix representation of T(1) + V0(1)
c                      vecin(n)   = vecin(n3,n2,n1)
c                      vecout(n)  = vecout(n3,n2,n1)
c
c                      the calling arguments must be passed appropriately
c                      for the one, two or three dimensional problem for
c                      proper operation.
c***references         
c***routines called    
c***end prologue       honv
      subroutine honv(h03,h02,h01,vt,vecin,vecout,n,n3,n2,n1,dim)
      implicit integer (a-z)
      real*8 h03, h02, h01, vt, vecin, vecout
      dimension h03(n3,n3), h02(n2,n2), h01(n1,n1)
      dimension vt(n), vecin(n), vecout(n)
      common/io/inp, iout
      if(dim.eq.1) then
c
c        in this case n3 = dimension of the single operator 
c
         call ebc(vecout,h03,vecin,n3,n3,1)
c
      elseif(dim.eq.2) then
c
c        in this case n3 and n2 dimension the vector
c      
         call ebc(vecout,h02,vecin,n3,n3,n2)
         call apbc(vecout,vecin,h01,n3,n2,n2)
      elseif(dim.eq.3) then             
         call ebc(vecout,h03,vecin,n3,n3,n2*n1)
         begin=1
         do 10 i=1,n1
            call apbc(vecout(begin),vecin(begin),h02,n3,n2,n2)
            begin=begin+n3*n2
 10      continue
         call apbc(vecout,vecin,h01,n3*n2,n1,n1)                       
      else
         call lnkerr('error in dimension')
      endif
c
c     now add in the diagonal
c
      call va2apbc(vecout,vt,vecin,n)                        
      return
      end       



