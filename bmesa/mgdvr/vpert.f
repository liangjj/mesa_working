*deck vpert.f
c***begin prologue     vpert
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential matrix
c***                   
c***description        calculate the dependent potential
c***                   matrix elements in the dvr representation.
c***                   
c***references         
c
c***routines called    
c***end prologue       vpert
      subroutine vpert(pq1,pq2,pq3,pv1,pv2,pv3,pv,addv,n,nd,dim,
     1                 vword,prnt)
      implicit integer (a-z)
      real*8 q1, q2, q3, v1, v2, v3, v
      real*8 a11, a12, a13, a21, a22, a23, a31, a32, a33
      real*8 b11, b12, b13, b21, b22, b23, b31, b32, b33
      real*8 d11, d12, d13, d21, d22, d23, d31, d32, d33
      real*8 fpkey
      real*8 f11, f22, f33, f21, f31, f32
      real*8 depth
      character*32 type
      character*1600 card
      character*80 cpass, chrkey, title
      logical dollar, addv, prnt
      dimension nd(3), ngot(6)      
      common/io/inp, iout
      pointer (pq1,q1(1))
      pointer (pq2,q2(1))
      pointer (pq3,q3(1))
      pointer (pv1,v1(1))
      pointer (pv2,v2(1))
      pointer (pv3,v3(1))
      pointer (pv,v(1))
      pointer (pf11,f11(1))
      pointer (pf22,f22(1))
      pointer (pf33,f33(1))
      pointer (pf21,f21(1))
      pointer (pf31,f31(1))
      pointer (pf32,f32(1))
      need=wptoin(n)
      call memory(need,pv,vword,'vint',0)
      call rzero(v,n)
      if(dollar('$vpert',card,cpass,inp) ) then
         type=chrkey(card,'interaction-potential','none',' ')
         write(iout,1) type
      endif
      if(type.eq.'none') then
         return
      elseif(type.eq.'well') then
         if(dim.eq.1) then
            d11=fpkey(card,'d11',0.d0,' ')
            call vwel1(v,q1,d11,nd(1))
         elseif(dim.eq.2) then
            d11=fpkey(card,'d11',0.d0,' ')
            d22=fpkey(card,'d22',0.d0,' ')
            d21=fpkey(card,'d21',0.d0,' ')
            d12=fpkey(card,'d12',0.d0,' ')
            if(d12.ne.d21) then
               call lnkerr('error in off-diagonal well parameters')
            endif
            need=wptoin(nd(1))
            call memory(need,pf11,ngot(1),'f11',0)
            need=wptoin(nd(2))
            call memory(need,pf22,ngot(2),'f22',0)
            need=wptoin(nd(1)*nd(2))
            call memory(need,pf21,ngot(3),'f21',0)
            call vwel1(f11,q1,d11,nd(1))
            call vwel1(f22,q2,d22,nd(2))
            call vwel12(f21,q1,q2,d21,nd(1),nd(2))
            call vfill2(v,f11,f22,f21,nd(1),nd(2))
            call memory(-ngot(1),pf11,idum,'f11',idum)
            call memory(-ngot(2),pf22,idum,'f22',idum)
            call memory(-ngot(3),pf21,idum,'f21',idum)
         elseif(dim.eq.3) then
            d11=fpkey(card,'d11',0.d0,' ')
            d12=fpkey(card,'d12',0.d0,' ')
            d13=fpkey(card,'d13',0.d0,' ')
            d21=fpkey(card,'d21',0.d0,' ')
            d22=fpkey(card,'d22',0.d0,' ')
            d23=fpkey(card,'d23',0.d0,' ')
            d31=fpkey(card,'d31',0.d0,' ')
            d32=fpkey(card,'d32',0.d0,' ')
            d33=fpkey(card,'d33',0.d0,' ')
            if(d12.ne.d21.or.
     1         d13.ne.d31.or.
     2         d23.ne.d32) then
               call lnkerr('error in off-diagonal well parameters')
            endif
            need=wptoin(nd(1))
            call memory(need,pf11,ngot(1),'f11',0)
            need=wptoin(nd(2))
            call memory(need,pf22,ngot(2),'f22',0)
            need=wptoin(nd(3))
            call memory(need,pf33,ngot(3),'f33',0)
            need=wptoin(nd(2)*nd(1))
            call memory(need,pf21,ngot(4),'f21',0)
            need=wptoin(nd(3)*nd(1))
            call memory(need,pf31,ngot(5),'f31',0)
            need=wptoin(nd(3)*nd(2))
            call memory(need,pf32,ngot(6),'f32',0)
            call vwel12(f21,q2,q1,d21,nd(2),nd(1))
            call vwel12(f31,q3,q1,d31,nd(3),nd(1))
            call vwel12(f32,q3,q2,d32,nd(3),nd(2))
            call vfill3(v,f11,f22,f33,f21,f31,f32,nd(1),nd(2),nd(3))
         endif
         if(prnt) then
            title='well'
            call prntrm(title,v,n,1,n,1,iout)
         endif
      elseif(type.eq.'exponential') then
         if(dim.eq.1) then
            a11=fpkey(card,'a11',0.d0,' ')
            b11=fpkey(card,'b11',0.d0,' ')
            call vexp(v,q1,a11,b11,nd(1))
             if(addv) then
                call vsub(v,v,v1,nd(1))
             endif
             if(prnt) then
                title='1d exponential'
                call prntrm(title,v,nd(1),1,nd(1),1,iout)
             endif
         elseif(dim.eq.2) then
            a11=fpkey(card,'a11',0.d0,' ')
            b11=fpkey(card,'b11',0.d0,' ')
            a22=fpkey(card,'a22',0.d0,' ')
            b22=fpkey(card,'b22',0.d0,' ')
            a21=fpkey(card,'a21',0.d0,' ')
            b21=fpkey(card,'b21',0.d0,' ')
            a12=fpkey(card,'a12',0.d0,' ')
            b12=fpkey(card,'b12',0.d0,' ')
            if(a12.ne.a21.or.b12.ne.b21) then
               call lnkerr('error in off-diagonal exponential '//
     1                     'parameters')
            endif
            need=wptoin(nd(1))
            call memory(need,pf11,ngot(1),'f11',0)
            need=wptoin(nd(2))
            call memory(need,pf22,ngot(2),'f22',0)
            need=wptoin(nd(1)*nd(2))
            call memory(need,pf21,ngot(3),'f21',0)
            call vexp(f11,q1,a11,b11,nd(1))
            call vexp(f22,q2,a22,b22,nd(2))
            call vpair(f21,a21,b21,q2,q1,nd(2),nd(1))
            call vfill2(v,f11,f22,f21,nd(1),nd(2))
            call memory(-ngot(1),pf11,idum,'f11',idum)
            call memory(-ngot(2),pf22,idum,'f22',idum)
            call memory(-ngot(3),pf21,idum,'f21',idum)
            if(addv) then
               call vsub2(v,v1,v2,nd(1),nd(2))
            endif
            if(prnt) then
               title='2d exponential'
            endif   
         elseif(dim.eq.3) then
            a11=fpkey(card,'a11',0.d0,' ')
            b11=fpkey(card,'b11',0.d0,' ')
            a22=fpkey(card,'a22',0.d0,' ')
            b22=fpkey(card,'b22',0.d0,' ')
            a33=fpkey(card,'a33',0.d0,' ')
            b33=fpkey(card,'b33',0.d0,' ')
            a21=fpkey(card,'a21',0.d0,' ')
            b21=fpkey(card,'b21',0.d0,' ')
            a21=fpkey(card,'a12',0.d0,' ')
            b21=fpkey(card,'b12',0.d0,' ')
            a31=fpkey(card,'a31',0.d0,' ')
            b31=fpkey(card,'b31',0.d0,' ')
            a31=fpkey(card,'a13',0.d0,' ')
            b31=fpkey(card,'b13',0.d0,' ')
            a32=fpkey(card,'a32',0.d0,' ')
            b32=fpkey(card,'b32',0.d0,' ')
            a32=fpkey(card,'a23',0.d0,' ')
            b32=fpkey(card,'b23',0.d0,' ')
            a12=a21
            b12=b21
            a13=a31
            b13=b31
            a23=a32
            b23=b32
            need=wptoin(nd(1))
            call memory(need,pf11,ngot(1),'f11',0)
            need=wptoin(nd(2))
            call memory(need,pf22,ngot(2),'f22',0)
            need=wptoin(nd(3))
            call memory(need,pf33,ngot(3),'f33',0)
            need=wptoin(nd(2)*nd(1))
            call memory(need,pf21,ngot(4),'f21',0)
            need=wptoin(nd(3)*nd(1))
            call memory(need,pf31,ngot(5),'f31',0)
            need=wptoin(nd(3)*nd(2))
            call memory(need,pf32,ngot(6),'f32',0)
            call vexp(f11,q1,a11,b11,nd(1))
            call vexp(f22,q2,a22,b22,nd(2))
            call vexp(f33,q3,a33,b33,nd(3))
            call vpair(f21,a21,b21,q2,q1,nd(2),nd(1))
            call vpair(f31,a31,b31,q3,q1,nd(3),nd(1))
            call vpair(f32,a32,b32,q3,q2,nd(3),nd(2))
            call vfill3(v,f11,f22,f33,f21,f31,f32,nd(1),nd(2),nd(3))
            if(addv) then
               call vsub3(v,v1,v2,v3,nd(1),nd(2),nd(3))
            endif   
            call memory(-ngot(1),pf11,idum,'f11',idum)
            call memory(-ngot(2),pf22,idum,'f22',idum)
            call memory(-ngot(3),pf33,idum,'f33',idum)
            call memory(-ngot(4),pf21,idum,'f21',idum)
            call memory(-ngot(5),pf31,idum,'f31',idum)
            call memory(-ngot(6),pf31,idum,'f32',idum)
         else
            call lnkerr('error in dimension')
         endif
      elseif(type.eq.'2d-model') then
            a11 = -3.d0
            b11 = 1.d0
            a22 = a11
            b22 = b11
            a21 = 10.d0
            b21 = 1.d0
            need=wptoin(nd(1))
            call memory(need,pf11,ngot(1),'f11',0)
            need=wptoin(nd(2))
            call memory(need,pf22,ngot(2),'f22',0)
            need=wptoin(nd(1)*nd(2))
            call memory(need,pf21,ngot(3),'f21',0)
            call vexp(f11,q1,a11,b11,nd(1))
            call vexp(f22,q2,a22,b22,nd(2))
            call vexp2(f21,a21,b21,q2,q1,nd(2),nd(1))
            call vfill2(v,f11,f22,f21,nd(1),nd(2))
            call memory(-ngot(1),pf11,idum,'f11',idum)
            call memory(-ngot(2),pf22,idum,'f22',idum)
            call memory(-ngot(3),pf21,idum,'f21',idum)
            if(addv) then
               call vsub2(v,v1,v2,nd(1),nd(2))
            endif
      else
         call lnkerr('error in potential type')
      endif
      return
 1    format(/,5x,'interaction potential = ',a32)
      end       





