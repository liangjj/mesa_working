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
      subroutine vpert(v4,x1,x2,x3,v1,v2,v3,n,dim,key,vword,
     1                 zeroit,prnt)
      implicit integer (a-z)
      real*8 x1, x2, x3, v1, v2, v3, v4
      real*8 a, b, d, len, dif
      real*8 fpkey, depth
      character*32 type
      character*1600 card
      character*80 cpass, chrkey, title
      character*(*) key
      character*1 itoc
      character*24 phr
      logical dollar, prnt, zeroit
      dimension n(dim+1)      
      dimension a(3,3), b(3,3), d(3,3), len(3,3), dif(3,3), phr(4)
      dimension x1(n(1)), x2(n(2)), x3(n(3))
      dimension v1(n(1)), v2(n(2)), v3(n(3)), v4(n(dim+1))
      common/io/inp, iout
      if(zeroit) then
         call rzero(v4,n(dim+1))
      endif
      pv(dim+1)=pv4
      if(dim.eq.1) then
         return
      endif 
c
c     assign variables and locations
c
      px1=px(1)
      pv1=pv(1)
      if(dim.gt.1) then
         px2=px(2)
         pv2=pv(2)
      endif
      if(dim.gt.2) then
         px3=px(3)
         pv3=pv(3)
      endif
      if(dollar(key,card,cpass,inp) ) then
         type=chrkey(card,'potential','none',' ')
          prnt=logkey(card,'print=potential',.false.,' ')
         write(iout,1) type
      endif
      if(type.eq.'none') then
         return
      elseif(type.eq.'well') then
         if(dim.ge.2) then
            dif(1,2)=abs( x2(n(2)) - x1(n(1)) )
            dif(2,1)=dif(1,2)
         endif
         if(dim.eq.3) then
            dif(1,3)=abs( x1(n(1)) - x3(n(3)) )
            dif(3,1)=dif(1,3)
            dif(2,3)=abs( x2(n(2)) - x3(n(3)) )
            dif(3,2)=dif(2,3)
         endif 
         do 10 i=1,dim
            phr(1)='d'//itoc(i)
            phr(2)='l'//itoc(i)
            do 20 j=1,dim
               phr(3)=phr(1)(1:2)//itoc(j)
               phr(4)=phr(2)(1:2)//itoc(j)
               d(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
               len(i,j)=fpkey(card,phr(4)(1:3),dif(i,j),' ')
 20         continue
 10      continue
         do 30 i=1,dim
            do 40 j=1,i
               if(d(i,j).ne.d(j,i).or.len(i,j).ne.len(j,i)) then
                  call lnkerr('error in off-diagonal well parameters')
               endif
 40         continue
 30      continue   
         if(dim.eq.2) then   
            call vwel12(v4,x1,x2,len(2,1),d(2,1),n(1),n(2),.false.)
         endif
         if(dim.eq.3) then
            call vwel123(v4,x1,x2,x3,len(2,1),len(3,1),len(3,2),
     1                   d(2,1),d(3,1),d(3,2),n(1),n(2),n(3),.false.)
         endif
         if(prnt) then
            title='well'
            call prntrm(title,v4,n(dim+1),1,n(dim+1),1,iout)
         endif
      elseif(type.eq.'exponential') then
         do 50 i=1,dim
            phr(1)='a'//itoc(i)
            phr(2)='b'//itoc(i)
            do 60 j=1,dim
               phr(3)=phr(1)(1:2)//itoc(j)
               phr(4)=phr(2)(1:2)//itoc(j)
               a(i,j)=fpkey(card,phr(3)(1:3),0.d0,' ')
               b(i,j)=fpkey(card,phr(4)(1:3),0.d0,' ')
 60         continue
 50      continue
         do 70 i=1,dim
            do 80 j=1,i   
               if(a(i,j).ne.a(j,i).or.b(i,j).ne.b(j,i)) then
                  call lnkerr('error in off-diagonal exponential '//
     1                        'parameters')
               endif
 80         continue
 70      continue   
         if(dim.eq.2) then
            call vexp12(v4,x1,x2,a(1,2),b(1,2),n(1),n(2),.false.)
            if(prnt) then
               title='2d exponential'
            endif   
         endif
         if(dim.eq.3) then
            call vexp123(v4,x1,x2,x3,a(1,2),a(1,3),a(2,3),
     1                               b(1,2),b(1,3),b(2,3),
     2                               n(1),n(2),n(3),.false.)
         endif
      else
         call lnkerr('error in potential type')
      endif
      return
 1    format(/,5x,'interaction potential = ',a32)
      end       














