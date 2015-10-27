*deck h02dvr.f
c***begin prologue     h02dvr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           transformation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       h02dvr
      subroutine h02dvr(resid,vec,u01,u02,u03,t1,t2,t3,nd,n,m,dim)
      implicit integer (a-z)
      real*8 resid, vec, u01, u02, u03, t1, t2, t3
      character*80 title
      dimension nd(3)
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension vec(n,m), resid(n,m), t1(nd(2),nd(1)), t2(nd(2),nd(1))
      dimension t3(nd(1),nd(2)) 
      common/io/inp, iout
      n1=nd(1)
      n2=nd(2)
      n3=nd(3) 
      if(dim.eq.1) then
         call ebc(resid,u01,vec,n,n1,m)
      elseif(dim.eq.2) then
c
c         transform residual from h0 to dvr representation
c 
         do 10 i=1,m
            do 20 j=1,n1
               jj=n2*(j-1)                
               do 30 k=1,n2
                  t1(k,j)=vec(jj+k,i)
 30            continue
 20         continue
            call ebc(t2,u02,t1,n2,n2,n1)
            call ebct(t3,u01,t2,n1,n1,n2)
            count=0
            do 40 j=1,n1
               do 50 k=1,n2
                  count=count+1
                  resid(count,i) = t3(j,k)
 50            continue
 40         continue   
 10      continue   
      elseif(dim.eq.3) then
         call lnkerr('not yet implimented')
      endif
      return
      end       



