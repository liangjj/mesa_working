*deck exact.f
c***begin prologue     exact
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            exact approximation to error.
c***references         
c
c***routines called    
c***end prologue       exact
      subroutine exact(resid,d,ham,root,temp,ipvt,n)
      implicit integer (a-z)
      real*8 resid, d, ham, root, temp, nrzero
      dimension resid(n), d(n), ham(n,n), temp(n,n), ipvt(n)
      common/io/inp, iout
      data nrzero / 1.0d-06 /
      call rzero(temp,n*n)
      do 10 i=1,n
         temp(i,i) = root-d(i)
         do 20 j=1,n
            temp(i,j) = - ham(i,j) + temp(i,j)
 20      continue
 10   continue
      call sgefa(temp,n,n,ipvt,info)
      if(info.ne.0) then
         call lnkerr('error in matrix solve')
      endif   
      call sgesl(temp,n,n,ipvt,resid,0)
      return
      end       
