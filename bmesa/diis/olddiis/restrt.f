*deck restrt.f
c***begin prologue     restrt
c***date written       961209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           restrt, pulay, extrapolation
c***author             schneider, barry (nsf)
c***source             math
c***purpose            restart the diis procedure keeping the latest
c***                   updates.
c***                   
c***                   
c***references         peter pulay
c
c***routines called    
c***end prologue       restrt
      subroutine restrt(pvec,dpvec,b,btmp,nobj,nerr,iter,keep,maxit)
      implicit integer (a-z)
      real*8 pvec, dpvec, b, btmp
      character*80 title
      dimension pvec(nobj,maxit+1), dpvec(nerr,maxit+1)
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1)
      common/io/inp, iout
      title='b matrix in'
      call prntrm(title,b,iter,iter,maxit+1,maxit+1,iout)
      write(iout,1) iter, keep
      delta = iter - keep
      if(delta.le.0) then
         return
      else
         count=0
         do 10 i=delta+1, iter+1
            count=count+1
            call copy(pvec(1,i),pvec(1,count),nobj)
   10    continue
         count = 0
         do 20 i=delta+1,iter
            count = count + 1
            call copy(dpvec(1,i),dpvec(1,count),nerr)
   20    continue
         icount = 0
         do 30 i=delta+1, iter
            icount = icount + 1
            jcount = 0
            do 40 j=delta+1, iter
               jcount = jcount +1
               btmp(icount,jcount) = b(i,j)
   40       continue
   30    continue
         do 50 i=1, keep
            do 60 j=1,i
               b(i,j) = btmp(i,j)
               b(j,i)=b(i,j)
   60       continue
   50    continue                
      endif
      title='b matrix out'
      call prntrm(title,b,keep,keep,maxit+1,maxit+1,iout)         
      return
 1    format(/,5x,'restarting diis at iteration = ',i3,
     2            ' with nvec = ',i3)
      end       
