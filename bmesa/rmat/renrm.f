*deck renrm
c***begin prologue     renrm
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           renormalize
c***author             schneider, barry (nsf)
c***source             
c***purpose            normalize a matrix and its basis
c***references       
c
c***routines called
c***end prologue       renrm
      subroutine renrm(s,f,df,ddf,flst,dflst,tmp,npts,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 s, f, df, ddf, flst, dflst, tmp
      character*80 title
      dimension s(n,n), f(npts,n), df(npts,n), ddf(npts,n), tmp(*)
      dimension flst(n), dflst(n) 
      logical prnt
      do 10 i=1,n
         tmp(i)=1.d0/sqrt(s(i,i))
   10 continue
      do 20 i=1,n
         do 30 j=1,i
            s(i,j)=s(i,j)*tmp(i)*tmp(j)
            s(j,i)=s(i,j)
   30    continue
   20 continue
      do 40 i=1,n
         call sscal(npts,tmp(i),f(1,i),1)
         call sscal(npts,tmp(i),df(1,i),1)
         call sscal(npts,tmp(i),ddf(1,i),1)
         flst(i)=flst(i)*tmp(i)
         dflst(i)=dflst(i)*tmp(i)
   40     continue
      if (prnt) then
          title='re-normalized overlap matrix'
          call prntrm(title,s,n,n,n,n,iout)
      endif
      return
      end
