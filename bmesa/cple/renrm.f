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
      subroutine renrm(s,fi,dfi,ddfi,fj,dfj,ddfj,flsti,flstj,dflsti,
     1                 dflstj,tmp,ni,nj,npts,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 s, fi, dfi, ddfi, fj, dfj, ddfj
      real*8 flsti, flstj, dflsti, dflstj, tmp
      character*80 title
      dimension s(n,n), fi(npts,ni), dfi(npts,ni), ddfi(npts,ni)
      dimension fj(npts,nj), dfj(npts,nj), ddfj(npts,nj), tmp(*)
      dimension flsti(ni), flstj(nj), dflsti(ni), dflstj(nj) 
      logical prnt
      nbeg=ni+1
      nfinal=ni+nj
      do 10 i=1,nfinal
         tmp(i)=1.d0/sqrt(s(i,i))
   10 continue
      do 20 i=1,nfinal
         do 30 j=1,i
            s(i,j)=s(i,j)*tmp(i)*tmp(j)
            s(j,i)=s(i,j)
   30    continue
   20 continue
      if (ni.ne.0) then
          do 40 i=1,ni
             call sscal(npts,tmp(i),fi(1,i),1)
             call sscal(npts,tmp(i),dfi(1,i),1)
             call sscal(npts,tmp(i),ddfi(1,i),1)
             flsti(i)=flsti(i)*tmp(i)
             dflsti(i)=dflsti(i)*tmp(i)
   40     continue
      endif
      if (nj.ne.0) then
          count=ni
          do 50 i=1,nj
             count=count+1
             call sscal(npts,tmp(count),fj(1,i),1)
             call sscal(npts,tmp(count),dfj(1,i),1)
             call sscal(npts,tmp(count),ddfj(1,i),1)
             flstj(i)=flstj(i)*tmp(count)
             dflstj(i)=dflstj(i)*tmp(count)
   50     continue
      endif                         
      if (prnt) then
          title='re-normalized overlap matrix'
          call prntrm(title,s,nfinal,nfinal,n,n,iout)
      endif
      return
      end
