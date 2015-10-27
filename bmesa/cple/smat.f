*deck smat
c***begin prologue     smat
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, overlap
c***author             schneider, barry (nsf)
c***source             
c***purpose            overlap of numerical basis functions
c***description        
c***references       
c
c***routines called
c***end prologue       smat
      subroutine smat(fi,fj,s,ni,nj,npts,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fi, fj, s
      character*80 title
      logical prnt
      dimension fi(npts,ni), fj(npts,nj), s(n,n)
      nbeg=ni+1
      nfinal=ni+nj  
c     calculate overlap integral
      if (ni.ne.0) then
          call ebtcxx(s,fi,fi,ni,npts,ni,n,npts,npts)
          if (nj.ne.0) then
              call ebtcxx(s(nbeg,1),fj,fi,nj,npts,ni,n,npts,npts)
          endif
      endif
      if (nj.ne.0) then
          call ebtcxx(s(nbeg,nbeg),fj,fj,nj,npts,nj,n,npts,npts)
      endif
      do 10 i=1,n
         do 20 j=1,i
            s(j,i)=s(i,j)
   20    continue
   10 continue                                
      if (prnt) then
          title='primitive overlap matrix'
          call prntrm(title,s,n,n,n,n,iout)
      endif
      return
      end
