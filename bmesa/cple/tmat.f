*deck tmat
c***begin prologue     tmat
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           bessel, kinetic-energy
c***author             schneider, barry (nsf)
c***source             
c***purpose            kinetic energy plus bloch operator integrals
c***                   of numerical basis functions
c***description        
c***references       
c
c***routines called
c***end prologue       tmat
      subroutine tmat(fi,dfi,ddfi,fj,dfj,ddfj,flsti,flstj,dflsti,dflstj,
     1                t,ni,nj,npts,n,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 fi, dfi, ddfi, fj, dfj, ddfj, flsti, flstj, dflsti, dflstj
      real*8 t
      character*80 title
      logical prnt
      dimension fi(npts,ni), dfi(npts,ni), ddfi(npts,ni)
      dimension fj(npts,nj), dfj(npts,nj), ddfj(npts,nj)
      dimension flsti(ni), flstj(nj), dflsti(ni), dflstj(nj)
      dimension t(n,n)
      nbeg=ni+1
      nfinal=ni+nj  
c     calculate kinetic energy plus bloch operator integral
      if (ni.ne.0) then
          call ebtcxx(t,fi,ddfi,ni,npts,ni,n,npts,npts)
          do 10 i=1,ni
             do 20 j=1,i
                t(i,j)=-.5d+00*( t(i,j) - flsti(i)*dflsti(j) )
 20          continue
 10       continue   
          if (nj.ne.0) then
              call ebtcxx(t(nbeg,1),fj,ddfi,nj,npts,ni,n,npts,npts)
              jcount=0
              do 30 j=nbeg,nfinal
                 jcount=jcount+1
                 do 40 i=1,ni
                    t(j,i)=-.5d+00*( t(j,i)-flstj(jcount)*dflsti(i) )
 40              continue
 30           continue   
          endif
      endif
      if (nj.ne.0) then
          call ebtcxx(t(nbeg,nbeg),fj,ddfj,nj,npts,nj,n,npts,npts)
          icount=0
          do 50 i=nbeg,nfinal
             icount=icount+1
             jcount=0
             do 60 j=nbeg,i
                jcount=jcount+1
                t(i,j)=-.5d0*( t(i,j) - 
     1                         flstj(icount)*dflstj(jcount) )
 60          continue
 50       continue   
      endif
      do 70 i=1,n
         do 80 j=1,i
            t(j,i)=t(i,j)
   80    continue
   70 continue                                
      if (prnt) then
          title='primitive kinetic energy matrix'
          call prntrm(title,t,n,n,n,n,iout)
      endif
      return
      end
