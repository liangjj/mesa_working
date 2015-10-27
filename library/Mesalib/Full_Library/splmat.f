*deck splmat
      subroutine splmat(n,xdata,k,sc)
      implicit real*8(a-h,o-z)
      dimension sc((4*n+1)*k), xdata(n)
      common /io/ inp, iout
c
c  set up matrix to solve for coefficients in b-spline representation.
c
      ia=n+k
      iq=ia+n
      idummy=iq+n*(3*k-2)
      iscr=idummy+k
      do 30 i=1,n
          call interv(sc(k),n+2-k,xdata(i),ileft,mflag)
          ileft=ileft+k-1
          if (mflag.lt.0) then
              call lnkerr('error in splmat')
          else
              if (mflag.gt.0) then
                  if(i .lt. n) then
                     call lnkerr('error in splmat')
                  endif
                  ileft=n
              endif
          endif
          call bsplvb(sc,k,1,xdata(i),ileft,sc(idummy+1))
          l=ileft-i
          do 16 j=1,k
              l=l+1
              sc(iq+i+n*(l-1))=sc(idummy+j)
   16     continue      
          if(sc(iq+i+n*(k-1)).eq.0.d0) then
             call lnkerr('error in splmat')
          endif
   30 continue      
      return
      end      
