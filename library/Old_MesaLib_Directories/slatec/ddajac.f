*deck ddajac
      subroutine ddajac (neq, x, y, yprime, delta, cj, h, ier, wt, e,
     *   wm, iwm, res, ires, uround, jac, rpar, ipar, ntemp)
c***begin prologue  ddajac
c***subsidiary
c***purpose  compute the iteration matrix for ddassl and form the
c            lu-decomposition.
c***library   slatec (dassl)
c***type      double precision (sdajac-s, ddajac-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------------
c     this routine computes the iteration matrix
c     pd=dg/dy+cj*dg/dyprime (where g(x,y,yprime)=0).
c     here pd is computed by the user-supplied
c     routine jac if iwm(mtype) is 1 or 4, and
c     it is computed by numerical finite differencing
c     if iwm(mtype)is 2 or 5
c     the parameters have the following meanings.
c     y        = array containing predicted values
c     yprime   = array containing predicted derivatives
c     delta    = residual evaluated at (x,y,yprime)
c                (used only if iwm(mtype)=2 or 5)
c     cj       = scalar parameter defining iteration matrix
c     h        = current stepsize in integration
c     ier      = variable which is .ne. 0
c                if iteration matrix is singular,
c                and 0 otherwise.
c     wt       = vector of weights for computing norms
c     e        = work space (temporary) of length neq
c     wm       = real work space for matrices. on
c                output it contains the lu decomposition
c                of the iteration matrix.
c     iwm      = integer work space containing
c                matrix information
c     res      = name of the external user-supplied routine
c                to evaluate the residual function g(x,y,yprime)
c     ires     = flag which is equal to zero if no illegal values
c                in res, and less than zero otherwise.  (if ires
c                is less than zero, the matrix was not completed)
c                in this case (if ires .lt. 0), then ier = 0.
c     uround   = the unit roundoff error of the machine being used.
c     jac      = name of the external user-supplied routine
c                to evaluate the iteration matrix (this routine
c                is only used if iwm(mtype) is 1 or 4)
c-----------------------------------------------------------------------
c***routines called  dgbfa, dgefa
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901010  modified three max calls to be all on one line.  (fnf)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c   901101  corrected purpose.  (fnf)
c***end prologue  ddajac
c
      integer  neq, ier, iwm(*), ires, ipar(*), ntemp
      double precision
     *   x, y(*), yprime(*), delta(*), cj, h, wt(*), e(*), wm(*),
     *   uround, rpar(*)
      external  res, jac
c
      external  dgbfa, dgefa
c
      integer  i, i1, i2, ii, ipsave, isave, j, k, l, lenpd, lipvt,
     *   lml, lmtype, lmu, mba, mband, meb1, meband, msave, mtype, n,
     *   npd, npdm1, nrow
      double precision  del, delinv, squr, ypsave, ysave
c
      parameter (npd=1)
      parameter (lml=1)
      parameter (lmu=2)
      parameter (lmtype=4)
      parameter (lipvt=21)
c
c***first executable statement  ddajac
      ier = 0
      npdm1=npd-1
      mtype=iwm(lmtype)
      go to (100,200,300,400,500),mtype
c
c
c     dense user-supplied matrix
100   lenpd=neq*neq
      do 110 i=1,lenpd
110      wm(npdm1+i)=0.0d0
      call jac(x,y,yprime,wm(npd),cj,rpar,ipar)
      go to 230
c
c
c     dense finite-difference-generated matrix
200   ires=0
      nrow=npdm1
      squr = sqrt(uround)
      do 210 i=1,neq
         del=squr*max(abs(y(i)),abs(h*yprime(i)),abs(wt(i)))
         del=sign(del,h*yprime(i))
         del=(y(i)+del)-y(i)
         ysave=y(i)
         ypsave=yprime(i)
         y(i)=y(i)+del
         yprime(i)=yprime(i)+cj*del
         call res(x,y,yprime,e,ires,rpar,ipar)
         if (ires .lt. 0) return
         delinv=1.0d0/del
         do 220 l=1,neq
220      wm(nrow+l)=(e(l)-delta(l))*delinv
      nrow=nrow+neq
      y(i)=ysave
      yprime(i)=ypsave
210   continue
c
c
c     do dense-matrix lu decomposition on pd
230      call dgefa(wm(npd),neq,neq,iwm(lipvt),ier)
      return
c
c
c     dummy section for iwm(mtype)=3
300   return
c
c
c     banded user-supplied matrix
400   lenpd=(2*iwm(lml)+iwm(lmu)+1)*neq
      do 410 i=1,lenpd
410      wm(npdm1+i)=0.0d0
      call jac(x,y,yprime,wm(npd),cj,rpar,ipar)
      meband=2*iwm(lml)+iwm(lmu)+1
      go to 550
c
c
c     banded finite-difference-generated matrix
500   mband=iwm(lml)+iwm(lmu)+1
      mba=min(mband,neq)
      meband=mband+iwm(lml)
      meb1=meband-1
      msave=(neq/mband)+1
      isave=ntemp-1
      ipsave=isave+msave
      ires=0
      squr=sqrt(uround)
      do 540 j=1,mba
         do 510 n=j,neq,mband
          k= (n-j)/mband + 1
          wm(isave+k)=y(n)
          wm(ipsave+k)=yprime(n)
          del=squr*max(abs(y(n)),abs(h*yprime(n)),abs(wt(n)))
          del=sign(del,h*yprime(n))
          del=(y(n)+del)-y(n)
          y(n)=y(n)+del
510       yprime(n)=yprime(n)+cj*del
      call res(x,y,yprime,e,ires,rpar,ipar)
      if (ires .lt. 0) return
      do 530 n=j,neq,mband
          k= (n-j)/mband + 1
          y(n)=wm(isave+k)
          yprime(n)=wm(ipsave+k)
          del=squr*max(abs(y(n)),abs(h*yprime(n)),abs(wt(n)))
          del=sign(del,h*yprime(n))
          del=(y(n)+del)-y(n)
          delinv=1.0d0/del
          i1=max(1,(n-iwm(lmu)))
          i2=min(neq,(n+iwm(lml)))
          ii=n*meb1-iwm(lml)+npdm1
          do 520 i=i1,i2
520         wm(ii+i)=(e(i)-delta(i))*delinv
530      continue
540   continue
c
c
c     do lu decomposition of banded pd
550   call dgbfa(wm(npd),meband,neq,
     *    iwm(lml),iwm(lmu),iwm(lipvt),ier)
      return
c------end of subroutine ddajac------
      end
