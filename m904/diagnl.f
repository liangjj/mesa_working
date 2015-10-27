*deck @(#)diagnl.f	5.1  11/6/94
      subroutine diagnl(vector,e,row,thresh,coeff,idomnt,p,q,u,y,almax)
c
c  set up initial guesses for ci vectors and other items for
c  diagonalization
c
      implicit real*8 (a-h,o-z)
      character*80 frmt
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
      common /args/ omega,amumu,bi,n,k,m,l,itmax,it,noconv,mu,idgwrt
c
      dimension vector(*), e(*), row(*), thresh(*), coeff(*),
     1  idomnt(*), p(*), q(*), u(*), y(*), almax(*)
c
      data thres/1.0d-07/
c
      write (iw,1000)
c
c  zero the vector array
c
      nlast = nroots*nsef
      do 25 i=1,nlast
   25 vector(i) = 0.0d0
c
c  read in the maximum number of iterations to be allowed
c
      read (5,1201) itmax,idgwrt,omega
 1201 format (2i5,f10.0)
      if( (itmax.le.0).or.(itmax.gt.maxit) ) itmax = maxit
      if( omega.eq.0.0d0 ) omega = 1.4d0
      write (iw,1300) itmax,idgwrt,omega
c
c  read in the initial guesses for the eigenvectors
c
      write (iw,1350)
c
      do 90 i=1,nroots
c
      read (5,1200) nread, thresh(i)
      if(thresh(i).eq.0.0d0) thresh(i)=thres
c
      if( nread-1 ) 30,40,80
c
c  default read option.  coefficient of i th element in the
c  vector set equal to one.
c
   30 limit = 1
      idomnt(1) = i
      coeff(1) = 1.0d0
      go to 50
c
c  read in the coefficients of the dominant configurations.
c
   40 read (5,1200) limit
      read (5,1200) (idomnt(j), coeff(j), j=1,limit)
c
   50 write (iw,1400) i, thresh(i), (idomnt(j), coeff(j), j=1,limit)
c
      do 60 j=1,limit
   60 vector((idomnt(j)-1)*nroots+i) = coeff(j)
      go to 90
c
   80 read (5,1200) msef
      if( (msef.le.0).or.(msef.gt.nsef) ) msef = nsef
      read (5,'(a80)') frmt
      klast = (msef - 1)*nroots + i
      read (5,frmt) (vector(k), k=i,klast,nroots)
      write (iw,1600) i, thresh(i), (vector(k), k=i,klast,nroots)
c
   90 continue
c
      if(idgwrt.gt.0) write (iw,1700)
c
      n=nsef
      k=nroots
      m=nsef
      l=nroots
      call simeig(e,vector,thresh,almax,u,p,q,y,row)
      if(noconv.eq.0) write (iw,1800) it
c
      return
c
 1000 format('1',4x,'ci diagonalization information'//5x,
     2'shavitt diagonalization method (interaction matrix version)')
 1200 format(5(i5,f10.0))
 1300 format (//5x,'maximum number of iterations =',i5//
     1          5x,'print parameter              =',i5//
     2          5x,'over-relaxation factor       =',f7.1)
 1350 format(//5x,36hinitial guesses for the eigenvectors)
 1400 format(/10x,6hroot =,i3,10x,13hconvergence =,1p,d10.1//
     2       (10x, 4(i9,0p,f12.6)) )
 1600 format(/10x,6hroot =,i3,10x,13hconvergence =,1p,d10.1//
     2       (15x,0p,5f15.7) )
 1700 format(//5x,32hdiagonalization of the ci matrix)
 1800 format(//5x,'convergence reached after',i4,' iterations')
c
      end
