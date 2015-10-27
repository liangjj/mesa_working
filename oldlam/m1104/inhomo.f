      subroutine inhomo (g1,g2,fold,fnew,pt,wt,ncst,norb,chnloc,jind
     1 ,nsts,ntchn,nptmx,ldim,nmo,ncmax,lplsmx,ops)
      implicit integer(a-z)
      logical logkey
      real *8g1, g2, fold, fnew, pt, wt
      character *(*) ops
      common /io/ inp, iout
      dimension g1(nptmx,ntchn), g2(nptmx,ntchn), fold(nptmx,ldim,nmo)
      dimension fnew(nptmx,ldim,nmo), pt(nptmx), wt(nptmx), norb(nsts)
      dimension ncst(nsts), chnloc(lplsmx,nsts)
      dimension jind(ncmax,nsts)
*
*
*          the minus sign in the definition of the inhomogeneity
*          comes from the greens function definition as the negative
*          of the g1 , g2 product. the 2. multiplying the fnew absorbs
*         the two's on the rhs of the integral equation.
*
      do 10 i=1,nmo
      do 10 j=1,ldim
      do 10 k=1,nptmx
   10 fnew(k,j,i)=0.d+00
      mo=0
      do 50 is=1,nsts
      nmos=norb(is)
      if (nmos.eq.0) go to 50
      nsch=ncst(is)
      do 40 i=1,nsch
      lpls=jind(i,is)
      ncl=chnloc(lpls,is)
      mocnt=mo
      do 30 j=1,nmos
      mocnt=mocnt+1
      call grophi (g1(1,ncl),g2(1,ncl),fold(1,lpls,mocnt),fnew(1,lpls
     1 ,mocnt),pt,wt,nptmx)
      do 20 k=1,nptmx
   20 fnew(k,lpls,mocnt)=-2.d+00*fnew(k,lpls,mocnt)
   30 continue
   40 continue
   50 mo=mo+nmos
      if (logkey(ops,'print=lam=newrhs',.false.,' ')) then
      write (iout,70)
      do 60 i=1,mocnt
      write (iout,80) i
      call matprt (fnew(1,1,i),nptmx,ldim,nptmx,ldim,0,0,0,0,0,0)
   60 continue
      endif
      return
c
   70 format (/,5x,'newrhs')
   80 format (/,5x,'mo',1x,i4)
      end
