*deck lmdcmp.f
c***begin prologue     lmdcmp
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           lmdcmp,y(l,m), projections
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            decompose a function into its legendre projections
c***description        a three dimensional function is decomposed into
c***                   its projection onto p(l,m) phi(m) angular
c***                   functions. the angular functions in the phi
c***                   variable are the normalized sin(m*phi) and
c***                   cos(m*phi). it is assumed that the function is
c***                   stored with no blank spaces in any column direction.
c***references         none
c
c***routines called
c***end prologue       lmdcmp
      subroutine lmdcmp (f,plm,phifn,ylm,flm,wtthe,wtphi,wtang,
     1                   scr,l,m,nr,nthet,nphi,nang,nonsep,prnt)
      implicit integer (a-z)
      real*8 f, plm, phifn, ylm, flm, wtthe, wtphi, wtang, scr
      logical prnt, nonsep
      character*2 itoc, tmp
      character*80 title
      dimension f(*), plm(nthet,*), phifn(nphi,*), ylm(nang,*)
      dimension scr(nr,nthet,2), flm(nr,*)
      dimension wtthe(nthet,2), wtphi(nphi,2), wtang(nang,2)
      common /io/ inp, iout
      if (nonsep) then
          cntlm=1
          do 10 mu=0,m
             dim=l-mu+1 
             tmp=itoc(mu)
             len=length(tmp)
             title='y(l,'//tmp(1:len)//') decomposition'
             no=2
             if ( mu.eq.0 ) then
                  no=1
             endif
             call toylm(ylm(1,cntlm),flm(1,cntlm),f,wtang,l,mu,
     1                  no,nang,nr,title,prnt)
             cntlm=cntlm+no*dim
   10     continue
      else      
          plmcnt=1
          flmcnt=1
          phicnt=1
          do 20 mu=0,m
             dim=l-mu+1 
             tmp=itoc(mu)
             len=length(tmp)
             title='p(l,'//tmp(1:len)//') decomposition'
             no=2
             if ( mu.eq.0) then
                  no=1
             endif
             call toplm(plm(1,plmcnt),phifn(1,phicnt),flm(1,flmcnt),f,
     1                  wtthe,wtphi,scr,l,mu,no,nthet,nphi,nr,
     2                  title,prnt) 
             phicnt=phicnt+no 
             plmcnt=plmcnt+dim
             flmcnt=flmcnt+no*dim
   20     continue
      endif   
      return
      end
