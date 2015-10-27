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
     1                   scr,l,m,nr,nth,nph,nang,nonsep,prnt)
      implicit integer (a-z)
      real *8 f, plm, phifn, ylm, flm, wtthe, wtphi, wtang, scr
      logical prnt, nonsep
      character*2 itoc, tmp
      character*80 title
      dimension f(nr,*), plm(nth,*), phifn(nph,*), ylm(nang,*)
      dimension scr(nr,nth,2), flm(nr,*), wtthe(nth,2), wtphi(nph,2)
      dimension wtang(nang,2)
      common /io/ inp, iout
      if (nonsep) then
          cntlm=1
          do 10 mu=0,m
             ndim=l-mu+1 
             tmp=itoc(mu)
             len=cskipb(tmp,' ')
             title='y(l,'//tmp(1:len)//') decomposition'
             no=2
             if ( mu.eq.0 ) then
                  no=1
             endif
             do 20 count=1,no
                call vmmul(wtang,ylm(1,cntlm),ylm(1,cntlm),nang,ndim)
                call ebc(flm(1,cntlm),f,ylm(1,cntlm),nr,nang,ndim)
                call setzro(flm(1,cntlm),nr*ndim)
                if (prnt) then
                    write(iout,1) count,mu
                    call prntfm(title,flm(1,cntlm),nr,ndim,nr,
     1                          ndim,iout)
                endif
c
c      scale functions by inverse of weights to return to original state.
                call vmmul(wtang(1,2),ylm(1,cntlm),ylm(1,cntlm),
     1                     nang,ndim)
                cntlm=cntlm+ndim
   20        continue
   10     continue
      else      
          plmcnt=1
          flmcnt=1
          phicnt=1
          do 30 mu=0,m
             ndim=l-mu+1 
             tmp=itoc(mu)
             len=cskipb(tmp,' ')
             title='p(l,'//tmp(1:len)//') decomposition'
             no=2
             if ( mu.eq.0) then
                  no=1
             endif
c
c           scale functions by weights to make p(l,m) projections
c           easily vectorized.
c 
c                    call scalfn(plm(1,plmcnt),wtthe,nth,l,mu)
c                    call scalfn(phifn(1,phicnt),wtphi,nph,no,1)
             call vmmul(wtthe,plm(1,plmcnt),plm(1,plmcnt),nth,ndim)
             call vmmul(wtphi,phifn(1,phicnt),phifn(1,phicnt),nph,no)
             call ebc(scr,f,phifn(1,phicnt),nr*nth,nph,no)
             do 40 count=1,no
                call ebc(flm(1,flmcnt),scr(1,1,count),plm(1,plmcnt),
     1                   nr,nth,ndim)
                call setzro(flm(1,flmcnt),nr*ndim)
                if (prnt) then
                    write(iout,1) count,mu
                    call prntfm(title,flm(1,flmcnt),nr,ndim,nr,
     1                          ndim,iout)
                endif
                flmcnt=flmcnt+ndim
   40        continue
c
c           scale functions by inverse of weights to return to original state.
c 
c                    call scalfn(plm(1,plmcnt),wtthei,nth,l,mu)
c                    call scalfn(phifn(1,phicnt),wtphii,nph,no,1)
             call vmmul(wtthe(1,2),plm(1,plmcnt),plm(1,plmcnt),nth,ndim)
             call vmmul(wtphi(1,2),phifn(1,phicnt),phifn(1,phicnt),
     1                  nph,no)    
             phicnt=phicnt+no 
             plmcnt=plmcnt+ndim
   30     continue
      endif   
    1 format(/,'the ',i2,' legendre decomposition matrix for m = '
     1        ,i3)
      return
      end
