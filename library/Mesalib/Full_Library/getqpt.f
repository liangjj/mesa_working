*deck getqpt
      subroutine getqpt(x,wt,lftept,rtept,qtyp,trns,
     1                  scr,nfix,nptsr,n,nreg,prnt)
c***begin prologue     getqpt
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            lagrange polynomials.
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       getqpt
c
      implicit integer (a-z)
      real*8 x, wt, lftept, rtept, scr, temp, sumwt, prefac
      real*8 dum
      character*(*) qtyp, trns
      character*80 title
      character*3 itoc
      logical prnt, nfix
      dimension x(n), wt(n), temp(2), scr(n)
      dimension nptsr(nreg), lftept(*), rtept(*), nfix(2,*)
      common /io/ inp, iout
      if(qtyp.eq.'hermite'.and.nreg.ne.1) then
         call lnkerr('hermite quadrature with nreg gt 1')
      endif
c
c     loop over the regions
c
      cntpt=0
      do 10 nr=1,nreg      
         nn=nptsr(nr)
c
c     get the points on the interval ( -1 , +1 )
c
         cnt=0
         if(nfix(1,nr)) then
            cnt=cnt+1
            temp(cnt) = -1.d0
         endif
         if(nfix(2,nr)) then
            cnt=cnt+1
            temp(cnt) = 1.d0
         endif
         if(qtyp.eq.'jacobi') then
            call gaussp(qtyp,nn,0.d0,1.d0,cnt,
     1                  temp,scr,x(cntpt+1),
     2                  wt(cntpt+1),dum,.false.)
         elseif(qtyp.eq.'legendre') then
            call gaussp(qtyp,nn,0.d0,0.d0,cnt,
     1                  temp,scr,x(cntpt+1),
     2                  wt(cntpt+1),dum,.false.)
         elseif(qtyp.eq.'simpson') then
            call gaussp(qtyp,nn,0.d0,0.d0,cnt,
     1                  temp,scr,x(cntpt+1),
     2                  wt(cntpt+1),dum,.false.)
         elseif(qtyp.eq.'hermite') then
            call gaussp(qtyp,nn,0.d0,0.d0,0,temp,scr,x,
     1                  wt,dum,.false.)
         elseif(qtyp.eq.'chebyshev-1') then
            call gaussp(qtyp,nn,0.d0,0.d0,cnt,
     1                  temp,scr,x(cntpt+1),
     2                  wt(cntpt+1),dum,.false.)
         elseif(qtyp.eq.'chebyshev-2') then
            call gaussp(qtyp,nn,0.d0,0.d0,cnt,
     1                  temp,scr,x(cntpt+1),
     2                  wt(cntpt+1),dum,.false.)
         else
            call lnkerr('error in quadrature type')
         endif
         if(prnt) then
            title='primitive gaussian points. quadrature type = '
     1             //qtyp//' region = '//itoc(nr)
            call prntfm(title,x(cntpt+1),nn,1,nn,1,iout)
            title='primitive gaussian weights. quadrature type = '
     1             //qtyp//' region = '//itoc(nr)
            call prntfm(title,wt(cntpt+1),nn,1,nn,1,iout)
         endif
         sumwt=0.d0
         do 30 i=1,nn
            sumwt = sumwt + wt(i+cntpt)
 30      continue
         if(prnt) then
            write(iout,1) nr, sumwt         
         endif
         if(qtyp.eq.'hermite') then
            return
         endif
         if(trns.eq.'before') then
            call chnvar(x(cntpt+1),wt(cntpt+1),-1.d0,1.d0,
     1                  lftept(nr),rtept(nr),prefac,nn)
         endif
         cntpt=cntpt+nn
 10   continue
      if(prnt) then
         title='transformed gaussian points. quadrature type = '
     1             //qtyp
         call prntfm(title,x,n,1,n,1,iout)
         title='transformed gaussian weights. quadrature type = '
     1             //qtyp
         call prntfm(title,wt,n,1,n,1,iout)
      endif
      sumwt=0.d0
      do 40 i=1,n
         sumwt=sumwt+wt(i)
 40   continue   
      if(prnt) then
         write(iout,2) sumwt         
      endif
      return
 1    format(/,5x,'region = ',i3,1x,'sum of the weights = ',e15.8)      
 2    format(/,5x,'full sum of the weights = ',e15.8)      
      end















