      program drvdif
      implicit integer(a-z)
      logical dollar, logkey, store, schroe
      character*1600 card
      character*4096 ops
      character*80 cpass
      pointer (p0,t(1))
      common/io/inp, iout
      real*8 t, fpkey, t0, tn, xin, xout, eps, h1
      real*8 hmin, dxsav, bcond, rpar, acc, stpsze, tstrt, tfin
      real*8 dum
      dimension eps(2), acc(2), info(15)
      external df
c
      call drum
      call manmem(0,idum,idum,'drvdif',idum)
      call iosys ('read character options from rwf',-1,0,0,ops)
      if ( dollar('$tprop',card,cpass,inp) ) then
         schroe=logkey(card,'real-time',.false.,' ')
         t0=fpkey(card,'initial-t-value',0.d0,' ')
         tn=fpkey(card,'final-t-value',5.d0,' ')
         h1=fpkey(card,'guess-step-size',.001d0,' ')
         h1=fpkey(card,'step-size',h1,' ')
         hmin=fpkey(card,'minimum-step-size',.1d0*h1,' ')
         words=(xlast-xfirst)/h1 + 100
         need=wptoin(words)
         call manmem(need,p0,ngot0,'drvdif',0)
         call grid(t(1),t0,tn,h1,nstp)
         eps(1)=1.d-06
         eps(2)=0.d0
         call fparr(card,'accuracy',eps,2,' ')
         bcond=fpkey(card,'boundary-condition',0.d0,' ')
         store=logkey(card,'store-wavefunction',.false.,' ')
         if(store) then
            kmax=nstp
         else
            kmax=0
         endif
         dxsav=fpkey(card,'storage-interval',h1,' ')
         write(iout,1) t0, tn, nstp, bcond
      endif

      y=1
      yp=y+neq
      rwork=yp+neq*nstp
      rwords=rwork+130+21*neq
      iwords=51
      iwork=wpadti(rwords)
      need=iwork+iwords
      call manmem(need,p1,ngot1,'drvdif',0)
      xin=xfirst
      xout=xlast
      acc(1)=eps(1)
      acc(2)=eps(2)
      call init(z(y),bcond,schroe,neq)
      call copy(z(y),z(yp),neq)
      call deabm(x,z(yp),z(y),xin,xout,h1,info,
     1           acc(1),acc(2),z(rwork),rwords,ia(iwork),
     2           iwords,neq,nstp-1,
     3           dum,dum,dum,dum,cpass,dum,dum,cpass,dum,dum,
     4           idum,idum,idum)     
      call prdiff(x,z(yp),schroe,neq,nstp)  
      call manmem(-ngot1,p1,idum,'drvdif',idum)
 1    format(/,5x,'values of input parameters:',/,5x,
     1            'starting value of t        = ',e15.8,/,5x,
     2            'ending   value of t        = ',e15.8,/,5x,
     3            'number of steps            = ',i5,/,5x,
     4            'initial boundary condition = ',e15.8)
      call manmem(-ngot0,p0,idum,'drvdif',idum)
      call chainx(0)
      stop
      end
