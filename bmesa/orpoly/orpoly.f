*deck orpoly.f 
c***begin prologue     orpol
c***date written       951230   (yymmdd)
c***revision date               (yymmdd)
c***keywords           orthogonal polynomials
c***author             schneider, b. i.(nsf)
c***source             orpoly
c***purpose            driver orthogonal polynomials
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       orpoly
      program orpoly
c
      implicit integer (a-z)
      character*4096 ops
      character*8 cpass
      character*16 chrkey, type, pottyp
      character*800 card
      character*128 fillam
      logical posinp, logkey, prnply, prnpwc, prnh, check, decomp, fixed
      logical diffeq
      common z(1)
      dimension ia(1), endpts(2), der(2)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, left, right, fpkey, endpts, alpha, beta, norm0, der
c
      call drum
      write(iout,*)
      write(iout,*) '          orthogonal polynomials'
      call iosys ('read character options from rwf',-1,0,0,ops)
      prnpwc=logkey(ops,'print=m6240=polynomial-information',
     1              .false.,' ')
      prnply=logkey(ops,'print=m6240=polynomials',.false.,' ')
      prnh=logkey(ops,'print=m6240=hamiltonian',.false.,' ')
      check=logkey(ops,'check-orthogonality',.false.,' ')
      decomp=logkey(ops,'decompose',.false.,' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$poly',cpass) ) then
           call cardin(card)
           type=chrkey(card,'type-polynomials','legendre',' ')
           diffeq=logkey(card,'solve-differential-equation',.false.,' ')
           pottyp=chrkey(card,'potential','none',' ')
           left=fpkey(card,'left-boundary',-1.d0,' ')
           der(1)=fpkey(card,'left-derivative',0.d0,' ')
           der(2)=fpkey(card,'right-derivative',0.d0,' ')
           right=fpkey(card,'right-boundary',1.d0,' ')
           nmax=intkey(card,'order-of-polynomials',10,' ')
           npts=intkey(card,'number-of-points',nmax,' ')
           fixed=logkey(card,'fix-end-points',.false.,' ')
           if (fixed) then
               nfixed=intkey(card,'number-of-fixed-points',2,' ')
               call fparr(card,'end-points',endpts,nfixed,' ')
           endif
           alpha=0.d0
           beta=0.d0
           if (type.eq.'jacobi'.or.type.eq.'laguerre') then
               alpha=fpkey(card,'alpha',0.d0,' ')
               beta=fpkey(card,'beta',0.d0,' ')
           endif
           write(iout,1) type, nmax, left, right
           if (fixed) then
               if(nfixed.eq.1) then
                  write(iout,2) nfixed, endpts(1)
               endif
               if(nfixed.eq.1) then
                  write(iout,3) nfixed, endpts(1), endpts(2)
               endif 
           endif
           if (type.eq.'jacobi') then
               write(iout,4) type, alpha, beta
           endif
           if (type.eq.'laguerre') then
               write(iout,5) type, alpha
           endif
      endif
      dim=max(npts,nmax)
      ioff=1
      do 10 i=1,2
         x=ioff
         wts=x+dim
         a=wts+dim
         b=a+nmax
         scr=b+nmax
         dum=scr+dim*dim
         pn=dum+nmax*dim
         dpn=pn+nmax*dim
         ddpn=dpn+nmax*dim
         f=ddpn+nmax*dim
         df=f+dim
         ddf=df+dim
         words=wpadti(ddf+dim)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'orpoly',0)
         endif
   10 continue
c
c           calculate the a and b coefficients, the points and
c           the weights.
c
      call cpoly(alpha,beta,endpts,left,right,nfixed,z(x),z(wts),z(a),
     1           z(b),norm0,z(scr),nmax,type,prnpwc)
      call gpoly(z(pn),z(dpn),z(ddpn),left,right,z(x),z(wts),z(a),
     1           z(b),norm0,z(scr),nmax,nmax,type,.false.,prnply)
c           calculate the polynomials and derivatives ( 1st and 2nd )
c           check orthogonality on discrete interval
      if(check) then
         call chk(z(pn),z(wts),z(dum),z(scr),nmax)
      endif
      if(decomp) then
         call fcoef(z(f),z(df),z(ddf),z(pn),z(dpn),z(ddpn),z(x),z(wts),
     1              z(scr),z(dum),nmax,nmax,.false.)
         call grid(z(x),npts)
         call gpoly(z(pn),z(dpn),z(ddpn),left,right,z(x),z(wts),z(a),
     1              z(b),norm0,z(scr),nmax,npts,type,.true.,prnply)
         call fcoef(z(f),z(df),z(ddf),z(pn),z(dpn),z(ddpn),z(x),z(wts),
     1              z(scr),z(dum),nmax,npts,.true.)
      endif
      if (diffeq) then
          call diff(z(pn),z(dpn),z(ddpn),z(x),z(wts),z(scr),z(dum),
     1              z(f),der,nmax,nmax,pottyp,prnh)
      endif
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'type polynomial = ',a16,/,5x,
     1            'order of polynomials = ',i3,/,5x,
     2            'left boundary = ',e15.8,/,5x,
     3            'right boundary = ',e15.8)
 2    format(/,5x,'number fixed end points = ',i1,/,5x,     
     1            'end point = ',e15.8)
 3    format(/,5x,'number fixed end points = ',i1,/,5x,     
     1            'left end point = ',e15.8,/,5x,
     2            'right end point = ',e15.8)
 4    format(/,5x,'for polynomial type  = ',a16,/,5x,
     1            'alpha = ',e15.8,/,5x,
     2            'beta = ',e15.8)
 5    format(/,5x,'for polynomial type  = ',a16,/,5x,
     1            'alpha = ',e15.8)
      end

