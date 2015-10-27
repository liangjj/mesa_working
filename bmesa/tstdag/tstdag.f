*deck tstdag.f 
c***begin prologue     tstdag
c***date written       951230   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            driver for testing davidson
c
c***references       
c
c***routines called    
c***end prologue       tstdag
      program tstdag
c
      implicit integer (a-z)
      character*4096 ops
      character*8 cpass, itdiag, prtflg
      character*800 card
      character*16 chrkey
      character*128 fillam
      logical posinp
      common z(1)
      dimension ia(1)
      equivalence (ia(1),z(1))
      common/io/inp, iout      
      common/memory/ioff
      real*8 z, thresh, cnverg, fpkey
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)
      itdiag=chrkey(ops,'iterative-diagonalization','false',' ')
      call iosys ('read character "linear algebraic filename" from rwf',
     1                 -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      if ( posinp('$tstdag',cpass) ) then
           call cardin(card)
           n=intkey(card,'matrix-size',2,' ')
      endif
      if(itdiag.eq.'davidson') then
         nroots=intkey(ops,'davidson=number-of-roots',n,' ')
         thresh=fpkey(ops,'davidson=tolerance',1.0d-06,' ')
         cnverg=fpkey(ops,'davidson=convergence',1.d-08,' ')
         maxvec=intkey(ops,'davidson=maximum-number-of-vectors',
     1                 n,' ')
         niter=intkey(ops,'davidson=maximum-number-of-iterations',
     1                n,' ')
         nattim=intkey(ops,'davidson=number-of-roots-at-a-time',
     1                 nroots,' ')
         maxvec=min(maxvec,n)
         niter=min(niter,n)
         write(iout,1) nroots, maxvec, niter, nattim, thresh, cnverg
      endif
      ioff=1
      do 10 i=1,2
         ham=ioff
         eig=ham+n*n
         work=eig+n
         dum=work+n
         vec=dum+n
         hvec=vec+n*maxvec
         root=hvec+n*maxvec
         dvdmat=root+n
         dvdvec=dvdmat+maxvec*maxvec
         cvn=dvdvec+maxvec*maxvec
         hmev=cvn+n
         words=wpadti(hmev+n*maxvec)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
         else
             call getscm(words,z,ngot,'tstdag',0)
         endif
   10 continue
      call mkeham(z(ham),n)
      if(itdiag.eq.'false') then
             call diag(z(ham),z(eig),z(work),z(dum),n)
      elseif(itdiag.eq.'davidson') then
             call david(z(ham),z(eig),z(vec),z(hvec),z(root),
     1                  z(dvdmat),z(dvdvec),z(cvn),z(hmev),thresh,
     2                  cnverg,ops,n,nroots,nattim,maxvec,niter,
     3                  prtflg)
      else
             call lnkerr('error in diagonalization technique')
      endif 
      call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
      call chainx(0)               
      stop
 1    format(/,5x,'davidson iterative diagonalization',/,5x,
     1            'number of roots                          = ',i3,/,5x,
     2            'maximum number of vectors stored in core = ',i3,/,5x,
     3            'maximum number of iterations             = ',i3,/,5x,
     4            'number of roots at a time                = ',i3,/,5x,
     5            'threshold criterion for vectors          = ',e15.8,
     6                                                           /,5x,
     7            'convergence criterion                    = ',e15.8)
      end

