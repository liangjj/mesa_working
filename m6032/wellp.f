*deck m6032
c***begin prologue     m6032
c***date written       940101   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6032 ,link 6032
c***author             schneider, barry (nsf)
c***source             m6032
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       m6032
      program wellp
      implicit integer (a-z)
      parameter ( nbasis=200, nquad=200 )
      real*8 z, rbox, fpkey, v0, vl
      common a(1)
      dimension z(1)
      equivalence (z,a)
      common /memory/ ioff
      common /io / inp, iout
      character*4096 ops
      character*1600 card
      character*10 cpass
      character*80 title
      logical logkey, pnorm
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      pnorm=logkey(ops,'calculate-normalization',.false.,' ')
      call posinp('$well',cpass)
      call cardin(card)
      lval=intkey(card,'l-value',0,' ')
      rbox=fpkey(card,'r-matrix-box-size',5.d0,' ')
      v0=fpkey(card,'well-depth',-2.d0,' ')
      vl=fpkey(card,'well-length',rbox,' ')
      nprim=intkey(card,'number-of-polynomials',5,' ')
      nq=intkey(card,'number-of-quadrature-points',2*nprim+1,' ')
      write(iout,1) nprim,nq,rbox,v0,vl
      ndim=nprim+1
      nplus=ndim+1
      nq=nq+1
      ioff=1
      do 20 i=1,2
         pt=ioff
         wt=pt+nq
         x=wt+nq
         ply=x+nq
         dply=ply+nplus*nq
         ddply=dply+ndim*nq
         norm=ddply+ndim*nq
         ekin=norm+ndim*ndim
         v=ekin+ndim*ndim
         ham=v+ndim*ndim
         eig=ham+ndim*ndim
         dum=eig+ndim
         words=wpadti(dum+ndim)
         if (i.eq.1) then
             call iosys ('read integer maxsiz from rwf',1,
     1                    canget,0,' ')
             if (words.gt.canget) then
                 call lnkerr('not enough memory. will quit')
             endif
             call iosys ('write integer maxsiz to rwf',1,
     1                    words,0,' ')
             call getscm(words,z,ngot,'m6032',0)
             write (iout,*) 'get ',words,' words of core'      
         endif
 20   continue
c
c
      ic=pt
      jc=wt
      do 30 i=1,nq-1
         call lgndrx (nq-1,i,z(jc),z(ic))
         z(jc)=z(jc)*rbox
         z(ic)=z(ic)*rbox
         ic=ic+1
         jc=jc+1
   30 continue
      z(ic)=rbox
      z(jc)=0.d0
      call polyab (z(ply),z(dply),z(ddply),z(pt),z(x),0.d0,rbox,nq,
     1             nprim,.false.)
      if (pnorm) then
          call tstnrm(z(ply),z(wt),z(norm),nq,ndim)
      endif
      call kinmat(z(ply),z(dply),z(ddply),z(pt),z(wt),z(ekin),rbox,
     1            lval,nprim,nq)      
      call vmat(z(ply),z(pt),z(wt),z(v),v0,vl,nprim,nq)
      call vadd(z(ham),z(ekin),z(v),ndim*ndim)
      title='hamiltonian matrix'
      call prntrm(title,z(ham),ndim,ndim,ndim,ndim,iout)
      call diag(z(ham),z(eig),z(dum),ndim)
      call iosys ('write integer maxsiz to rwf',1,
     1             canget,0,' ')
      call chainx(0)
    1 format(/,2x,'nprim = ',i3,1x,'nq = ',i3,1x,'rbox = ',e15.8,
     1       /,2x,'well depth = ',e15.8,1x,'well length = ',e15.8)
      stop
      end


