*deck  @(#)pm819.f	1.5 7/30/91
      program m819
      implicit integer(a-z)
      real*8 z
      integer a
      pointer(p,z(1)),(p,a(1))
      character*4096 ops
      character*8 errtyp
      common /io/ inp,iout
c
      call drum
      write(iout,99)
  99  format(' m819: zeroing integrals ',/)
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nbf,0,' ')
      nnp=nbf*(nbf+1)/2
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     check on existence of nsmall
c
      nsmall=intkey(ops,'scattering=nsmall',0,' ')
      if(nsmall.eq.0) then
         nsmall=intkey(ops,'internal-orbitals',0,' ')
      endif
      if(nsmall.eq.0) then
	 nsmall=intkey(ops,'m819=nsmall',0,' ')
      endif
      if (nsmall.eq.0 ) then
         write(iout,*) 'No flag set for number of internal '
         write(iout,*) 'orbitals. Must quit'
         call lnkerr(' m819: input error')
      else
          write(iout,1) nsmall	 
      end if
    1 format(/,5x,'number of internal orbitals = ',i3)      
      call iosys('write integer "internal orbitals" to rwf',
     1            1,nsmall,0,' ')
      norbp=nbf-nsmall
      call iosys('write integer "p space orbitals" to rwf',
     1            1,norbp,0,' ')
c
      call getmem(0,p,ngot,'first',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
      s=1
      v=s+nnp
      t=v+nnp
      values=t+nnp
      left2=iadtwp(maxcor)-values
      ntriang=left2/nnp
      if(ntriang.lt.1) then
         write(iout,*)' maxcor left2 nnp values ',
     $                  iadtwp(maxcor),left2,nnp,values
         call lnkerr(' m819: insufficient core ')
      end if
c
      ntriang=min(nnp,ntriang)
      need = wptoin(values + ntriang*nnp)
      call getmem(need,p,ngot,'m819',0)
      
c
      call fixint(z(s),z(t),z(v),z(values),nnp,nbf,nsmall,ntriang,
     $            ops)
c
      call getmem(-ngot,p,idum,'m819',idum)
      call chainx(0)
c
c
      stop
      end
