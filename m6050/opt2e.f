*deck m6050
c***begin prologue     m6050
c***date written       920403   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m6050, link 6050, optical potential
c***author             Schneider, Barry (NSF)
c***source             m6050
c***purpose            driver for two electron optical potential
c***description        data on the eigenvalues and v(r) functions
c***                                                lamda
c***                   are read in, the functions constructed on 
c***                   the grid and all the information placed on
c***                   a file for later matrix element construction.
c***                   
c***                   
c*** 
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       m6050
      program opt2e
      implicit integer (a-z)
      logical logkey, grdtyp, logky, uvec, model
      character *4096 ops
      character *1600 card
      character *13 grdnam
      character *8 chrkey, spin
      character *80 title
      character *20 cpass
      character *128 filkne, filgrd, filopt
      character *3 itoc, ans
      real *8 z, pi, charge, fpkey
c----------------------------------------------------------------------c
c                unicos memory management                              c
      common a(1)
      dimension z(1)
      common /memory / ioff 
      equivalence (z,a)
c----------------------------------------------------------------------c
      dimension logky(10)
      common /io/ inp,iout
      data pi /3.14159265358979323846d+00/
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      model=logkey(ops,'m6050=non-s-wave',.false.,' ')
      call posinp ('$opt2e',cpass) 
      call cardin (card)
      angmom=intkey(card,'l-value',0,' ')
      nobfn=intkey(card,'no-basis-functions-in-optical-potential',1,' ')
      nolam=intkey(card,'no-terms-in-optical-potential',1,' ')
      spin=chrkey(card,'spin','singlet',' ')
      charge=fpkey(card,'charge',1.d0,' ')
      uvec=logkey(card,'unit-matrix',.false.,' ')
      call iosys ('read character "kohn data filename" from rwf',-1,0,0,
     1             filkne)
      call iosys ('read character "grid filename" from rwf',-1,0,0,
     1             filgrd)
      call iosys ('read character "optical potential filename" '//
     1            'from rwf',-1,0,0,filopt)
c----------------------------------------------------------------------c
c                    set options                                       c
c----------------------------------------------------------------------c 
      logky(1)=logkey(ops,'print=m6050=eigenvalues',.false.,' ')
      logky(2)=logkey(ops,'print=m6050=eigenvectors',.false.,' ')
      logky(3)=logkey(ops,'print=m6050=basis',.false.,' ')
      logky(4)=logkey(ops,'m6050=check-data-only',.false.,' ')
      logky(5)=logkey(ops,'print=m6050=vlamdas',.false.,' ')
      logky(6)=logkey(ops,'print=m6050=grid',.false.,' ')
      logky(7)=logkey(ops,'print=m6050=vlm',.false.,' ')
      if (logky(4)) then
          logky(1)=.true.
          logky(2)=.true.          
          logky(3)=.true.
      endif
      write (iout,1020)
      write (iout,1030) nobfn, nolam, spin
      grdtyp=logkey(card,'untransformed-grid',.false.,' ')
      grdnam='"trns grid"'
      if (grdtyp) then
          grdnam='"untrns grid"'
      endif
      if (.not.logky(4)) then
c----------------------------------------------------------------------c
c               open grid file and get no. of points                   c
c               need four arrays of dimension pntbuf                   c
c----------------------------------------------------------------------c
          call iosys ('open grid as old',0,0,0,filgrd)
          call iosys ('read integer "no. grid pts" from grid',1,npnts,0,
     1                 ' ')
          call iosys ('read integer "point buffer" from grid',1,
     1                 pntbuf,0,' ')
          call iosys ('open kohndt as unknown',0,0,0,filkne)
          call iosys ('open optint as new on ssd',262144,0,0,filopt)
      else
c----------------------------------------------------------------------c
c                just checking data so pntbuf not needed               c
c                set it to a small number. its not used.               c
c----------------------------------------------------------------------c
          pntbuf=500
      endif
c----------------------------------------------------------------------c
c                 calculate memory requirements                        c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,canget,0,' ')
      wreal=4*nobfn+2*nolam+nolam*2*nobfn+4*pntbuf+
     1     12*nolam*pntbuf+4*pntbuf+110
      wint=iadtwp(nobfn+nobfn)
      words=wreal+wint
      if (words.gt.canget) then
          call lnkerr('insufficient core')
      endif
      call iosys ('write integer maxsiz to rwf',1,words,0,' ')
      call getscm(z,words,ngot,'m6050',0)
      eigval=ioff
      eigvec=eigval+2*nolam
      gam=eigvec+nolam*2*nobfn
      del=gam+2*nobfn
      grid=del+2*nobfn
      vlamda=grid+4*pntbuf
      wlamda=vlamda+nolam*10*pntbuf
      scr=wlamda+nolam*2*pntbuf
      fact=scr+4*pntbuf
      xbuf=fact+101+101
      leg=xbuf+pntbuf
      lval=wpadti(leg+pntbuf)
      mval=lval+nobfn
      call posinp('$basis',cpass)
      call cardin(card)
      call intarr(card,'l-values',a(lval),nobfn,' ')
      call intarr(card,'m-values',a(mval),nobfn,' ')
      call fparr(card,'gamma-values',z(gam),2*nobfn,' ')
      call fparr(card,'delta-values',z(del),2*nobfn,' ')
      if (.not.uvec) then
          l1=eigval
          l2=eigvec
          call posinp('$eigen',cpass)
          read(inp,1000) title
          read(inp,1000) title
          read(inp,1000) title
          do 20 i=1,nolam
             read(inp,*) (z(jj),jj=l1,l1+1)
             l1=l1+2
             read(inp,*) ii
             read(inp,*) (z(jj),jj=l2,l2+2*nobfn-1)
             l2=l2+2*nobfn
   20     continue
      else
          call putone(z(eigval),z(eigvec),nobfn,nolam)
      endif
      if (logky(3)) then
          call prntbs(a(lval),a(mval),z(gam),z(del),nobfn)
      endif
c----------------------------------------------------------------------c
c             convert energies to Hartrees                             c
c----------------------------------------------------------------------c
      call sscal(2*nolam,.5d0,z(eigval),1)
      if (logky(1)) then
          title='complex eigenvalues of optical potential'
          call wrcvec(title,z(eigval),nolam,iout)
          title='complex eigenvectors of optical potential'
          call prntcmn(title,z(eigvec),nobfn,nolam,nobfn,nolam,
     1		       iout,'e')
      endif
      if (.not.logky(4)) then
c----------------------------------------------------------------------c
c                read in charge of first center                        c
c----------------------------------------------------------------------c
          call iosys ('does "nuclear charges" exist on kohndt',0,0,
     1                 0,ans)
          if (ans.ne.'no') then
              call iosys ('read real "nuclear charges" from kohndt',1,
     1                     charge,0,' ')
          endif
          write(iout,8000) charge
          call factl(z(fact),100)
          if (angmom.ne.0) then
              call dfactl(z(fact+101),100)
          endif 
          call iosys ('write integer "no. vlamdas" to optint',1,
     1                 nolam,0,' ')
          call iosys ('write real "complex eigenvalues" to optint',
     1                 2*nolam,z(eigval),0,' ')
c----------------------------------------------------------------------c
c                 create file to hold vlamdas                          c
c----------------------------------------------------------------------c
          filsiz=2*nolam*npnts
          call iosys ('does "complex vlamdas" exist on optint',-1,0,
     1                 0,ans) 
          if (ans.eq.'no') then
              call iosys ('create real "complex vlamdas" on optint',
     1                     filsiz,0,0,' ')
          endif
c----------------------------------------------------------------------c
c           calculate number of passes of grid file needed             c
c----------------------------------------------------------------------c
          nwds=0
          nwrite=0
          npass=npnts/pntbuf
          nolst=npnts-npass*pntbuf
          if (nolst.ne.0) then
              npass=npass+1
          else
              nolst=pntbuf
          endif
          write(iout,3000) pntbuf,npass,words
          do 200 ipass=1,npass
             noptrg=pntbuf
             if (ipass.eq.npass) then
                 noptrg=nolst
             endif
c----------------------------------------------------------------------c
c           read in this block of grid points                          c
c----------------------------------------------------------------------c
             call iosys ('read real '//grdnam//' from grid without '//
     1                   'rewinding',4*noptrg,z(grid),0,' ')   
             if (logky(6)) then
                 call grdprn(z(grid),noptrg,ipass)
             endif
c----------------------------------------------------------------------c
c               zero the potential array for this block                c
c----------------------------------------------------------------------c
             call rzero(z(vlamda),noptrg*nolam*10)
             call rzero(z(wlamda),noptrg*nolam*2)
c----------------------------------------------------------------------c
c                 calculate the four terms in the vlamda               c
c                            expression                                c
c----------------------------------------------------------------------c
             if (model) then
c
c                s-wave calculation
c
                 do 300 ncall=1,4         
                    call svlam(z(vlamda),z(grid),z(eigvec),z(scr),
     1                        z(gam),z(del),charge,z(fact),a(lval),
     2                        a(mval),nobfn,nolam,noptrg,ncall,spin,
     3                        logky(5),logky(7))
  300            continue
                 call swlam(z(wlamda),z(grid),z(eigvec),z(scr),z(gam),
     1                      z(del),charge,z(fact),a(lval),a(mval),nobfn,
     2                      nolam,noptrg,spin,logky(5),logky(7))
                 call sumvw(z(vlamda),z(wlamda),noptrg,nolam,logky(5))
                 call mkwlm(z(vlamda),z(grid),noptrg,nolam)
                 call wrtvlm(z(vlamda),noptrg,nolam,nwrite,nwds,ipass)
             else
c
c            non-s wave calculation
c
                 call mkthet(z(grid),z(xbuf),noptrg) 
                 call legend(z(leg),z(xbuf),z(fact),z(fact+101),noptrg,
     1                       angmom,0,100,'one')
                 call nsvlam(z(vlamda),z(leg),z(grid),z(eigvec),z(scr),
     1                       z(gam),z(del),z(fact),a(lval),a(mval),
     2                       nobfn,angmom,nolam,noptrg,ncall,spin,
     3                       logky(5),logky(7))
             endif
  200     continue
          write (iout,5000) nwrite
          write (iout,6000) nwds
          if (logky(5)) then
              call prnvlm(z(vlamda),nolam,npass,pntbuf,nolst)
          endif 
          call iosys ('rewind all on grid read-and-write',0,0,0,' ')
          call iosys ('rewind all on optint read-and-write',0,0,0,' ')
          call iosys ('rewind all on kohndt read-and-write',0,0,0,' ')
          call iosys ('close grid',0,0,0,' ')
          call iosys ('close optint',0,0,0,' ')
          call iosys ('close kohndt',0,0,0,' ')
          call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
          call chainx(0)
      else
          call iosys ('write integer maxsiz to rwf',1,canget,0,' ')
          call chainx(0)
      endif
      stop
 1000 format(a80)
 1020 format(//,20x,'*****  m6050:complex optical potential code *****')
 1030 format(/,5x,'no. basis functions',1x,i3,1x,'no. optical potential
     1terms',1x,i3,/,20x,'spin',1x,a8)
 3000 format (//,5x,'point buffer',2x,i6,2x,'no. passes of grid file',
     1        1x,i3,1x,'no. words memory required',1x,i8)
 5000 format(/,5x,'no. disk writes',1x,i8)
 6000 format(/,5x,'no. words written',1x,i8)
 8000 format(/,5x,'nuclear charge',1x,f12.6)
      end
