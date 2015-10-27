      program tchank 
      implicit integer(a-z)
      parameter (dime=100, mxchn=50, mxpchn=10)
      real *8 z, zdum, echan, fpkey, egrnd
      real *8 energy, eau, ek
      character *4096 ops
      character *8 cpass
      character *128 filtmt, blktmt, filkne
      character *1600 card
      character *16 chrkey
      character *24 filnm
      character *3 ans, itoc
      character *16 fptoc
      character *80 title
      character *11 fil(mxpchn*mxpchn)
c----------------------------------------------------------------------c
c                unicos memory management                              c
      common a(1)
      dimension z(1)
      common /memory / ioff 
      equivalence (z,a)
c----------------------------------------------------------------------c
      common /io/ inp, iout
      dimension energy(dime), echan(mxpchn), nlm(mxpchn)
      dimension lch(mxchn,mxpchn), mch(mxchn,mxpchn), ek(2)
      equivalence (z,a)
      call drum
      wds=1+4*mxchn*mxchn
      call iosys('write integer maxsiz to rwf',1,wds,0,' ')
      call getscm (wds,z(1),ngot,'tmat',0)
      call iosys ('read character "kohn data filename" from rwf',-1,
     1             0,0,filkne)
      call iosys('open kohndt as old',0,0,0,filkne)
      call iosys('read integer "no. channels" from kohndt',1,
     1            nchan,0,' ')
      call iosys('read real "chan energies" from kohndt',nchan,echan,0,
     1   	 ' ')
      egrnd=echan(1)
      call iosys('read integer "no. energies" from kohndt',1,nen,0,' ')
      call iosys ('read real "scatt energies" from kohndt',nen,
     1		   energy,0,' ')
      call iosys('rewind all on kohndt read-and-write',0,0,0,' ')
      call iosys('close kohndt',0,0,0,' ')
      call iosys('read character "kohn tmatrix filename" from rwf',-1,
     1            0,0,filtmt)
      call iosys('open tmat as old',0,0,0,filtmt)
      call iosys('read character "kohn full tmatrix filename" from '//
     1		 'rwf',-1,0,0,blktmt)
      call iosys('open blkt as new',262144,0,0,blktmt)
c----------------------------------------------------------------------c
c                  position input file                                 c 
c----------------------------------------------------------------------c
      write (iout,2)
      call posinp ('$readkn',cpass)
c----------------------------------------------------------------------c
c                  recover options string                              c
c----------------------------------------------------------------------c
      call iosys ('read character options from rwf',-1,0,0,ops)
c----------------------------------------------------------------------c
c            we process nsym files.                                    c
c            each file has a t-matrix over all channels(i,l,m)         c
c            for that overall symmetry at a number of                  c
c            different energies.                                       c
c----------------------------------------------------------------------c
      call cardin (card)
      nsym=intkey(card,'no-symmetries',1,' ')
c----------------------------------------------------------------------c
c              write relevant information to rwf file                  c
c----------------------------------------------------------------------c
         call iosys ('write integer "no channels" to rwf',1,nchan,
     1               0,' ')
         call iosys ('write integer "no energies" to rwf',1,nen,0,' ')
         call iosys ('write real "scattering energies" to rwf',nen,
     1               energy,0,' ')
      do 100 isym=1,nsym
         call posinp('$sym-'//itoc(isym),cpass)
         call cardin(card)
         ntchn=intkey(card,'total-no-channels',1,' ')
         if (ntchn.gt.mxchn) then
             call lnkerr('not enough memory')
         endif
         do 200 i=1,nchan
            nlm(i)=intkey(card,'no-lm-ch-'//itoc(i),1,' ')
            call intarr(card,'l-values-ch-'//itoc(i),lch(1,i),
     1                  nlm(i),' ')
            call intarr(card,'m-values-ch-'//itoc(i),mch(1,i),
     1                  nlm(i),' ')
  200    continue
         call iosys ('write integer "total channels" to blkt',1,ntchn,
     1               0,' ')
         call iosys ('write integer "no lm" to blkt',nchan,nlm,0,' ')
         do 300 i=1,nchan
            call iosys ('write integer "l vals ch-'//itoc(i)//
     1                  '" to blkt',nlm(i),lch(1,i),0,' ')
            call iosys ('write integer "m vals ch-'//itoc(i)//
     1                  '" to blkt',nlm(i),mch(1,i),0,' ')
  300    continue
         count=0
         do 400 i=1,nchan
            do 500 j=1,nchan
               count=count+1
               fil(count)='"t  '//itoc(i)//itoc(j)//'"'
               nwords=nen*(2+2*nlm(i)*nlm(j))
               call iosys ('create real '//fil(count)//' on blkt',
     1                     nwords,0,0,' ')
  500       continue
  400    continue
         tmat=ioff
         tij=tmat+2*ntchn*ntchn
c----------------------------------------------------------------------c
c                  begin energy dependent step                         c
c----------------------------------------------------------------------c
         do 20 ene=1,nen
            eau=.5d+00*energy(ene)
            write (iout,1) eau
c----------------------------------------------------------------------c
c              read and print variational t matrices                   c
c----------------------------------------------------------------------c
            filnm='tmat-'//fptoc(energy(ene))
            title='t-matrix'
            call wrmat(z(tmat),ntchn,filnm,title)
c----------------------------------------------------------------------c
c         block according to channels and output                       c
c----------------------------------------------------------------------c
            count=0
            do 30 i=1,nchan
               do 40 j=1,nchan
                  count=count+1
                  ek(1)=sqrt(energy(ene)+2.d+00*(egrnd-echan(i)))
                  ek(2)=sqrt(energy(ene)+2.d+00*(egrnd-echan(j)))
                  call iosys ('write real '//fil(count)//' to blkt '//
     1                        'without rewinding',2,ek,0,' ')
                  n1=nlm(i)
                  n2=nlm(j)
                  call tout (z(tmat),z(tij),ek(1),i,j,nchan,ntchn,
     1                       n1,n2,nlm) 
                  call iosys ('write real '//fil(count)//' to blkt '//
     1                        'without rewinding',2*n1*n2,z(tij),0,' ')
   40          continue
   30       continue
   20    continue
      call iosys ('rewind all on tmat read-and-write',0,0,0,' ')
      call iosys ('close tmat',0,0,0,' ')
      call iosys ('rewind all on blkt read-and-write',0,0,0,' ')
      call iosys ('close blkt',0,0,0,' ')
  100 continue
      call chainx (0)
      stop
    1 format (/,5x,'incident electron energy in hartrees',1x,f15.8)
    2 format (//,10x,'m6016: reformat and print t-matrices')
      end
