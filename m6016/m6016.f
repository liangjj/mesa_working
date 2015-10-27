*deck m6016
c***begin prologue     m6016
c***date written       900703   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6016, link m6016, kohn variational
c***author             schneider, barry (lanl)
c***source             m6016
c***purpose            read matrices from m6009 and reformat
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       m6016
      program tchank (tape55)
      implicit integer(a-z)
      parameter (dime=100, mxchn=50, mxpchn=10)
      real *8 z, zdum, echan, fpkey, egrnd
      real *8 energy, eau, ek
      character *4096 ops
      character *8 cpass, filtmt, blktmt
      character *1600 card
      character *16 chrkey, filnm
      character *3 ans, itoc
      character *10 fptoc
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
      nchan=intkey(card,'no-channels',1,' ')
      call fparr(card,'channel-energies',echan,nchan,' ')
      egrnd=fpkey(card,'reference-energy',echan(1),' ')
      nen=intkey(card,'no-energies',1,' ')
      call fparr(card,'scattering-energies',energy,nen,' ')
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
         filtmt=chrkey(card,'kohn-t-matrix-file-name','tmat',' ')            
         blktmt=chrkey(card,'blocked-t-matrix-file-name','kntmat',' ')
         do 200 i=1,nchan
            nlm(i)=intkey(card,'no-lm-ch-'//itoc(i),1,' ')
            call intarr(card,'l-values-ch-'//itoc(i),lch(1,i),
     1                  nlm(i),' ')
            call intarr(card,'m-values-ch-'//itoc(i),mch(1,i),
     1                  nlm(i),' ')
  200    continue
         call iosys ('open tmat as old',0,0,0,filtmt)
         call iosys ('open blkt as new',262144,0,0,blktmt) 
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
            write (55,1) eau
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
*deck wrmat
c***begin prologue     wrmat
c***date written       900118   (yymmdd)
c***revision date               (yymmdd)
c***keywords           wrmat, link m6015, kohn variational
c***author             schneider, barry (lanl)
c***source             m6015
c***purpose            write kohn matrices
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       wrtmat
      subroutine wrmat(a,n,tit,title)
      implicit integer (a-z)
      real *8  a
      character *16 tit
      character *3 ans
      character *80 title
      dimension a(*)
      common /io/ inp,iout
      call iosys ('does '//tit//' exist on tmat',0,0,0,ans)
      if (ans.eq.'yes') then
          call iosys ('read real '//tit//' from tmat',2*n*n,a,0,' ')
          call prntcm(title,a,n,n,n,n,iout)
      else
          write (iout,10) tit
      endif
      return
   10 format(/,5x,'file:',a16,' not on tmat')
      end
*deck tout
c***begin prologue     tout
c***date written       900703   (yymmdd)
c***revision date               (yymmdd)
c***keywords           tout, link m6016, kohn variational
c***author             schneider, barry (lanl)
c***source             m6016
c***purpose            write channel kohn matrices
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       tout
      subroutine tout(t,tij,ki,i,j,nchan,ntchn,ni,nj,nlm)
      implicit integer (a-z)
      complex*16  t, tij
      character *3 itoc
      character *80 title
      real *8 cross, ki, pi
      dimension t(ntchn,ntchn), tij(ni,nj), nlm(nchan)
      common /io/ inp,iout
      data pi / 3.141592653589793  /
      title='tij-matrix for i= '//itoc(i)//' j= '//itoc(j)
      write (55,1) title
      iloc=0
      do 10 ii=1,i-1
         iloc=iloc+nlm(ii)
   10 continue
      jloc=0
      do 20 jj=1,j-1
         jloc=jloc+nlm(jj)
   20 continue
      icnt=iloc
      do 30 ii=1,ni
         icnt=icnt+1
         jcnt=jloc
         do 40 jj=1,nj
            jcnt=jcnt+1
            tij(ii,jj)=t(icnt,jcnt)
   40    continue
   30 continue
      cross=0.d+00
      do 50 i=1,ni
         do 60 j=1,nj
            cross=cross+abs(tij(i,j))**2
   60    continue
   50 continue
      cross=4.d+00*pi*cross/(ki*ki)  
      call prntcm(title,tij,ni,nj,ni,nj,iout) 
      write (iout,70) cross
      write (55,70) cross
   70 format (//,5x,'total cross section',1x,f20.10)
    1 format (/,5x,a80)
      return
      end
 
 
 
 
 
 
 
