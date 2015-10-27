*deck @(#)detdrv.f	5.1  11/6/94
      subroutine detdrv(joutfg,iopn,jdbl,ms,ifiga,ifigb,ifig,isuma,
     1  isumb,nphase,idstor,msymtp)
c
c  loop over the spatial configurations and check their spatial
c  symmetry against the double-group symmetry to determine which
c  ms determinant routine to call.
c
c  variables to be calculated:
c
c     iopn       list of singly occupied orbitals
c     jdbl       list of doubly occupied orbitals
c     ndeti      number of determinants for configuration i
c     ndet       number of determinants generated
c     nsefi      number of double-group functions for configuration i
c     nsef       number of double-group functions generated
c     ncpi       start of the transformation of determinants to
c                double-geoup-adapted functions
c     ms         spin projection for a given determinant, in the
c                form (mxopn+2*ms)/2
c     ifiga      list of alpha-spin orbitals in a given determinant,
c                doubly occupied given first in list
c     ifigb      list of beta-spin orbitals in a given determinant,
c                doubly occupied given first in list
c     isuma      sum of the orbital numbers for the alpha list for a
c                given determinant
c     isumb      sum of the orbital numbers for the beta list for a
c                given determinant
c     nphase     phase of a given determinant
c
      character*80 rtitle, ctitle, blabel
c
      common /headng/ rtitle, ctitle, blabel
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c2/ msbas,jsym,nopen,kdbl,ndeti,nsefi,ncpi
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
c
      dimension joutfg(nbf),iopn(*),jdbl(*),ms(*),ifiga(mxocc,*),
     1  ifigb(mxocc,*),ifig(*),isuma(*),isumb(*),nphase(*),idstor(*)
c
      rewind iunts1
      write (iw,800) ctitle
c
      do 100 ii = 1,maxdet
      do 100 i = 1,mxocc
  100 ifiga(i,ii) = 0
      do 110 ii = 1,maxdet
      do 110 i = 1,mxocc
  110 ifigb(i,ii) = 0
c
      msbas = mxopn/2
      idspc = 1
      idspcu = 1
      ndet = 0
      nsef = 0
      do 260 ispc = 1,ntotfg
c
c  read in spatial configuration
c
      read (iunt1a,'(80i1)') joutfg, jsym
c
c     compile the list of singly and doubly occupied orbitals
c
      nelu = 0
      nopen = 0
      kdbl = 0
      do 140 i=1,nbf
      if( joutfg(i) - 1 )  140,120,130
c
c  the orbital i is singly occupied
c
  120 nopen = nopen + 1
      iopn(nopen) = i
      nelu = nelu + 1
      go to 140
c
c  the orbital i is doubly occupied
c
  130 kdbl = kdbl + 1
      jdbl(kdbl) = i
      nelu = nelu + 2
  140 continue
c
c  check to see if the configuration has the correct number of
c  electrons and that the number of open shells is within bounds
c
      if( nelu.ne.nel ) then
        write (iw,910) joutfg
        write (iw,920) nelu, nel
        stop
      endif
      if( nopen.gt.mxopn ) then
        write (iw,910) joutfg
        write (iw,930) mxopn
        stop
      endif
      if( mod((nel - nopen),2) .ne. 0 ) then
        write (iw,940) ispc, joutfg
        write (iw,950) nopen
        stop
      endif
c
      if(iciwrt.ge.3) then
        write (iw,820) ispc, jsym, joutfg
        write (iw,830)
      endif
c
      if( msymtp.eq.1 ) then
        ndeti = 2**nopen
      else
        ndeti = 2**max(nopen-1,0)
      endif
      if( msymtp.le.2 ) then
        nsefi = ndeti
      else
        nsefi = 4**((nopen-1)/2)
      endif
      nsef = nsef + nsefi
c
      if(ihamrd.ge.1) go to 220
c
c  select the spin functions for this spatial symmetry
c
      if(mel.eq.1) then
        ncpi=1
        if(mod((jsym+1)/2,2).eq.0) then
          call msfdet(iopn,jdbl,ms,ifigb,ifiga,-1)
        else
          call msfdet(iopn,jdbl,ms,ifiga,ifigb,1)
        endif
      else
        if(msymtp.ge.2.and.mod(jsym,2).ne.(mod(ksym,2))) then
          ncpi=2
        else
          ncpi=1
        endif
        if(msymtp.eq.1) then
          call msidet(iopn,jdbl,ms,ifiga,ifigb)
        elseif(msymtp.eq.2) then
          if(jsym.ne.ksym) then
            call msodet(iopn,jdbl,ms,ifiga,ifigb)
          else
            call msedet(iopn,jdbl,ms,ifiga,ifigb)
          endif
        else
          if((jsym+1)/2.ne.((ksym+1)/2)) then
            call msodet(iopn,jdbl,ms,ifiga,ifigb)
          else
            call msedet(iopn,jdbl,ms,ifiga,ifigb)
          endif
        endif
      endif
c
      call sumphz(ms,ifiga,ifigb,ifig,isuma,isumb,nphase)
c
  220 call detout(joutfg,iopn,ms,ifiga,ifigb,isuma,isumb,nphase,
     2  idstor,ispc)
c
      idspcu = max(idspcu,idspc)
      if(ispc.eq.1) idfrst = idspc
c
      if( (iciwrt.lt.3).or.(ihamrd.ge.1) ) go to 260
c
      if(nopen.eq.0) then
        write (iw,840) nopen
      else
        write (iw,850) nopen, (iopn(ii), ii=1,nopen)
      endif
      write (iw,860) ncpi, nsefi, ndeti
      do 250 i=1,ndeti
      nbeta = mxocc - ms(i)
      nalpha = nel - nbeta
      iphase = 1 - (nphase(i) + nphase(i))
      write (iw,870) (ifiga(jj,i), jj=1,nalpha)
      write (iw,880) (ifigb(jj,i), jj=1,nbeta)
  250 write (iw,890) isuma(i), isumb(i), iphase
c
      write (iw,*) ' '
c
  260 ndet=ndet+ndeti
c
      idstor(idspc) = 0
      write (iunts1) (idstor(i), i=1,idspc)
      rewind iunts1
c
      return
c
  800 format(///5x,27hconfiguration tape label - ,a80)
  820 format('1'/5x,'configurations'/
     1           5x,'number  spatial sym.  occupation numbers'//
     2           i11,i8,8x,80i1/ (27x,80i1) )
  830 format(23x,'list of determinants '/)
  840 format(/23x,'number of open shells =',i5)
  850 format(/23x,23hnumber of open shells =,i5,5x,
     2             21hopen-shell orbitals =,10i5)
  860 format(/23x,42hstart of double-group function expansion =, i5
     2       /23x,34hnumber of double-group functions =, i5
     3       /23x,34hnumber of determinants           =, i5)
  870 format(/23x,12halpha list =,20i3/ (35x,20i3) )
  880 format(23x,12hbeta list  =,20i3/ (35x,20i3) )
  890 format(/23x,7hisuma =,i5,5x,7hisumb =,i5,5x,8hnphase =,i5)
  910 format(///5x,'*****  deleted basic configuration  *****'
     2        // (5x,80i1) )
  920 format(//8x,'number of electrons found    =',i5
     2        /8x,'number of electrons expected =',i5)
  930 format(//6x,'number of open shells too large',i3)
  940 format(///5x,20h*****  configuration,i5,7h  *****// (80i1) )
  950 format(//5x,37hhas the wrong number of open shells (,i2,1h) )
c
      end
