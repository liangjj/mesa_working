*deck @(#)pm904.f	5.1  11/6/94
c
c*******************************************************************
c
c   this computer program contains work performed partially at
c   lawrence livermore national laboratory under
c   the auspices of the office of basic energy sciences,
c   division of chemical sciences, u.s. department of energy,
c   under contract w-7405-eng-48.
c
c   these programs may not be (re)distributed without the
c   written consent of the columbus program group.
c
c   since these programs are under development, correct results
c   are not guaranteed.
c
c*******************************************************************
c
      subroutine pm904
c
c pitzer's spin-orbit configuration interaction modified to 
c fit into mesa, m. braunstein nov. 91
c
c
c  spin-orbit configuration-interaction program
c
c  version log:
c  03-feb-90 btime improvements, xuflow (rmp,rbr,jwh)
c  18-aug-89 mdc version (rmp)
c  1984-1986 written, starting with cit/lasl ci program (rmp,nww)
c
c  cmdc info:
c  keyword   description
c  -------   -----------
c  vax       vax code.
c  cray      cray code. includes "d" exponents in constants.
c  ibm       ibm code.
c  sun       sun code.
c  stellar   stellar code.
c  ibmmvs    mvs specific code.
c  crayctss  ctss specific code.
c
c  dimension requirements
c
c      nbf       number of orbitals
c     nroots     number of roots of ci matrix
c      nsef      number of double-group functions
c     mxopn      maximum number of open shells
c     mxocc      maximum number of occupied orbitals
c     maxdet     maximum number of determinants/configuration
c     ndbgu      maximum number of double-group functions/configuration
c     maxne0     maximum number of nonzero hamiltonian matrix
c                elements per row
c
c     --------------------------------------------------------------
c
c     valone     nbf*(nbf + 1)/2
c     valtwo     maximum number of two-electron integrals (canonical)
c
c     vector     nroots*nsef
c       e        nroots
c
c     joutfg     2*nbf
c      iopn      2*mxopn
c      jdbl      mxocc
c       ms       2*maxdet
c     ifiga      2*mxocc*maxdet
c     ifigb       (same as ifiga)
c      ifig      2*nbf*maxdet
c     isuma       (same as ms)
c     isumb       (same as ms)
c     nphase      (same as ms)
c
c     idstor     (nbf + n2 - 1)/n2 + (mxopn + n8 - 1)/n8 +
c                (maxdet + n4 - 1)/n4 + 2*((maxdet + n16 - 1)/n16) +
c                (maxdet + n1 - 1)/n1 +
c                2*((maxdet*mxocc + n8 - 1)/n8) + 6
c
c     number     ndbgu
c      coef      maxdet*ndbgu or 2*maxdet*ndbgu depending on msymtp
c     detdet     maxdet**2
c     detsef     maxdet*ndbgu
c     sefsef     ndbgu**2
c      vnt       max(3,((mxopn - 1)*mxopn)/2 + 1)
c
c      row       maxne0 + 1
c     cimat      ndbgu*maxne0
c
c     thresh      (same as e)
c     coeff      maximum number of dominant configurations/root (10)
c     idomnt      (same as coeff)
c       p        nroots*(nroots + 1)/2
c       q         (same as p)
c       u         (same as vector)
c       y        nroots**2
c     almax       (same as e)
c
c     istsef     nsef + 1
c      indx      nsef
c     jfigs      (nbf + n2 - 1)/n2*maximum number of energy-sorted
c                configurations
c     jsyms      maximum number of energy-sorted configurations
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      character*24 dattim
      character*80 rtitle, ctitle, blabel, blank
c
      common /headng/ rtitle, ctitle, blabel
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
c  mesa
c
      common/io/inp,iout
c
c  iaup      maximum total array space, in working-precision words
      parameter (iaup = 690884)
      dimension a(iaup), ia(8)
      equivalence (a,ia), (a(2),ib)
c
c  integer word-length divided by 2
c
*mdc*if cray
*      parameter (n2=32)
*mdc*else
      parameter (n2=16)
*mdc*endif
c
      data blank/'        '/, dattim/'        '/
c
c  unit and file assignments
c
*mdc*if crayctss
*      call link('unit99=tty//')
*mdc*endif
c     5         unit for input
      open (5,file='ciin',status='old')
c     iw        unit for output
c     iw     =  6
c  mesa
      iw=iout
c
*mdc*if ibmmvs
*c  use preconnected unit 6
*mdc*else
c     open (iw,file='ciout',status='unknown')
*mdc*endif
c     iunt1a    unit on which the configurations reside
c     iunt1a =  (defaults to 1)
c     iunt1b    unit on which the molecular integrals reside
      iunt1b =  3
c     iunt2a    unit onto which the hamiltonian matrix is written
      iunt2a =  23 
c     icipun    unit onto which the ci vectors are written
c     icipun =  (value read in, if any)
c     iunts1    scratch unit for storing initial configurational data
      iunts1 = 24
      open (iunts1,status='scratch',form='unformatted')
c     iunts2    scratch unit for storing selected configurational data
      iunts2 = 25
c
*mdc*if ibm
*      call xuflow(0)
*mdc*endif
c
c  ratio of working-precision word-length to integer word-length
c
      do 10 i=1,8
   10 ia(i) = i
      irt = ib - 1
c
c  conversion factor from hartrees to electron volts
c
      ev = 27.2113961d0
c
c  program limitations
c
c     maxopn       maximum number of open-shell electrons
      maxopn =      8
c     maxit        maximum nunber of diagonalization iterations
      maxit  =    100
c     maxesc       maximum number of energy-sorted configurations
      maxesc =    500
c     iau          maximum total array space, in working-precision words
      iau    =   iaup
c
c  read in a title for the job
c
   15 read (5,'(a80)') rtitle
c
c  check for a blank indicating the end of the data
c
c mesa
c      if(rtitle.eq.blank) then
c        write (iw,1050)
c        stop 'end of cidbg'
c      endif
c
      write (iw,20)
   20 format('1',28x,'program "cidbg" 4.0b2'/
     1  29x,'columbus program system'/
     2  29x,'spin-orbit configuration-interaction program'//
     3  29x,'version date: 03-feb-90'//
     &  ' maintained by:'/
     &  t5,'russell m. pitzer'/
     &  t5,'department of chemistry'/
     &  t5,'the ohio state university'/
     &  t5,'columbus, oh  43210'/
     &  t5,'bitnet: ts0775 at ohstmvsa')
      call btime(dattim,bsec)
      write (iw,'(/5x,a80,5x,a24)') rtitle, dattim
c
c  input data (excluding configuration, diagonalization, and
c  analysis threshold data.  see detdrv, diagnl, and analyz.)
c
c     nbf        number of orbitals
c     nel        number of electrons
c     ksym       double-group symmetry of state

c                     ksym = 1 to  2 for c1
c                     ksym = 1 to  3 for cs, c2
c                     ksym = 1 to  5 for c2v, d2
c                     ksym = 1 to 10 for d2h
c     ntotfg     number of spatial configurations
c     nroots     number of roots of the ci matrix
c     iciwrt     ci print option
c                     iciwrt = 0     production run, minimum output
c                            = 1     prints all of ci vectors
c                            = 2     prints the ci matrix
c                            = 3     prints the complete list of ci data
c     icipun     ci vector write option
c                     icipun = 0     do not write the ci vectors
c                            > 0     write the ci vectors on unit icipun
c     intwrt     integral write option
c                     intwrt = 0     production run, minimum output
c                            = 1     writes out the one-electron
c                                    integrals
c                            = 2     writes out the one- and two-
c                                    electron integrals
c     ihamrd     hamiltonian read option
c                     ihamrd = 0     the hamiltonian matrix is computed
c                            = 1     the hamiltonian matrix is read
c                                    in from iunt2a
c     ianalz     vector and energy statistics option
c                     ianalz = 0     no vector or energy breakdown
c                            = 1     provide vector and energy breakdown
c     maxne0     maximum number of nonzero hamiltonian matrix elements
c                per row; ignored unless ihamrd = 1, in which case the
c                value is obtained from the previous run's output
c
      read (5,'(16i5)') nroots, iciwrt, icipun, intwrt, ihamrd,
     1  ianalz, iunt1a, maxne0
      if( (ihamrd.ge.1).and.(maxne0.le.0) ) then
        write (iw,30)
   30     format(/'provide value of maxne0, from previous output, ',
     1            'when ihamrd > 0')
        stop
      endif
      if( iunt1a .eq. 0 ) iunt1a = 1
      open (iunt1a,file='configs')
      read (iunt1a,'(a80)') ctitle
      read (iunt1a,'(16i5)') nbf, nel, ksym, ntotfg, mxopn, msymtp
      write (iw,1300) nbf,iciwrt,nel,icipun,ksym,intwrt,ntotfg,
     2  ihamrd,nroots
c
      if( (nbf.le.0).or.(ntotfg.le.0).or.(nroots.le.0) ) then
        write (iw,35)
   35     format(/'bad input data.  see values printed above.')
        stop
      endif
      if( mxopn.gt.maxopn ) then
        write (iw,40) mxopn
   40     format(/'number of open shells too large',i3)
        stop
      endif
      nw = (nbf-1)/n2 + 1
      mel = mod(nel,2)
c
c  calculate maximum number of occupied orbitals, maximum number
c  of determinants/configuration
c
      mxocc = (nel + mxopn)/2
      if( msymtp .eq. 1 ) then
        maxdet = 2**mxopn
      else
        maxdet = 2**max(mxopn-1,0)
      endif
c
      call btime(dattim,time1)
c
      koutfg = 1
      jopn   = koutfg + nbf
      idbl   = jopn   + mxopn
      ns     = idbl   + mxocc
      jfiga  = ns     + maxdet
      jfigb  = jfiga  + (mxocc*maxdet)
      jfig   = jfigb  + (mxocc*maxdet)
      jsuma  = jfig   + nbf + nbf
      jsumb  = jsuma  + maxdet
      jphase = jsumb  + maxdet
      jdstor = jphase + maxdet
      iauu   = (jdstor+irt-1)/irt
      if(iauu.gt.iau) then
        write (iw,45) iauu
   45     format(/'more space needed for detdrv.  iaup at least',i8)
        stop
      endif
      inextu = irt*iau-jdstor+1
      call detdrv(ia(koutfg),ia(jopn),ia(idbl),ia(ns),ia(jfiga),
     1  ia(jfigb),ia(jfig),ia(jsuma),ia(jsumb),ia(jphase),ia(jdstor),
     2  msymtp)
      close (iunt1a)
      ijwrd = (jdstor + idspcu + irt - 2)/irt
c
      call btime(dattim,time2)
c
c  prepare masks for column-index packing and unpacking
c
      fn = (nsef + nsef + 1)
      nbcol = (log(fn)/log(2.0d0))
      mask = 2**nbcol - 1
*mdc*if cray
*      nmask = compl(mask)
*mdc*else
      nmask = not(mask)
*mdc*endif
c
c  calculate maximum number of double-group functions/configuration
c
      if( msymtp .le. 2 ) then
        ndbgu = maxdet
      else
        ndbgu = 4**((mxopn-1)/2)
      endif
c
      write (iw,1400) ntotfg, nsef, ndet, nbcol
      if(iciwrt.ge.3) write (iw,1350) mxopn, ndbgu
c
      if(ihamrd.ge.1) go to 60
c
      call btime(dattim,time3)
c
c      open (iunt1b,file='moints',status='old',form='unformatted')
c
c  read the tape title and core energy from tape iunt1b
c
c      read (iunt1b) blabel, repnuc, efrzc, intmxo, ihso
c mesa read
c
      call iosys('read real "frozen core energy" from rwf',
     $            1,efrzc,0,' ')
      call iosys('read real "nuclear repulsion energy" from rwf',
     $            1,repnuc,0,' ')
      
      ecore=repnuc+efrzc
      write (iw,*) 'core energy, repnuc + efrzc = ', ecore
c
c      ival1  = 1
c      limit1 = (nbf*(nbf+1))/2
c      ival2  = ival1 + limit1
c      limit2 = (limit1*(limit1+1))/2
c      iv     = ival2  + limit2
c      jbl    = irt*(iv + intmxo -1) + 1
c      iauu   = (jbl + intmxo - 2)/irt  + 1
c
c new memory allocation for mesa
c
      ival1  = 1
      limit1 = (nbf*(nbf+1))/2
      ival2  = ival1 + limit1
      limit2 = (limit1*(limit1+1))/2
      ival3  = ival2 + limit2
      ival4  = ival3 + limit1
      ival5  = ival4 + limit1
      iv     = ival5 + limit1
      jbl    = irt*(iv + intmxo -1) + 1
      iauu   = (jbl + intmxo - 2)/irt  + 1
      if(iauu.gt.iau) then
        write (iw,50) iauu
   50     format(/'more space needed for rdints.  iaup at least',i8)
        stop
      endif
c      call rdints(a(ival1),a(ival2),a(iv),ia(jbl))
c
c mesa
c
      call rdints(a(ival1),a(ival2),a(ival3),a(ival4),a(ival5))
c
c      close (iunt1b)
      jkwrd = iauu
c
      call btime(dattim,time4)
c
      koutfg = irt*(iv - 1) + 1
      jopn   = koutfg + nbf + nbf
      ns     = jopn   + mxopn + mxopn
      jfiga  = ns     + (maxdet + maxdet)
      jfigb  = jfiga  + (mxocc*(maxdet + maxdet))
      jfig   = jfigb  + (mxocc*(maxdet + maxdet))
      jsuma  = jfig   + nbf*(maxdet + maxdet)
      jsumb  = jsuma  + (maxdet + maxdet)
      jphase = jsumb  + (maxdet + maxdet)
      jdstor = jphase + (maxdet + maxdet)
      mumber = jdstor + idspcu + idspcu
      icoef  = (mumber + ndbgu + irt - 2)/irt + 1
      idtdt  = icoef  + ndbgu*maxdet
      if(msymtp.ge.2) idtdt = idtdt + ndbgu*maxdet
      ivnt   = idtdt  + maxdet*maxdet
      if(mel.eq.0) ivnt = ivnt + ndbgu*(maxdet + ndbgu)
      irow   = ivnt   + max(2,((mxopn-1)*mxopn)/2) + 2
crlm
c     irow=irow+10
crlm
      iauu   = irow
      if(iauu.gt.iau) then
        write (iw,55) iauu
   55     format(/'more space needed for spinci.  iaup at least',i8)
        stop
      endif
      intmax = iau - irow + 1
      open (iunt2a,file='hmatrix',form='unformatted')
      call spinci(a(ival1),a(ival2),a(ival3),a(ival4),a(ival5),
     1  ia(koutfg),ia(jopn),ia(ns),ia(jfiga),ia(jfigb),ia(jfig),
     2  ia(jsuma),ia(jsumb),ia(jphase),ia(jdstor),ia(mumber),
     3  a(icoef),a(idtdt),a(ivnt),a(irow),msymtp)
      klwrd = irow + intmax - 1
c
   60 irow   = 1
      call hrdwrt(a(irow))
c
      call btime(dattim,time5)
c
      ivect  = 1
      ie     = ivect  + (nroots*nsef)
      irow   = ie     + nroots
      ithrsh = irow   + maxne0
      icoeff = ithrsh + nroots
      jdomnt = irt*(icoeff + 9) + 1
      ip     = icoeff
      iq     = ip     + ((nroots*(nroots + 1))/2)
      iu     = iq     + ((nroots*(nroots + 1))/2)
      iy     = iu     + (nroots*nsef)
      ialmax = iy     + nroots*nroots
      iauu   = max((jdomnt + 9)/irt,ialmax + nroots - 1)
      if(iauu.gt.iau) then
        write (iw,65) iauu
   65     format(/'more space needed for diagnl.  iaup at least',i8)
        stop
      endif
      call diagnl(a(ivect),a(ie),a(irow),a(ithrsh),a(icoeff),
     1  ia(jdomnt),a(ip),a(iq),a(iu),a(iy),a(ialmax))
      write(iw,*)'vector', (a(ivect+i-1),i=1,2)
      write(iw,*)'eigenv',(a(ie+i-1),i=1,2)
      lmwrd = iauu
c
      call btime(dattim,time6)
c
      ijsec = (time2 - time1)
      ijmin = ijsec/60
      ijsec = mod(ijsec,60)
      write (iw,1500) ijmin, ijsec, ijwrd
c
      if(ihamrd.ge.1) go to 70
c
      jksec = (time4 - time3)
      jkmin = jksec/60
      jksec = mod(jksec,60)
c
      klsec = (time5 - time4)
      klmin = klsec/60
      klsec = mod(klsec,60)
c
      write (iw,1525) jkmin, jksec, jkwrd, klmin, klsec, klwrd
c
   70 lmsec = (time6 - time5)
      lmmin = lmsec/60
      lmsec = mod(lmsec,60)
      write (iw,1550) lmmin, lmsec, lmwrd
c
      if(iciwrt.eq.0) go to 80
      koutfg = irt*(irow - 1) + 1
      jdstor = koutfg + nbf
      iauu   = (jdstor + nw + irt)/irt
      if(iauu.gt.iau) then
        write (iw,75) iauu
   75     format(/'more space needed for output.  iaup at least',i8)
        stop
      endif
      call output(a(ivect),a(ie),ia(koutfg),ia(jdstor))
c
   80 ih     = ithrsh
      idele  = ih     + nsef
      koutfg = irt*(idele + ntotfg - 1) + 1
      jdstor = koutfg + nbf
      jstsef = jdstor + nw
      jndx   = jstsef + ntotfg + 1
      kfigs  = jndx   + ntotfg
      ksyms  = kfigs  + nw*maxesc
      iauu   = (ksyms + maxesc + irt - 2)/irt
      if(iauu.gt.iau) then
        write (iw,85) iauu
   85     format(/'more space needed for analyz.  iaup at least',i8)
        stop
      endif
      call analyz(a(ivect),a(ie),a(irow),a(ih),a(idele),ia(koutfg),
     2  ia(jdstor),ia(jstsef),ia(jndx),ia(kfigs),ia(ksyms))
c
c extra code to analyze spin-orbit matrix elements
c
c check for memory
c
      read(5,*)ido
      if(ido.gt.-1)then
        ic=1
        id=ic+(nsef*nroots)
        ihb=id+(nroots*nroots)
        ihr=ihb+(nsef*nsef)
        icthc=ihr+(nroots*nroots)
        iscr=icthc+(nroots*nroots)
        iroot=iscr+5*nroots
        itemp=iroot+nroots
        itemp2=itemp+(nsef*nroots)
        irow=itemp2+(nroots*(nroots+1))/2
        ibfl=irow+nsef
        iauu=ibfl+nsef 
        if(iauu.gt.iau) then
          write(iw,*)'more space needed soanal'
          stop
        end if

        write(iw,*)'doing spin-orbit analysis, step ',ido

        call soanal(a(irow),ido,a(ic),a(id),a(ihb),a(ihr),
     >          a(icthc),a(iscr),a(iroot),a(itemp),a(itemp2),a(ibfl))
      end if
      close (iunt2a)
c
c exit gracefully
c
      call chainx(0)
c
 1000 format(///5x,'integral tape label - ',a80
     1       ///5x,'core energy =',f16.8)
 1050 format(///5x,34h*****  execution terminated  *****,
     2       //16x,11hend of data)
 1300 format( //5x,37hnumber of orbitals...................,i5,
     1         10x,37hci print option......................,i5,
     2        //5x,37hnumber of electrons..................,i5,
     3         10x,37hci vector write option...............,i5,
     4        //5x,37hstate symmetry.......................,i5,
     5         10x,37hintegral write option................,i5,
     6        //5x,37hnumber of spatial configurations.....,i5,
     7         10x,37hhamiltonian read option..............,i5,
     8        //5x,37hnumber of roots of ci matrix.........,i5)
 1350 format(///5x,'maximum number of open shells found =',i5/
     1  5x,'number of double-group functions/configuration =',i5)
 1400 format(//10x,48hnumber of configurations used in the calculation,
     2       //15x,37hnumber of spatial configurations    =,i6,
     3       //15x,37hnumber of double-group functions    =,i6,
     4       //15x,37hnumber of determinants              =,i6,
     5       //15x,37hnumber of bits for column indices   =,i6)
 1500 format('1',4x,'time, space information (minutes, seconds, ',
     1  'working-precision words)'///
     2  10x,'set up determinants  ',i6,' min',i5,' sec',i8,' words')
 1525 format(/
     1  10x,'read in integrals    ',i6,' min',i5,' sec',i8,' words'/
     2 /10x,'set up ci matrix     ',i6,' min',i5,' sec',i8,' words')
 1550 format(/
     1  10x,'diagonalize ci matrix',i6,' min',i5,' sec',i8,' words')
c
      return 
      end
