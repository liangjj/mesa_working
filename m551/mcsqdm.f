*deck @(#)mcsqdm.f	1.3  8/3/91
      subroutine mcsqdm (nsym,nobt,nob,nco,ipair,nft,dab,den,buf,
     $     ncor,mblkd,lpqrs,nblock,idcor,linear,ipfdm,iter,ecore,enuc,
     $     ops)
c
c***begin prologue     mcsqdm
c***date written       871022   (yymmdd)
c***revision date      871216   (yymmdd)
c   16 december 1987   bhl      brl
c   code added to compute the true one-electron energy
c   ( no 2j-k terms )
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcsqdm.f	1.3   8/3/91
c
c***purpose
c
c***description
c
c
c-------------------------------------------------
c  a header record for the scaled density matrices
c  is output by this program
c-------------------------------------------------
c
c     mblkd( 1,ipqrs)=idtype  integral block type
c     mblkd( 2,ipqrs)=nu      integral label
c     mblkd( 3,ipqrs)=ml                "    "     "
c     mblkd( 4,ipqrs)=mi                "    "     "
c     mblkd( 5,ipqrs)=mj                "    "     "
c     mblkd( 6,ipqrs)=mk                "    "     "
c     mblkd( 7,ipqrs)=ml                "    "     "
c     mblkd( 8,ipqrs)=njdm     number of coulomb  dm blocks (max of 2)
c     mblkd( 9,ipqrs)=nkdm     number of exchange dm blocks (max of 3)
c     mblkd(10,ipqrs)=mmk1     mmk1 and mml1 pt. to the hessian charge
c     mblkd(11,ipqrs)=mml1      distribution that this dm contributes
c     mblkd(12,ipqrs)=id1      location of this dm
c     mblkd(13,ipqrs)=ntot1    number of elements in this dm
c     mblkd(14,ipqrs)=ncol1    length of a column in the dm
c     mblkd(15,ipqrs)=nrow1    length of a row in the dm
c     mblkd(16,ipqrs)=isqr1    storage mode for dm
c     mblkd(17,ipqrs)=         repeat
c                               "
c                               "
c-----------------------------------
c
c***references
c
c***routines called    (none)
c
c***end prologue       mcsqdm
c
      implicit real*8 (a-h,o-z)
      integer wpadti
      logical debug
      character*(*) ops
      dimension nob(2),ipair(2),dab(2),den(2),buf(2),mblkd(51,2),
     $     lpqrs(2)
      character*8 mcscf
c
      data debug/.false./
c
      common /io/ inp,iout
c
      call iosys('read character mcscf_type from mcscr',-1,0,0,mcscf)
c
c----------------------------
c    build the ipair  pointer
c----------------------------
c
      ix=0
      do 5 i=1,nobt
         ipair(i)=ix
         ix=ix+i
    5 continue
      ipair(nobt+1)=ix
c
c
      norb1=(nobt*(nobt+1))/2
c
c-------------------------------------c
c  square the 1-electron density
c   and double the diagonals
c-------------------------------------c
c
      call iosys('read real "mo 1pdm" from mcscr',norb1,buf,0,' ')
c
      if(debug) then
            write (iout,9101)
 9101       format (5x,'mcscf one-particle density matrix:')
            call print(buf,norb1,nobt,iout)
      end if
c
c    scale the density by .5 to make guga code consistent with ibm
c
      do 9011 i=1,norb1
         buf(i)=buf(i)*.5d+00
 9011 continue
c
c---
c     scale the density temporarily
c---
      ix=0
      do 9001 i=1,nobt
         ix=ix+i
         buf(ix)=buf(ix)*.5d+00
 9001 continue
      do 9002 i=1,norb1
         buf(i)=buf(i)*2.d+00
 9002 continue
cc
      call iosys('write real mcscf_mo_1pdm to rwf',
     $     norb1,buf,0,' ')
c
c     compute the true one-electron energy ( no 2j-k terms )
c
      call iosys('read real mc_h_int from rwf',norb1,den,0,' ')
      e1h=2.d+00*sdot(norb1,den,1,buf,1)
c.h_int
c      e1h=0.0d+00
c      do 6011 i=1,norb1
c         e1h=e1h+den(i)*buf(i)
c 6011 continue
c      e1h=e1h*2.d+00
c.h_int
      call iosys('read real mc_tint1 from mcscr',norb1,den,0,' ')
c
      if(debug) then
          write (iout,9102)
 9102     format (5x,'mcscf one-electron integrals:')
          call print(den,norb1,nobt,iout)
      end if
c
c     trace the density with the transformed integrals
c     as a check
c
      e1=0.0d+00
      do 6001 i=1,norb1
         e1=e1+den(i)*buf(i)
 6001 continue
      e1=e1*2.
c
      ejk=e1-e1h
c
      ix=0
      trace=0.0d+00
      do 7001 i=1,nobt
         ix=ix+i
         trace=trace+buf(ix)
 7001 continue
c
      ia=1
      ib=1
      do 55 i=1,nsym
         if(nob(i).eq.0)go to 55
         if(ipfdm.lt.2) go to 54
         write(iout,53) i
 53      format(/,'  one-particle density for symmetry ',i5)
         call printm(buf(ia),nob(i),1)
 54      continue
         call mcdmsq(buf(ia),dab(ib),nob(i))
         ia=ia+(nob(i)*(nob(i)+1))/2
         ib=ib+nob(i)*nob(i)
 55   continue
c
c-----------------------------------------------
c     read in 2-electron density matrix elements
c-----------------------------------------------
c
      idcor=1
      ipqrs=0
c
c
      m1 = 1
      m2 = 1
      m3 = 1
      m4 = 1
      n1 = nobt
      n2 = nobt
      n3 = nobt
      n4 = nobt
      n12=norb1
      n34=n12
c
ccc   if (n12*n34.eq.0) go to 100
c
      n1234=n12*n34
c
ccc
c
c     ncodi=1 means triangle
c
c     case m1=m2=m3=m4   iiii  idtype=1
c
ccc
 200  continue
c
      ll=0
      ic=1
c
      if(n1234.eq.0) call lnkerr(' ')
c
      ntot= n12*(n12+1)/2
c
      call iosys('read real "mo 2pdm" from mcscr',n12**2,den,0,' ')
c
      if(debug) then
          write (iout,9103)
 9103     format (5x,'mcscf two-particle density matrix:')
          call matout(den,n12,n12,n12,n12,iout)
      end if
c
      call dcan(junk,den,buf,nobt,n12,(n12+1)*n12/2)
c
      call iosys('write real mcscf_mo_2pdm to rwf',ntot,buf,0,' ')
c
c
c    scale the density by .5
c
      do 9022 i=1,ntot
         buf(i)=buf(i)*.5d+00
 9022 continue
c
      ntott=n12*n12
c
      call iosys('read real mc_tint2 from mcscr',ntott,den,0,' ')
c
      if(debug) then
          write (iout,9104)
 9104     format (5x,'mcscf two-electron integrals:')
          call matout(den,n12,n12,n12,n12,iout)
      end if
c
      e2=0.0d+00
      ix=0
      jx=0
      do 6004 i=1,n12
         do 6003 j=1,i
            e2=e2+buf(ix+j)*den(jx+j)
 6003    continue
         ix=ix+i
         jx=jx+n12
 6004 continue
      e2=e2*2.d+00
      etota=e1+e2
      elec=etota+ecore
      etot=elec+enuc
      if(ipfdm.gt.0)then
         call iosys('read real "mcscf core 1e energy" from rwf',
     $   1,e1c,0,' ')
         call iosys('read real "mcscf core 2e energy" from rwf',
     $   1,e2c,0,' ')
c
       e1t=e1c+e1h
       e2t=e2c+e2+ejk
c
         write(iout,6006) iter,e1,e1h,ejk,e2,etota,ecore,e1c,e2c,elec,
     $         e1t,e2t,enuc,etot
c
 6006    format(//,' **** energy computed from the density ',
     $        'in iteration ',i4,//,
     $        ' active  one-electron energy ',2x,f20.12,/,
     $        '         h(1)         energy ',2x,f20.12,/,
     $        '         2j-k         energy ',2x,f20.12,/,
     $        ' active  two-electron energy ',2x,f20.12,/,
     $        ' active  electronic   energy ',2x,f20.12,/,
     $        ' core    electronic   energy ',2x,f20.12,/,
     $        '         h(1)         energy ',2x,f20.12,/,
     $        '         2j-k         energy ',2x,f20.12,/,
     $        ' total   electronic   energy ',2x,f20.12,/,
     $        '         h(1)         energy ',2x,f20.12,/,
     $        '         g(2)         energy ',2x,f20.12,/,
     $        ' nuclear repulsion    energy ',2x,f20.12,/,
     $        ' *total* mcscf        energy ',2x,f20.12)
      endif
c
cps      write (iout,16006) etot,iter
cps16006 format (/,' *total* mcscf energy:',f16.8,' on iteration ',
cps     1     i3,/)
cps      endfile iout
cps      backspace iout
c
c
      call iosys('write real energy to rwf',1,etot,0,' ')
      call iosys('write real "mcscf energy" to rwf',1,etot,0,' ')
      call iosys('read real "mcscf core 1e energy" from rwf',
     $     1,ecore,0,' ')
      call iosys('write real "mcscf one-electron energy" to rwf',
     $     1,e1+ecore,0,' ')
      call iosys('read real "mcscf core 2e energy" from rwf',
     $     1,ecore,0,' ')
      call iosys('write real "mcscf two-electron energy" to rwf',
     $     1,e2+ecore,0,' ')
c
      ndex=(m1*(m1+1))/2
      ndex=(ndex*(ndex+1))/2
c
      ipqrs=ipqrs+1
      lpqrs(ndex)=ipqrs
c
      n1=nobt
      nii=(n1*(n1+1))/2
      nij=n1*n1
c
      ntot1=nii*nii
      ntot3=nii*nij
c
      id1=idcor
      idk=id1+ntot1
      id3=idk+ntot3
      idx=id3+ntot1
      ineed=wpadti(idx+ntot1)
c
      idcor=id3
c
c     ----- dynamic core allocation by pws -----
c
      call getscm(ineed,den,ncor,'mc square up density',0)
c
cpsdyn      if(iadtwp(ineed).gt.ncor) go to 3000
cc
      if(n1234.eq.0)go to 293
c
c      write(iout,66006)id1,idk,id3,idx,n1,nii
66006 format(' id1 idk id3 idx n1 nii ',/,6(1x,i5))
      call mciiii(den(id1),den(idk),den(id3),den(idx),n1,nii,buf,ipair)
c
 293  continue
cc
      if(ipfdm.gt.1)then
         write(iout,7002) trace
 7002    format(/,'  trace of one-particle density matrix ',f8.3)
      endif
c
      if(ipfdm.lt.3)go to 9801
      write(iout,9901) ipqrs
 9901 format(/,' second order dm:  mciiii ',i4)
      call printm(buf,nii,1)
c
      write(iout,*)' scaled 2pdm..d1 '
      call printr(den(id1),nii,nii)
      write(iout,*)' scaled 2pdm..d3 '
      call printr(den(id3),nii,nii)
      write(iout,*)' scaled 2pdm..dx '
      call printr(den(idx),nii,nii)
      write(iout,*)' scaled 2pdm..dk '
      call printr(den(idk),nii,nij)
 9801 continue
c label
c..rlm 
c     nu and ml are undefined at this point. this array may be used and
c     so i am setting them explicitly to zero.
      nu=0
      ml=0
c
      mblkd( 1,ipqrs)=1
      mblkd( 2,ipqrs)=nu
      mblkd( 3,ipqrs)=ml
      mblkd( 4,ipqrs)=m1
      mblkd( 5,ipqrs)=m2
      mblkd( 6,ipqrs)=m3
      mblkd( 7,ipqrs)=m4
      mblkd( 8,ipqrs)=1
      mblkd( 9,ipqrs)=1
c dj1
      mblkd(10,ipqrs)=m1
      mblkd(11,ipqrs)=m1
      mblkd(12,ipqrs)=id1
      mblkd(13,ipqrs)=ntot1
      mblkd(14,ipqrs)=nii
      mblkd(15,ipqrs)=nii
      mblkd(16,ipqrs)=00
c dj2
      mblkd(17,ipqrs)=0
      mblkd(18,ipqrs)=0
      mblkd(19,ipqrs)=0
      mblkd(20,ipqrs)=0
      mblkd(21,ipqrs)=0
      mblkd(22,ipqrs)=0
      mblkd(23,ipqrs)=0
c dk1
      mblkd(24,ipqrs)=m1
      mblkd(25,ipqrs)=m1
      mblkd(26,ipqrs)=idk
      mblkd(27,ipqrs)=ntot3
      mblkd(28,ipqrs)=nii
      mblkd(29,ipqrs)=nij
      mblkd(30,ipqrs)=0
c dk2
      mblkd(31,ipqrs)=0
      mblkd(32,ipqrs)=0
      mblkd(33,ipqrs)=0
      mblkd(34,ipqrs)=0
      mblkd(35,ipqrs)=0
      mblkd(36,ipqrs)=0
      mblkd(37,ipqrs)=0
c dk3
      mblkd(38,ipqrs)=0
      mblkd(39,ipqrs)=0
      mblkd(40,ipqrs)=0
      mblkd(41,ipqrs)=0
      mblkd(42,ipqrs)=0
      mblkd(43,ipqrs)=0
      mblkd(44,ipqrs)=0
c dk4
      mblkd(45,ipqrs)=0
      mblkd(46,ipqrs)=0
      mblkd(47,ipqrs)=0
      mblkd(48,ipqrs)=0
      mblkd(49,ipqrs)=0
      mblkd(50,ipqrs)=0
      mblkd(51,ipqrs)=0
c
 2000 continue
c
c      call jtime(it2)
c      time = (it2 - it1) / 100.d0
c      write (iout,9099) time
 9099 format('0   ****** time in this section ',10x,f10.2,' secs')
c
      return
 3000 continue
      write(iout,3010) idcor,ncor
 3010 format(//,'  insufficient core for mcscdm  ineed ihave',2(2x,i8))
      call lnkerr(' ')
      end
