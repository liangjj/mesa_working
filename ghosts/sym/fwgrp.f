*deck fwgrp
      subroutine fwgrp(maxap3,pgrp,natoms,icharg,multip,ian,nop,
     $     maxop,nperm,a,b,d,molfor,lenfor,fwg,length,
     $     iatflg)
      implicit real*8(a-h,o-z)
c
c     given the schonflies symbol for the point group, pgrp,
c     and the re-oriented coordinates, a, determine
c     1--  the molecule's stoichiometry
c     2--  the molecule's framework group
c
      dimension a(maxap3,3), ian(1), b(maxap3,3), d(maxap3,3),
     $     t(3,3), v(3), iatflg(1), nperm(maxap3,1)
      parameter (one=1.0d+00,four=4.0d+00)
      character*2 itoc
      character molfor*(*), fwg*(*)
      character*4 pgrp
      character*1 ichar,jchar
      integer ctoi
c
c     call rtrace(6hfwgrp ,1)
      pi = four * atan(one)
c
c     obtain the stoichiometric formula.
c
c     initialize the buffers which are used to accumulate
c     the character strings and set all atom flags.
c
      length = 0
      lenfor = 0
      molfor=' '
      fwg=' '
      do 21 i=1,natoms
 21      iatflg(i)=2
c
c     get the stoichiometric formula.
c
         call stoich(maxap3,iatflg,ian,a,natoms,lenfor,molfor,0,3)
         call stoich(maxap3,iatflg,ian,a,natoms,lenfor,molfor,0,0)
c
c     remove punctuation.
c
         lenfor = lenfor - 3
         do 23 i=1,lenfor
            molfor(i:i) = molfor(i+1:i+1)
 23      continue
c
c     append the charge, if not zero.
c
         if (icharg .eq. 0) goto 25
         call append('(',molfor,lenfor)
         num = iabs(icharg)
         call append(itoc(num),molfor,lenfor)
         if (icharg .lt. 0) call append('-',molfor,lenfor)
         if (icharg .gt. 0) call append('+',molfor,lenfor)
         if (multip .eq. 1) goto 29
         call append(',',molfor,lenfor)
         goto 27
c
c     append the multiplicity, if not 1.
c
 25      continue
         if (multip .eq. 1) goto 33
         call append('(',molfor,lenfor)
 27      call append(itoc(multip),molfor,lenfor)
 29      call append(')',molfor,lenfor)
c
 33      continue
c
c     length is the current number of characters in the buffer used
c     to accumulate the framework group, fwg.
c
         fwg=pgrp
         length=4
         call append(' ',fwg,length)
         call append('[',fwg,length)
c
c     initialize the atom flags.  three values are possible during the
c     course of things
c     0--  the atom has not been assigned to a symmetric sub-space.
c     1--  the atom has been assigned to a symmetric sub-space.
c     2--  the atom is in the current symmetric sub-space.
c
         do 40 iat=1,natoms
            iatflg(iat) = 0
 40      continue
c
c     branch on the first character of the schonflies symbol.
c
         ichar = pgrp(1:1)
         if(ichar.eq.'c') goto 100
         if(ichar.eq.'s') goto 200
         if(ichar.eq.'d') goto 300
         if(ichar.eq.'t') goto 400
         if(ichar.eq.'o') goto 500
         if(ichar.eq.'i') goto 600
         if(ichar.eq.'k') goto 700
         goto 700
c
c     the first character in the schonflies symbol is c.
c
 100     continue
         jchar = pgrp(2:2)
c
c     test for the special groups c*v, c1, cs, ci, and c2v.
c
c     c*v.
c
         if(jchar .ne. '*') goto 110
         call append('c',fwg,length)
         call append('*',fwg,length)
         do 105 iat=1,natoms
            iatflg(iat) = 2
 105     continue
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,1)
         goto 2000
c
c     c1.
c
 110     continue
         if(pgrp(1:2).ne.'c1') goto 120
         goto 2000
c
c     cs.
c
 120     continue
         if(jchar.ne.'s') goto 130
         call ssssig(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 2000
         call append('sg',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     ci.
c
 130     continue
         if(jchar.ne.'i') goto 140
         call ssso(maxap3,natoms,iatflg,a,itst)
         if(itst.eq.0) goto 2000
         call append('o',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     c2v.
c
 140     continue
         if(pgrp.ne.'c2v') goto 150
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 143
         call append('c',fwg,length)
         call append('2',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,1)
 143     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 147
         call append('sg',fwg,length)
         call append('v',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 147     call ssssig(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 2000
         call append('sg',fwg,length)
         call append('v',fwg,length)
         call append('''',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     get norder and do cn, cnv, and cnh.
c
 150     continue
         norder = ctoi(pgrp(2:3))
         if(norder.lt.10) then
            jchar=pgrp(3:3)
         else
            jchar=pgrp(4:4)
         endif
         if(jchar.ne.' ') goto 160
c
c     cn.
c
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 2000
         call append('c',fwg,length)
         call append(itoc(norder),fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,1)
         goto 2000
c
c     cnh.
c
 160     continue
         if(jchar.ne.'h') goto 170
         call ssso(maxap3,natoms,iatflg,a,itst)
         if(itst.eq.0) goto 163
         call append('o',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 163     call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 167
         call append('c',fwg,length)
         call append(itoc(norder),fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
 167     call ssssig(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 2000
         call append('sg',fwg,length)
         call append('h',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     cnv (n odd).
c
 170     continue
         if(mod(norder,2).eq.0) goto 180
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 175
         call append('c',fwg,length)
         call append(itoc(norder),fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,1)
 175     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 2000
         call append(itoc(norder),fwg,length)
         call append('sg',fwg,length)
         call append('v',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     cnv (n even, > 2).
c
 180     continue
         nhalf = norder / 2
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 183
         call append('c',fwg,length)
         call append(itoc(norder),fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,1)
 183     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 187
         call append(itoc(nhalf),fwg,length)
         call append('sg',fwg,length)
         call append('v',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 187     theta = pi / norder
         call rotate(maxap3,a,b,natoms,t,3,theta)
         call ssssig(maxap3,natoms,iatflg,b,1,itst)
         if(itst.eq.0) goto 2000
         call append(itoc(nhalf),fwg,length)
         call append('sg',fwg,length)
         call append('d',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     sn.
c
 200     continue
         call ssso(maxap3,natoms,iatflg,a,itst)
         if(itst.eq.0) goto 210
         call append('c',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 210     norder = ctoi(pgrp(2:3)) / 2
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 2000
         call append('c',fwg,length)
         call append(itoc(norder),fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
         goto 2000
c
c     the first character in the schonflies symbol is d.
c
 300     continue
         call ssso(maxap3,natoms,iatflg,a,itst)
         if(itst.eq.0) goto 310
         call append('o',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
c
 310     if(pgrp(2:2).ne.'*') goto 320
c
c     d*h.
c
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         call append('c',fwg,length)
         call append('*',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
         goto 2000
c
 320     continue
         norder = ctoi(pgrp(2:3))
         if(norder.lt.10) then
            jchar=pgrp(3:3)
         else
            jchar=pgrp(4:4)
         endif
         nhalf = norder / 2
         if(pgrp.ne.'d2h') goto 338
c
c     d2h.
c
         call sssc(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 323
         call append('c',fwg,length)
         call append('2',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,1,2)
 323     call sssc(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 326
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('''',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
 326     call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 329
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('"',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
 329     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 332
         call append('sg',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 332     call ssssig(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 335
         call append('sg',fwg,length)
         call append('''',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 335     call ssssig(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 2000
         call append('sg',fwg,length)
         call append('"',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
 338     continue
         call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 340
         call append('c',fwg,length)
         call append(itoc(norder),fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
c
 340     continue
         if(jchar.ne.' ') goto 360
         if(mod(norder,2).eq.0) goto 350
c
c     dn (n odd).
c
         call sssc(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 2000
         call append(itoc(norder),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
         goto 2000
 350     continue
c
c     dn (n even).
c
         call sssc(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 355
         if(nhalf.ne.1) call append(itoc(nhalf),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('''',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
 355     theta = pi / norder
         call rotate(maxap3,a,b,natoms,t,3,theta)
         call sssc(maxap3,natoms,iatflg,b,2,itst)
         if(itst.eq.0) goto 2000
         if(nhalf.ne.1) call append(itoc(nhalf),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('"',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,b,natoms,length,fwg,2,2)
         goto 2000
 360     continue
c
c     dnd.
c
         if(jchar.ne.'d') goto 370
         theta = pi / norder
         call rotate(maxap3,a,b,natoms,t,3,theta)
         call sssc(maxap3,natoms,iatflg,b,2,itst)
         if(itst.eq.0) goto 365
         call append(itoc(norder),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('''',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
 365     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 2000
         call append(itoc(norder),fwg,length)
         call append('sg',fwg,length)
         call append('d',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
 370     continue
c
c     dnh (n odd).
c
         if(mod(norder,2).eq.0) goto 380
         call sssc(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 373
         call append(itoc(norder),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
 373     call ssssig(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 377
         call append('sg',fwg,length)
         call append('h',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 377     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 2000
         call append(itoc(norder),fwg,length)
         call append('sg',fwg,length)
         call append('v',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
 380     continue
c
c     dnh (n even, > 2).
c
         call sssc(maxap3,natoms,iatflg,a,2,itst)
         if(itst.eq.0) goto 384
         call append(itoc(nhalf),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('''',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
 384     theta = pi / norder
         call rotate(maxap3,a,b,natoms,t,3,theta)
         call sssc(maxap3,natoms,iatflg,b,2,itst)
         if(itst.eq.0) goto 388
         call append(itoc(nhalf),fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call append('"',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
 388     call ssssig(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 392
         call append('sg',fwg,length)
         call append('h',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 392     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 396
         call append(itoc(nhalf),fwg,length)
         call append('sg',fwg,length)
         call append('v',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 396     call ssssig(maxap3,natoms,iatflg,b,1,itst)
         if(itst.eq.0) goto 2000
         call append(itoc(nhalf),fwg,length)
         call append('sg',fwg,length)
         call append('d',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,b,natoms,length,fwg,0,0)
         goto 2000
c
c     the first character in the schonflies symbol is t.
c
 400     continue
         if(pgrp(2:2).ne.'h') goto 405
         length = length - 2
         do 402 i=3,15
            fwg(i:i) = ' '
 402     continue
         return
 405     call ssso(maxap3,natoms,iatflg,a,itst)
         if(itst.eq.0) goto 410
         call append('o',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
c
c     rotation operations c3 and c2.
c
 410     continue
         call move(maxap3,a,b,natoms)
         v(1) = one
         v(2) = one
         v(3) = one
         call putsym(maxap3,b,d,t,v,natoms,3)
         call sssc(maxap3,natoms,iatflg,b,3,itst)
         if(itst.eq.0) goto 420
         call append('4',fwg,length)
         call append('c',fwg,length)
         call append('3',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
 420     call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 430
         call append('3',fwg,length)
         call append('c',fwg,length)
         call append('2',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
 430     if(pgrp(2:2).eq.' ') goto 2000
c
c     dihedral planes for td.
c
         theta = atan(one)
         call rotate(maxap3,a,b,natoms,t,3,theta)
         call ssssig(maxap3,natoms,iatflg,b,1,itst)
         if(itst.eq.0) goto 2000
         call append('6',fwg,length)
         call append('sg',fwg,length)
         call append('d',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     the first character in the schonflies symbol is o.
c
 500     continue
         call ssso(maxap3,natoms,iatflg,a,itst)
         if(itst.eq.0) goto 510
         call append('o',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
c
c     rotation operations c4, c3, and c2.
c
 510     call sssc(maxap3,natoms,iatflg,a,3,itst)
         if(itst.eq.0) goto 520
         call append('3',fwg,length)
         call append('c',fwg,length)
         call append('4',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
 520     call move(maxap3,a,b,natoms)
         v(1) = one
         v(2) = one
         v(3) = one
         call putsym(maxap3,b,d,t,v,natoms,3)
         call sssc(maxap3,natoms,iatflg,b,3,itst)
         if(itst.eq.0) goto 530
         call append('4',fwg,length)
         call append('c',fwg,length)
         call append('3',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,3,2)
 530     theta = atan(one)
         call rotate(maxap3,a,b,natoms,t,3,theta)
         call sssc(maxap3,natoms,iatflg,b,2,itst)
         if(itst.eq.0) goto 540
         call append('4',fwg,length)
         call append('c',fwg,length)
         call append('3',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,2,2)
c
c     dihedral and horizontal planes of oh.
c
 540     if(pgrp(2:2).ne.'h') goto 2000
         call ssssig(maxap3,natoms,iatflg,b,1,itst)
         if(itst.eq.0) goto 550
         call append('6',fwg,length)
         call append('sg',fwg,length)
         call append('d',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 550     call ssssig(maxap3,natoms,iatflg,a,1,itst)
         if(itst.eq.0) goto 2000
         call append('3',fwg,length)
         call append('sg',fwg,length)
         call append('h',fwg,length)
         call ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
         goto 2000
c
c     the first character in the schonflies symbol is i.
c
 600     continue
         length = length - 2
         i1 = length + 1
         do 620 i=i1,15
            fwg(i:i) = ' '
 620     continue
         return
c
c     the first character in the schonflies symbol is k
c     or is not recognized.
c
 700     continue
         length = length - 2
         i1 = length + 1
         do 720 i=i1,15
            fwg(i:i) = ' '
 720     continue
         return
c
c     fill the asymmetric sub-space.
c
 2000    continue
         num = 0
         do 2020 iat=1,natoms
            if(iatflg(iat).eq.1) goto 2020
            num = num + 1
            iatflg(iat) = 2
 2020    continue
         if(num.eq.0) goto 2040
         call append('x',fwg,length)
         call stoich(maxap3,iatflg,ian,a,natoms,length,fwg,0,0)
 2040    continue
         length=length-1
         call append(']',fwg,length)
         return
         end
