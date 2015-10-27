*deck @(#)group.f	5.1  11/6/94
      subroutine group(grp,ngen,nirrep,lirrep,gengrp,naxis,axis,nop,
     #     lambda,labop,prtsym,redsym,subgrp)
c
c***begin prologue     group.f
c***date written       850601  
c***revision date      11/6/94      
c
c   7 march    1991    rlm at lanl
c      fixing bug in reduction of cn symmetry to something abelian.
c
c***keywords           
c***author             martin,richard and saxe,paul (lanl) 
c***source             @(#)group.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       group.f
      implicit integer (a-z)
c
      character*(*) grp,lirrep(*),gengrp,labop(*),axis
      character*4 subgrp,ptgrp
      character*4 grperr
      integer ctoi,naxis,lambda(*)
      logical prtsym,redsym,nodd
c
      common/io/inp,iout
 1000 format(5x,'utilizing group :',2x,a4)
c
c     ----- determine the generic group type -----
c           the parameter redsym specifies if we are instructed to run
c           in reduced symmetry.  this is primarily for the ci program
c           which is not set up to handle degenerate representations.
c           note that there is an additional parameter which is
c           sometimes modified in the reduction to the subgroup, and
c           that is "axis", the principal rotation axis.  for example,
c           in d3h, the principal c3 axis is aligned along z in the
c           full group, but this cannot be taken coincident with the
c           c2 principal axis in c2v.  in order to get around this,
c           certain reductions modify the flag which controls the
c           principal rotation axis known to the rest of the code.
c
      ptgrp=grp
      subgrp=grp
      if(redsym) then
c        check for some special groups.
         if     (ptgrp.eq.'c*v') then
            ptgrp='c2v'
         else if(ptgrp.eq.'d*h') then
            ptgrp='d2h'
         else if(ptgrp.eq.'t') then
            ptgrp='c2v'
         else if(ptgrp.eq.'th') then
            ptgrp='d2h'
         else if(ptgrp.eq.'td') then
            ptgrp='c2v'
         else if(ptgrp.eq.'o') then
            ptgrp='c2v'
         else if(ptgrp.eq.'oh') then
            ptgrp='d2h'
         else if(ptgrp.eq.'i') then
            ptgrp='c2'
         else if(ptgrp.eq.'ih') then
            ptgrp='c2h'
         else if(ptgrp.eq.'d2d') then
            ptgrp='c2v'
         else if(ptgrp.eq.'d3d') then
            ptgrp='c2h'
         else if(ptgrp.eq.'d4d') then
            ptgrp='c2v'
         else if(ptgrp.eq.'d5d') then
            ptgrp='c2'
         else if(ptgrp.eq.'d6d') then
            ptgrp='c2v'
         else if(ptgrp.eq.'d3') then
            ptgrp='c2'
         else if(ptgrp.eq.'d4') then
            ptgrp='c2'
         else if(ptgrp.eq.'d5') then
            ptgrp='c2'
         else if(ptgrp.eq.'d6') then
            ptgrp='c2v'
         else if(ptgrp.eq.'s4'.or.ptgrp.eq.'s8') then
            ptgrp='c2'
         else if(ptgrp.eq.'s6') then
            ptgrp='ci'
         else if(ptgrp.eq.'c1') then
         else if(ptgrp.eq.'cs') then
         else if(ptgrp.eq.'ci') then
         else
c           the only remaining groups should be cn,cnh,cnv,dnh.
            naxis=ctoi(ptgrp(2:2))
            nodd=(mod(naxis,2).ne.0)
            if(naxis.ge.3) then
               if(nodd) then
                  if(ptgrp(1:1).eq.'c') then
                     if(ptgrp(3:3).eq.'h'.or.ptgrp(3:3).eq.'v') then
c                       cnh,cnv reduce to cs if n is odd.
                        ptgrp='cs'
                     else
c                       cn reduces to c1 if n is odd.
                        ptgrp='c1'
                     endif
                  else if(ptgrp(1:1).eq.'d') then
c                    dnh reduces to c2v if n is odd.
c                    choose the principal axis to lie along y.
c                    if the principal axis is not the default, then
c                    the user has overridden the default, so do nothing.
                     ptgrp='c2v'
                     if (axis.eq.'z') then
                        axis='y'
                     end if
                  endif
               else
c                 cn,cnv,cnh,dnh reduce to c2,c2v,c2h,d2h if n is even. 
                  ptgrp(2:2)='2' 
               end if 
            end if
         end if
c
c        check to see that the reduced point group is something we recognize.
         if(    ptgrp.eq.'d2h'
     $      .or.ptgrp.eq.'d2'
     $      .or.ptgrp.eq.'c2v'
     $      .or.ptgrp.eq.'c2h'
     $      .or.ptgrp.eq.'c2'
     $      .or.ptgrp.eq.'cs'
     $      .or.ptgrp.eq.'ci'
     $      .or.ptgrp.eq.'c1') then
         else
            call lnkerr('error in point group reduction')
         end if
         subgrp=ptgrp
      end if
c
c
      if(prtsym) write(iout,1000) subgrp
c
c
    1 continue
      if (ptgrp.eq.'c1') then
         gengrp='c1'
         ngen=1
         nop=1
         nirrep=1
         lirrep(1)='a'
         labop(1)='     e'
      else if (ptgrp.eq.'cs') then
c
         gengrp='cs'
         ngen=1
         nop=2
         nirrep=2
         lirrep(1)='a '''
         lirrep(2)='a "'
         labop(1)='     e'
         labop(2)='  sigma.h'
      else if (ptgrp.eq.'ci') then
         gengrp='ci'
         ngen=1
         nop=2
         nirrep=2
         lirrep(1)='a g'
         lirrep(2)='a u'
         labop(1)='     e'
         labop(2)='     i'
      else if (ptgrp.eq.'c*v') then
c        for now, default to c2v group.
         gengrp='cnv'
         naxis=2
         ngen=2
         nop=4
         nirrep=4
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         labop(1)='      e'
         labop(2)='      c2'
         labop(3)=' sigma.v(xz)'
         labop(4)=' sigma.v(yz)'
      else if (ptgrp.eq.'d*h') then
c        for now, default to d2h group.
c        gengrp='d*h'
         gengrp='dnh'
         naxis=2
         ngen=3
         nop=8
         nirrep=8
         lirrep(1)='a g'
         lirrep(2)='b1g'
         lirrep(3)='b2g'
         lirrep(4)='b3g'
         lirrep(5)='a u'
         lirrep(6)='b1u'
         lirrep(7)='b2u'
         lirrep(8)='b3u'
      else if (ptgrp.eq.'t') then
         gengrp='t'
         ngen=0
         nop=12
         nirrep=3
         lirrep(1)='a'
         lirrep(2)='e'
         lirrep(3)='t'
         labop(1)='     e'
         labop(2)='   c3(1)'
         labop(3)='   c3(2)'
      else if (ptgrp.eq.'th') then
         gengrp='th'
         ngen=0
         nop=24
         nirrep=6
         lirrep(1)='a g'
         lirrep(2)='e g'
         lirrep(3)='t g'
         lirrep(4)='a u'
         lirrep(5)='e u'
         lirrep(6)='t u'
      else if (ptgrp.eq.'td') then
         gengrp='td'
         ngen=0
         nop=24
         nirrep=5
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e'
         lirrep(4)='t1'
         lirrep(5)='t2'
      else if (ptgrp.eq.'o') then
         gengrp='o'
         naxis=4
         ngen=1
         nop=24
         nirrep=5
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e'
         lirrep(4)='t1'
         lirrep(5)='t2'
      else if (ptgrp.eq.'oh') then
         gengrp='oh'
         naxis=4
         ngen=1
         nop=48
         nirrep=10
         lirrep(1)='a1g'
         lirrep(2)='a2g'
         lirrep(3)='e g'
         lirrep(4)='t1g'
         lirrep(5)='t2g'
         lirrep(6)='a1u'
         lirrep(7)='a2u'
         lirrep(8)='e u'
         lirrep(9)='t1u'
         lirrep(10)='t2u'
      else if (ptgrp.eq.'i') then
         gengrp='i'
         ngen=0
         nop=60
         nirrep=5
         lirrep(1)='a'
         lirrep(2)='t1'
         lirrep(3)='t2'
         lirrep(4)='g'
         lirrep(5)='h'
      else if (ptgrp.eq.'ih') then
         gengrp='ih'
         ngen=0
         nop=120
         nirrep=10
         lirrep(1)='a g'
         lirrep(2)='t1g'
         lirrep(3)='t2g'
         lirrep(4)='g g'
         lirrep(5)='h g'
         lirrep(6)='a u'
         lirrep(7)='t1u'
         lirrep(8)='t2u'
         lirrep(9)='g u'
         lirrep(10)='h u'
      else if (ptgrp(1:1).eq.'c') then
         if (len(ptgrp).lt.2) then
            gengrp=ptgrp
            grperr=ptgrp
            call lnkerr('bad point group symbol: '//grperr)
         end if
         naxis=ctoi(ptgrp(2:2))
         if (len(ptgrp).eq.2) then
            gengrp='cn'
         else
            if (ptgrp(3:3).eq.'h') then
               if (naxis.eq.1) then
                  ptgrp='cs'
                  go to 1
               end if
               gengrp='cnh'
            else if (ptgrp(3:3).eq.'v') then
               gengrp='cnv'
            else if (ptgrp(3:3).eq.' ') then
               gengrp='cn'
            else
               gengrp=ptgrp
               grperr=ptgrp
               call lnkerr('bad point group symbol: '//grperr)
            end if
         end if
c
         if (ptgrp.eq.'c2') then
            ngen=1
            nop=2
            nirrep=2
            lirrep(1)='a'
            lirrep(2)='b'
            labop(1)='     e'
            labop(2)='     c2'
         else if (ptgrp.eq.'c3') then
            ngen=1
            nop=3
            nirrep=2
            lirrep(1)='a'
            lirrep(2)='e'
            labop(1)='     e'
            labop(2)='    c3*1'
            labop(3)='    c3*2'
         else if (ptgrp.eq.'c4') then
            ngen=1
            nop=4
            nirrep=3
            lirrep(1)='a'
            lirrep(2)='b'
            lirrep(3)='e'
            labop(1)='     e'
            labop(2)='    c4*1'
            labop(3)='     c2'
            labop(4)='    c4*3'
         else if (ptgrp.eq.'c5') then
            ngen=1
            nop=5
            nirrep=3
            lirrep(1)='a'
            lirrep(2)='e1'
            lirrep(3)='e2'
            labop(1)='     e'
            labop(2)='    c5*1'
            labop(3)='    c5*2'
            labop(4)='    c5*3'
            labop(5)='    c5*4'
cdir$ block
         else if (ptgrp.eq.'c6') then
            ngen=1
            nop=6
            nirrep=4
            lirrep(1)='a'
            lirrep(2)='b'
            lirrep(3)='e1'
            lirrep(4)='e2'
            labop(1)='      e'
            labop(2)='      c6'
            labop(3)='      c3'
            labop(4)='      c2'
            labop(5)='     c3*2'
            labop(6)='     c6*5'
         else if (ptgrp.eq.'c7') then
            ngen=1
            nop=7
            nirrep=4
            lirrep(1)='a'
            lirrep(2)='e1'
            lirrep(3)='e2'
            lirrep(4)='e3'
            labop(1)='     e'
            labop(2)='     c7'
            labop(3)='    c7*2'
            labop(4)='    c7*3'
            labop(5)='    c7*4'
            labop(6)='    c7*5'
            labop(7)='    c7*6'
         else if (ptgrp.eq.'c8') then
            ngen=1
            nop=8
            nirrep=5
            lirrep(1)='a'
            lirrep(2)='b'
            lirrep(3)='e1'
            lirrep(4)='e2'
            lirrep(5)='e3'
            labop(1)='     e'
            labop(2)='     c8'
            labop(3)='     c4'
            labop(4)='    c8*3'
            labop(5)='     c2'
            labop(6)='    c8*5'
            labop(7)='    c4*3'
            labop(8)='    c8*7'
         else if (ptgrp.eq.'c2v') then
            ngen=2
            nop=4
            nirrep=4
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            labop(1)='      e'
            labop(2)='      c2'
            labop(3)=' sigma.v(xz)'
            labop(4)=' sigma.v(yz)'
         else if (ptgrp.eq.'c3v') then
            ngen=2
            nop=6
            nirrep=3
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='e'
            labop(1)='      e'
            labop(2)='      c3'
            labop(3)='     c3*2'
            labop(4)='   sigma.v(1)'
            labop(5)='   sigma.v(2)'
            labop(6)='   sigma.v(3)'
         else if (ptgrp.eq.'c4v') then
            ngen=2
            nop=8
            nirrep=5
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e'
         else if (ptgrp.eq.'c5v') then
            ngen=2
            nop=10
            nirrep=4
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='e1'
            lirrep(4)='e2'
         else if (ptgrp.eq.'c6v') then
            ngen=2
            nop=12
            nirrep=6
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e1'
            lirrep(6)='e2'
         else if (ptgrp.eq.'c2h') then
            ngen=2
            nop=4
            nirrep=4
            lirrep(1)='a g'
            lirrep(2)='b g'
            lirrep(3)='a u'
            lirrep(4)='b u'
         else if (ptgrp.eq.'c3h') then
            ngen=2
            nop=6
            nirrep=4
            lirrep(1)='a '''
            lirrep(2)='a "'
            lirrep(3)='e '''
            lirrep(4)='e "'
         else if (ptgrp.eq.'c4h') then
            ngen=2
            nop=8
            nirrep=6
            lirrep(1)='a g'
            lirrep(2)='b g'
            lirrep(3)='a u'
            lirrep(4)='b u'
            lirrep(5)='e g'
            lirrep(6)='e u'
         else if (ptgrp.eq.'c5h') then
            ngen=2
            nop=10
            nirrep=6
            lirrep(1)='a '''
            lirrep(2)='a "'
            lirrep(3)='e1'''
            lirrep(4)='e2'''
            lirrep(5)='e1"'
            lirrep(6)='e2"'
         else if (ptgrp.eq.'c6h') then
            ngen=2
            nop=12
            nirrep=8
            lirrep(1)='a g'
            lirrep(2)='b g'
            lirrep(3)='a u'
            lirrep(4)='b u'
            lirrep(5)='e1g'
            lirrep(6)='e2g'
            lirrep(7)='e1u'
            lirrep(8)='e2u'
         else
            grperr=ptgrp
            call lnkerr('cannot handle point group '//grperr)
         end if
cdir$ block
      else if (ptgrp(1:1).eq.'s') then
         if (len(ptgrp).lt.2) then
            gengrp=ptgrp
            grperr=ptgrp
            call lnkerr('bad point group symbol: '//grperr)
         end if
         naxis=ctoi(ptgrp(2:2))
         if (naxis.eq.2) then
            ptgrp='ci'
            go to 1
         end if
         gengrp='sn'
         if (ptgrp.eq.'s4') then
            ngen=1
            nop=4
            nirrep=3
            lirrep(1)='a'
            lirrep(2)='b'
            lirrep(3)='e'
         else if (ptgrp.eq.'s6') then
            ngen=1
            nop=6
            nirrep=4
            lirrep(1)='a g'
            lirrep(2)='a u'
            lirrep(3)='e g'
            lirrep(4)='e u'
         else if (ptgrp.eq.'s8') then
            ngen=1
            nop=8
            nirrep=5
            lirrep(1)='a'
            lirrep(2)='b'
            lirrep(3)='e1'
            lirrep(4)='e2'
            lirrep(5)='e3'
         else
            grperr=ptgrp
            call lnkerr('cannot handle point group '//grperr)
         end if
      else if (ptgrp(1:1).eq.'d') then
         if (len(ptgrp).lt.2) then
            gengrp=ptgrp
            grperr=ptgrp
            call lnkerr('bad point group symbol: '//grperr)
         end if
         naxis=ctoi(ptgrp(2:2))
         if (len(ptgrp).eq.2) then
            gengrp='dn'
         else
            if (ptgrp(3:3).eq.'h') then
               gengrp='dnh'
            else if (ptgrp(3:3).eq.'d') then
               gengrp='dnd'
            else if (ptgrp(3:3).eq.' ') then
               gengrp='dn'
            else
               gengrp=ptgrp
               grperr=ptgrp
               call lnkerr('bad point group symbol: '//grperr)
            end if
         end if
c
         if (ptgrp.eq.'d2') then
            ngen=2
            nop=4
            nirrep=4
            lirrep(1)='a'
            lirrep(2)='b1'
            lirrep(3)='b2'
            lirrep(4)='b3'
         else if (ptgrp.eq.'d3') then
            ngen=2
            nop=6
            nirrep=3
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='e'
         else if (ptgrp.eq.'d4') then
            ngen=2
            nop=8
            nirrep=5
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e'
         else if (ptgrp.eq.'d5') then
            ngen=2
            nop=10
            nirrep=4
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='e1'
            lirrep(4)='e2'
         else if (ptgrp.eq.'d6') then
            ngen=2
            nop=12
            nirrep=6
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e1'
            lirrep(6)='e2'
         else if (ptgrp.eq.'d2h') then
            ngen=3
            nop=8
            nirrep=8
            lirrep(1)='a g'
            lirrep(2)='b1g'
            lirrep(3)='b2g'
            lirrep(4)='b3g'
            lirrep(5)='a u'
            lirrep(6)='b1u'
            lirrep(7)='b2u'
            lirrep(8)='b3u'
         else if (ptgrp.eq.'d3h') then
            ngen=3
            nop=12
            nirrep=6
            lirrep(1)='a1'''
            lirrep(2)='a2'''
            lirrep(3)='a1"'
            lirrep(4)='a2"'
            lirrep(5)='e '''
            lirrep(6)='e "'
            labop(1)='      e'
            labop(2)='      c3'
            labop(3)='     c3*2'
            labop(4)='      c2'
            labop(5)='      c2'''
            labop(6)='      c2"'
            labop(7)='   sigma.h(xy)'
            labop(8)='      s3'
            labop(9)='     s3*2'
            labop(10)='   sigma.v'
            labop(11)='   sigma.v'''
            labop(12)='   sigma.v"'
         else if (ptgrp.eq.'d4h') then
            ngen=3
            nop=16
            nirrep=10
            lirrep(1)='a1g'
            lirrep(2)='a2g'
            lirrep(3)='b1g'
            lirrep(4)='b2g'
            lirrep(5)='a1u'
            lirrep(6)='a2u'
            lirrep(7)='b1u'
            lirrep(8)='b2u'
            lirrep(9)='e g'
            lirrep(10)='e u'
         else if (ptgrp.eq.'d5h') then
            ngen=3
            nop=20
            nirrep=8
            lirrep(1)='a1'''
            lirrep(2)='a2'''
            lirrep(3)='a1"'
            lirrep(4)='a2"'
            lirrep(5)='e1'''
            lirrep(6)='e2'''
            lirrep(7)='e1"'
            lirrep(8)='e2"'
         else if (ptgrp.eq.'d6h') then
            ngen=3
            nop=24
            nirrep=12
            lirrep(1)='a1g'
            lirrep(2)='a2g'
            lirrep(3)='b1g'
            lirrep(4)='b2g'
            lirrep(5)='a1u'
            lirrep(6)='a2u'
            lirrep(7)='b1u'
            lirrep(8)='b2u'
            lirrep(9)='e1g'
            lirrep(10)='e2g'
            lirrep(11)='e1u'
            lirrep(12)='e2u'
         else if (ptgrp.eq.'d8h') then
            ngen=3
            nop=32
            nirrep=14
            lirrep(1)='a1g'
            lirrep(2)='a2g'
            lirrep(3)='b1g'
            lirrep(4)='b2g'
            lirrep(5)='a1u'
            lirrep(6)='a2u'
            lirrep(7)='b1u'
            lirrep(8)='b2u'
            lirrep(9)='e1g'
            lirrep(10)='e2g'
            lirrep(11)='e3g'
            lirrep(12)='e1u'
            lirrep(13)='e2u'
            lirrep(14)='e3u'
cdir$ block
         else if (ptgrp.eq.'d2d') then
            ngen=2
            nop=8
            nirrep=5
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e'
         else if (ptgrp.eq.'d3d') then
            ngen=3
            nop=12
            nirrep=6
            lirrep(1)='a1g'
            lirrep(2)='a2g'
            lirrep(3)='a1u'
            lirrep(4)='a2u'
            lirrep(5)='e g'
            lirrep(6)='e u'
         else if (ptgrp.eq.'d4d') then
            ngen=2
            nop=16
            nirrep=7
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e1'
            lirrep(6)='e2'
            lirrep(7)='e3'
         else if (ptgrp.eq.'d5d') then
            ngen=3
            nop=20
            nirrep=8
            lirrep(1)='a1g'
            lirrep(2)='a2g'
            lirrep(3)='a1u'
            lirrep(4)='a2u'
            lirrep(5)='e1g'
            lirrep(6)='e2g'
            lirrep(7)='e1u'
            lirrep(8)='e2u'
         else if (ptgrp.eq.'d6d') then
            ngen=2
            nop=24
            nirrep=9
            lirrep(1)='a1'
            lirrep(2)='a2'
            lirrep(3)='b1'
            lirrep(4)='b2'
            lirrep(5)='e1'
            lirrep(6)='e2'
            lirrep(7)='e3'
            lirrep(8)='e4'
            lirrep(9)='e5'
         else
            grperr=ptgrp
            call lnkerr('cannot handle point group '//grperr)
         end if
      else
         gengrp=ptgrp
         grperr=ptgrp
         call lnkerr('bad point group symbol: '//grperr)
      end if
c
c     ----- tabulate degeneracies of irreps -----
c
      do 20 irrep=1,nirrep
         if (lirrep(irrep)(1:1).eq.'a') then
            lambda(irrep)=1
         else if (lirrep(irrep)(1:1).eq.'b') then
            lambda(irrep)=1
         else if (lirrep(irrep)(1:1).eq.'e') then
            lambda(irrep)=2
         else if (lirrep(irrep)(1:1).eq.'t') then
            lambda(irrep)=3
         else if (lirrep(irrep)(1:1).eq.'g') then
            lambda(irrep)=4
         else if (lirrep(irrep)(1:1).eq.'h') then
            lambda(irrep)=5
         else
            call lnkerr('error with degeneracy of irrep')
         end if
 20   continue
c
c
      return
      end
