*deck @(#)grpsym.f
      subroutine grpsym(gengrp,lirrep,maxrep)
c
c***begin prologue     grpsym.f
c***date written       850601  
c***revision date      11/6/94      
c
c***keywords           
c***author             Schneider, Barry(NSF) 
c***source             sym
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       grpsym.f
      implicit none
c
      integer inp, iout, maxrep
      character*(*) gengrp, lirrep(maxrep)
      common/io/inp,iout
      write(iout,1000) gengrp
 1000 format(5x,'utilizing group :',2x,a4)
      if(gengrp.eq.'cs') then
           lirrep(1)='a'''
           lirrep(2)='a"'
      elseif(gengrp.eq.'c1') then
           lirrep(1)='a'
           lirrep(2)='e'
      elseif(gengrp.eq.'c2') then
           lirrep(1)='a'
           lirrep(2)='b'
      elseif(gengrp.eq.'c3') then
           lirrep(1)='a'
           lirrep(2)='e'
      else if (gengrp.eq.'c4') then
         lirrep(1)='a'
         lirrep(2)='b'
         lirrep(3)='e'
      else if (gengrp.eq.'c5') then
         lirrep(1)='a'
         lirrep(2)='e1'
         lirrep(3)='e2'
      else if (gengrp.eq.'c6') then
         lirrep(1)='a'
         lirrep(2)='b'
         lirrep(3)='e1'
         lirrep(4)='e2'
      else if (gengrp.eq.'c7') then
         lirrep(1)='a'
         lirrep(2)='e1'
         lirrep(3)='e2'
         lirrep(4)='e3'
      else if (gengrp.eq.'c8') then
         lirrep(1)='a'
         lirrep(2)='b'
         lirrep(3)='e1'
         lirrep(4)='e2'
         lirrep(5)='e3'
      elseif(gengrp.eq.'ci') then
           lirrep(1)='ag'
           lirrep(2)='au'
      elseif(gengrp.eq.'t') then
         lirrep(1)='a'
         lirrep(2)='e'
         lirrep(3)='t'
      else if (gengrp.eq.'th') then
         lirrep(1)='ag'
         lirrep(2)='eg'
         lirrep(3)='tg'
         lirrep(4)='au'
         lirrep(5)='eu'
         lirrep(6)='tu'
      else if (gengrp.eq.'td') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e'
         lirrep(4)='t1'
         lirrep(5)='t2'
      else if (gengrp.eq.'o') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e'
         lirrep(4)='t1'
         lirrep(5)='t2'
      else if (gengrp.eq.'oh') then
         lirrep(1)='a1g'
         lirrep(2)='a2g'
         lirrep(3)='eg'
         lirrep(4)='t1g'
         lirrep(5)='t2g'
         lirrep(6)='a1u'
         lirrep(7)='a2u'
         lirrep(8)='eu'
         lirrep(9)='t1u'
         lirrep(10)='t2u'
      else if (gengrp.eq.'i') then
         lirrep(1)='a'
         lirrep(2)='t1'
         lirrep(3)='t2'
         lirrep(4)='g'
         lirrep(5)='h'
      else if (gengrp.eq.'ih') then
         lirrep(1)='ag'
         lirrep(2)='t1g'
         lirrep(3)='t2g'
         lirrep(4)='gg'
         lirrep(5)='hg'
         lirrep(6)='au'
         lirrep(7)='t1u'
         lirrep(8)='t2u'
         lirrep(9)='gu'
         lirrep(10)='hu'
      else if (gengrp.eq.'c2v') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
      else if (gengrp.eq.'c3v') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e'
      else if (gengrp.eq.'c4v') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         lirrep(5)='e'
      else if (gengrp.eq.'c5v') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e1'
         lirrep(4)='e2'
      else if (gengrp.eq.'c6v') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         lirrep(5)='e1'
         lirrep(6)='e2'
      else if (gengrp.eq.'c2h') then
         lirrep(1)='ag'
         lirrep(2)='bg'
         lirrep(3)='au'
         lirrep(4)='bu'
      else if (gengrp.eq.'c3h') then
         lirrep(1)='a'''
         lirrep(2)='a"'
         lirrep(3)='e'''
         lirrep(4)='e"'
      else if (gengrp.eq.'c4h') then
         lirrep(1)='ag'
         lirrep(2)='bg'
         lirrep(3)='au'
         lirrep(4)='bu'
         lirrep(5)='eg'
         lirrep(6)='eu'
      else if (gengrp.eq.'c5h') then
         lirrep(1)='a'''
         lirrep(2)='a"'
         lirrep(3)='e1'''
         lirrep(4)='e2'''
         lirrep(5)='e1"'
         lirrep(6)='e2"'
      else if (gengrp.eq.'c6h') then
         lirrep(1)='ag'
         lirrep(2)='bg'
         lirrep(3)='au'
         lirrep(4)='bu'
         lirrep(5)='e1g'
         lirrep(6)='e2g'
         lirrep(7)='e1u'
         lirrep(8)='e2u'
      elseif (gengrp.eq.'s4') then
         lirrep(1)='a'
         lirrep(2)='b'
         lirrep(3)='e'
      else if (gengrp.eq.'s6') then
         lirrep(1)='ag'
         lirrep(2)='au'
         lirrep(3)='eg'
         lirrep(4)='eu'
      else if (gengrp.eq.'s8') then
         lirrep(1)='a'
         lirrep(2)='b'
         lirrep(3)='e1'
         lirrep(4)='e2'
         lirrep(5)='e3'
      else if(gengrp.eq.'d2') then
         lirrep(1)='a'
         lirrep(2)='b1'
         lirrep(3)='b2'
         lirrep(4)='b3'
      else if (gengrp.eq.'d3') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e'
      else if (gengrp.eq.'d4') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         lirrep(5)='e'
      else if (gengrp.eq.'d5') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='e1'
         lirrep(4)='e2'
      else if (gengrp.eq.'d6') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         lirrep(5)='e1'
         lirrep(6)='e2'
      else if (gengrp.eq.'d2h') then
         lirrep(1)='ag'
         lirrep(2)='b1g'
         lirrep(3)='b2g'
         lirrep(4)='b3g'
         lirrep(5)='au'
         lirrep(6)='b1u'
         lirrep(7)='b2u'
         lirrep(8)='b3u'
      else if (gengrp.eq.'d3h') then
         lirrep(1)='a1'''
         lirrep(2)='a2'''
         lirrep(3)='a1"'
         lirrep(4)='a2"'
         lirrep(5)='e'''
         lirrep(6)='e"'
      else if (gengrp.eq.'d4h') then
         lirrep(1)='a1g'
         lirrep(2)='a2g'
         lirrep(3)='b1g'
         lirrep(4)='b2g'
         lirrep(5)='a1u'
         lirrep(6)='a2u'
         lirrep(7)='b1u'
         lirrep(8)='b2u'
         lirrep(9)='eg'
         lirrep(10)='eu'
      else if (gengrp.eq.'d5h') then
         lirrep(1)='a1'''
         lirrep(2)='a2'''
         lirrep(3)='a1"'
         lirrep(4)='a2"'
         lirrep(5)='e1'''
         lirrep(6)='e2'''
         lirrep(7)='e1"'
         lirrep(8)='e2"'
      else if (gengrp.eq.'d6h') then
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
      else if (gengrp.eq.'d8h') then
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
      else if (gengrp.eq.'d2d') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         lirrep(5)='e'
      else if (gengrp.eq.'d3d') then
         lirrep(1)='a1g'
         lirrep(2)='a2g'
         lirrep(3)='a1u'
         lirrep(4)='a2u'
         lirrep(5)='eg'
         lirrep(6)='eu'
      else if (gengrp.eq.'d4d') then
         lirrep(1)='a1'
         lirrep(2)='a2'
         lirrep(3)='b1'
         lirrep(4)='b2'
         lirrep(5)='e1'
         lirrep(6)='e2'
         lirrep(7)='e3'
      else if (gengrp.eq.'d5d') then
         lirrep(1)='a1g'
         lirrep(2)='a2g'
         lirrep(3)='a1u'
         lirrep(4)='a2u'
         lirrep(5)='e1g'
         lirrep(6)='e2g'
         lirrep(7)='e1u'
         lirrep(8)='e2u'
      else if (gengrp.eq.'d6d') then
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
         call lnkerr('error in point group')
      endif
c
c
      return
      end
