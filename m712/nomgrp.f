*deck @(#)nomgrp.f	5.1  11/6/94
      subroutine nomgrp(noprim,mgrpno,nat,nbtype,nmgrp)
      implicit integer (a-z)
c
      integer noprim(nat,nbtype),mgrpno(nat,nbtype)
c
c     ----- give each momentum group present a unique number from 1 to nmgrp
c
      nmgrp=0
      do 2 atom=1,nat
         do 1 type=1,nbtype
            if (noprim(atom,type).gt.0) then
               nmgrp=nmgrp+1
               mgrpno(atom,type)=nmgrp
            else
               mgrpno(atom,type)=0
            end if
    1    continue
    2 continue
c
c
      return
      end
