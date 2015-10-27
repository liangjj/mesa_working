*deck @(#)tmat.f	5.1  11/6/94
      subroutine tmat(t,nop,gengrp,naxis,axis)
c
      implicit integer (a-z)
c
      real*8 t(3,3,nop)
      character*(*) gengrp,axis
c
c     ----- form the transformation matrices for the coordinates
c     under each operation. the principal rotation axis is given by "axis".
c
      call runit(t,3)
      if (gengrp.eq.'c1') then
      else if (gengrp.eq.'cs') then
         call sigmah(t(1,1,2),1,axis,t)
      else if (gengrp.eq.'ci') then
         call ci(t(1,1,2),1,t)
      else if (gengrp.eq.'cn') then
         call cnaxis(t,naxis,axis)
      else if (gengrp.eq.'sn') then
         call snaxis(t,naxis)
      else if (gengrp.eq.'cnh') then
         call cnaxis(t,naxis,axis)
         if(mod(naxis,2).eq.0) then
c           naxis even..
            call ci(t(1,1,naxis+1),naxis,t)
         else
c           naxis odd..	
            call sigmah(t(1,1,naxis+1),naxis,axis,t)
         endif
      else if (gengrp.eq.'cnv') then
         call cnaxis(t,naxis,axis)
         call sigmav(t(1,1,naxis+1),naxis,axis,t)
      else if (gengrp.eq.'dn') then
         call cnaxis(t,naxis,axis)
         call c2prim(t(1,1,naxis+1),naxis,axis,t)
      else if (gengrp.eq.'dnh') then
         call cnaxis(t,naxis,axis)
         call c2prim(t(1,1,naxis+1),naxis,axis,t)
         if(mod(naxis,2).eq.0) then
c           naxis even..
            call ci(t(1,1,2*naxis+1),2*naxis,t)
         else
c           naxis odd..
            call sigmah(t(1,1,2*naxis+1),2*naxis,axis,t)
         endif
      else if (gengrp.eq.'dnd') then
         if(mod(naxis,2).eq.0) then
c           naxis even..           
            call snaxis(t,2*naxis)
            call c2prim(t(1,1,2*naxis+1),2*naxis,axis,t) 
         else
c           naxis odd..
            call cnaxis(t,naxis,axis)
            call c2prim(t(1,1,naxis+1),naxis,axis,t)
            call ci(t(1,1,2*naxis+1),2*naxis,t)
         endif
      else if (gengrp.eq.'o') then
c        fills the 24 operations of the octahedral group
         call octa(t)
      else if (gengrp.eq.'oh') then
         call octa(t)
c        oh = o x i
         call ci(t(1,1,25),24,t)
      endif
c
c
      return
      end
