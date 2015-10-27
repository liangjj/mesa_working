*deck @(#)dordr.f	5.1  11/6/94
      subroutine dordr(dmo,d,nbf,norbs,bf2orb,trans)
      implicit integer(a-z)
      real*8 dmo(nbf,nbf),d(norbs,norbs)
      real*8 zero,two
      integer bf2orb(nbf)
      logical trans
      parameter (zero=0.0d+00,two=2.0d+00)
c
c     routine to reintroduce the frozen-core/frozen-core contribution
c     to the density matrix, and reorder from the drt list to the
c     original transformation matrix sequence.
c
c
      call rzero(dmo,nbf*nbf)
c
c     introduce the frozen core/frozen core block.
      do 10 mo=1,nbf
         if(bf2orb(mo).lt.0) then
            if(trans) then
               dmo(mo,mo)=zero
            else
               dmo(mo,mo)=two
            endif
         endif
   10 continue
c
c     the frozen core/drt space and frozen core/frozen virtual space
c     is already zero.
c     transpose the drt/drt space.
      do 30 imo=1,nbf
         do 20 jmo=1,nbf
            if(bf2orb(imo).gt.0.and.bf2orb(jmo).gt.0) then
               dmo(imo,jmo)=d(bf2orb(imo),bf2orb(jmo))
            endif
   20    continue
   30 continue
c
c     the drt/frozen virtual and frozen virtual/frozen virtual blocks
c     are already zero.
c
      return
      end
