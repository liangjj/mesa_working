*deck @(#)funcnm.f	5.1  11/6/94
      function funcnm(powx,powy,powz)
c***begin prologue     funcnm
c***date written       861001  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           cartesian, names
c***author             martin, richard (lanl)
c***source
c***purpose            a character function which returns the name of a
c                      cartesian polynomial.
c***description
c                      name=funcnm(powx,powy,powz)
c                        powx      the power of x.
c                        powy      the power of y.
c                        powz      the power of z.
c
c***references
c***routines called    none
c***end prologue       funcnm
      implicit integer (a-z)
      character*(*) funcnm
c
c     convert the powers of a polynomial term into a label.
      funcnm=' '
      pos=1
      if(powx.ne.0) then
         do 10 i=pos,pos+powx-1
            funcnm(i:i)='x'
   10    continue
         pos=pos+powx
      endif
      if(powy.ne.0) then
         do 20 i=pos,pos+powy-1
            funcnm(i:i)='y'
   20    continue
         pos=pos+powy
      endif
      if(powz.ne.0) then
         do 30 i=pos,pos+powz-1
            funcnm(i:i)='z'
   30    continue
         pos=pos+powz
      endif
c
c
      return
      end
