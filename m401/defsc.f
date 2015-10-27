*deck @(#)defsc.f	5.1  11/6/94
      function defsc(ia,j)
c***begin prologue     defsc
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           scale factor
c***author             martin, richard (lanl)
c***source             @(#)defsc.f	5.1   11/6/94
c***purpose            returns a default scale factor for a specific
c                      orbital and atom.
c***description
c     defsc is a real function used as:
c       scalef=defsc(ia,j)
c         ia      atomic number.
c         j       orbital type.
c***references
c***routines called    (none)
c***end prologue       defsc
c
      implicit real*8(a-h,o-z)
c
      real*8 defsc
c
c
      dimension x(19), y(19), n(19)
      data zero/0.0d0/, one/1.0d0/, f54/54.0d0/
      data n/1,2,2,3,3,4,4,3,5,5,4,6,6,5,4,7,7,6,5/
      data x/19*0.0d0/
      data y/0.017d0,0.097d0,0.075d0,0.279d0,0.266d0,0.662d0,0.678d0,
     $       0.224d0,1.648d0,1.865d0,0.745d0,2.679d0,3.024d0,1.465d0,
     $       0.510d0,4*0.0d0/
      save zero,one,f54,n
      save x,y
c
      defsc = zero
      if(j.eq.0.or.iabs(j).gt.15.or.ia.lt.1) return
      defsc = one / float(ia)
      if(j.eq.1) return
      rn = float(n(iabs(j)))
      rx = x(iabs(j))
      ry = y(iabs(j)) * f54
      ria = float(ia)
      defsc = rn*(rn-1)/(rx+(ry/ria))
c
c
      return
      end
