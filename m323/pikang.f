*deck @(#)pikang.f	5.1  11/6/94
      function pikang(angmom,symcen,cenang,nsymcn,nbtype,cencod,angcod)
c
c***module to step in a canonical fashion through the angular momentum
c   shells of the four centres in symcen. the angular momenta are stored
c   in angmom. this function returns a value of .false. until it has
c   completed the list, when it returns .true.. initialisation is signaled
c   by a negative value of angcod.
c
c         paul saxe          28 june 1984                   lanl
c
      implicit integer (a-z)
c
      integer angmom(4),symcen(4),cenang(nsymcn,nbtype)
      logical pikang,nxtang,temp
c
c     ----- check if initialising -----
c
      if (angcod.lt.0) then
         pikang=.true.
         if (nxtang(angmom,symcen,cenang,nsymcn,nbtype,1,0)) return
         if (nxtang(angmom,symcen,cenang,nsymcn,nbtype,2,0)) return
         if (nxtang(angmom,symcen,cenang,nsymcn,nbtype,3,0)) return
         if (nxtang(angmom,symcen,cenang,nsymcn,nbtype,4,0)) return
         angcod=0
         pikang=.false.
      else
         pikang=.false.
         temp=nxtang(angmom,symcen,cenang,nsymcn,nbtype,4,1)
         if (temp.or.
     #      ((symcen(1).eq.symcen(3).and.angmom(1).eq.angmom(3)).and.
     #      ((symcen(2).eq.symcen(4).and.angmom(4).gt.angmom(2))
     #        .or.symcen(4).gt.symcen(2))).or.
     #       (symcen(3).eq.symcen(4).and.angmom(4).gt.angmom(3)))then
            temp=nxtang(angmom,symcen,cenang,nsymcn,nbtype,4,0)
            temp=nxtang(angmom,symcen,cenang,nsymcn,nbtype,3,1)
            if (temp.or.(symcen(1).eq.symcen(3).and.angmom(3).gt.
     #                                              angmom(1))) then
               temp=nxtang(angmom,symcen,cenang,nsymcn,nbtype,3,0)
               temp=nxtang(angmom,symcen,cenang,nsymcn,nbtype,2,1)
               if (temp.or.(symcen(1).eq.symcen(2).and.
     #                      angmom(2).gt.angmom(1))) then
                  temp=nxtang(angmom,symcen,cenang,nsymcn,nbtype,2,0)
                  pikang=nxtang(angmom,symcen,cenang,nsymcn,nbtype,1,1)
               end if
            end if
         end if
      end if
c
c
      return
      end
