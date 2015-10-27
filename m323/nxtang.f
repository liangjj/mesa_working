*deck @(#)nxtang.f	5.1  11/6/94
      function nxtang(angmom,symcen,cenang,nsymcn,nbtype,centre,
     #                switch)
c
c***module to either find the lowest angular momentum function
c   on symcen(centre) or to find the next angular momentum function
c   on symcen(centre), depending on whether switch is 0 or 1.
c   the function returns a .false. if a function is found or a .true.
c   if there are no higher functions on the centre. the value of the
c   angular momentum is returned in angmom(centre).
c
c      paul saxe         28 june 1984                  lanl
c
      implicit integer (a-z)
c
      integer angmom(4),symcen(4),cenang(nsymcn,nbtype)
      logical nxtang
c
c
c
      nxtang=.false.
c
c     ----- check switch to see if incrementing or starting -----
c
      if (switch.eq.0) then
         min=1
      else
         min=angmom(centre)+1
      end if
c
      do 1 angm=min,nbtype
         if (cenang(symcen(centre),angm).ne.0) then
            angmom(centre)=angm
            return
         end if
    1 continue
c
c     ----- have run out of functions on this centre -----
c
      nxtang=.true.
c
      return
      end
