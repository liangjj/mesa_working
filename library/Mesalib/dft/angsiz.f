*deck @(#)angsiz.f	5.3 4/17/95
      function angsiz(l)
c***begin prologue     angsiz.f
c***date written       950120  
c***revision date      4/17/95      
c
c***keywords           spherical harmonics, quadrature
c***author             martin, richard (lanl) 
c***source             @(#)angsiz.f	5.3   4/17/95
c***purpose            returns the number of angular points associated
c                      with quadratures on a sphere.
c***description
c                      this function returns the number of angular
c                      points associated with quadratures on the unit
c                      sphere designed to integrate exactly all spherical 
c                      harmonics through angular momentum l. 
c
c                      in the cases in which lebedev quadratures are
c                      available, they are used. note that for l=5,
c                      the icosahedral formula is returned.
c
c                      if a lebedev quadrature is not available, the 
c                      routine dies.
c
c                              the correspondence between order and npts is:
c                              l     npts
c                              ----------
c                              2        4    tetrahedral
c                              3        6    octahedral
c                              5       12    icosahedral formula
c                              5       16    abramowitz and stegun
c                              7       26    abramowitz and stegun
c                              9       38    lebedev
c                             11       50    lebedev
c                             13       74    lebedev
c                             15       86    lebedev
c                             17      110    lebedev
c                             19      146    lebedev
c                             23      194    lebedev
c                             29      302    lebedev
c                             35      434    lebedev
c                             41      590    lebedev
c                             47      852    lobatto
c***references
c 
c
c***routines called
c
c***end prologue       angsiz.f
      implicit none
      integer angsiz
c     --- input variables -----
      integer l
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer angpts(47)
      data angpts/0,4,6,0,12,0,26,0,38,0,50,0,74,0,86,0,110,0,146,0,
     $            0,0,194,0,0,0,0,0,302,0,0,0,0,0,434,5*0,590,5*0,852/
c
      save angpts
c
c
      if (l.gt.41) then
         call lnkerr('requested angular quadrature not available')
      endif
c
      angsiz=angpts(l)
c
      if (angsiz.eq.0) then
         call lnkerr('requested angular quadrature not available')
      endif
c
c
      return
      end
