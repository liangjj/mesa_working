*deck @(#)newton.f	1.1  4/25/95
      subroutine newton(ian,nrings,rpts,rwts,edges,nr,mxrings,
     $                  nrtot,nwtot)
c***begin prologue     newton.f
c***date written       940304   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard(lanl)
c***source             @(#)newton.f	1.1   4/25/95
c***purpose            
c***description
c   for each atom, we define a number of rings which cover the integration
c   region from the origin to some finite distance rmax.
c   within each ring, we are free to define the number of points to
c   use in the quadrature (the effective step size), and the order of the
c   newton-cotes quadrature to use in that interval.
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       newton.f
      implicit none
c     --- input variables -----
      integer ian,mxrings
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      real*8 edges(0:mxrings)
c     --- output arrays ---
      integer nr(mxrings)
c     --- output variables ---
      integer nrings,nrtot,nwtot
      real*8 rpts(*),rwts(*)
c     --- scratch arrays ---
c     --- local variables ---
      integer ring
      integer ptr,ptwt
c
c     --- get the number of rings to use for this atom
      nrings=8
c
c     --- set up the edges and number of points in each ring
      edges(0)=1.0d-10
      do 10 ring=1,nrings
         edges(ring)=edges(ring-1)+0.4d+00
         nr(ring)=5
   10 continue
c
c     last ring
      nrings=nrings+1
      edges(nrings)=edges(nrings-1)+150.d00
      nr(nrings)=5
c
c
      ptr=1
      ptwt=1
      nrtot=0
      nwtot=0
      do 20 ring=1,nrings
         call necote(edges(ring-1),edges(ring),rpts(ptr),rwts(ptwt),
     $               nr(ring),.false.)
c        --- the biasing by -1 is to avoid including the last point
c            twice.  for right now, we have a little sleaze going on here.
c            in order to do the newton-cotes quadrature over the intervals
c            it is helpful to have points duplicated -- i.e. [a,b], [b,c],
c            the common endpoints in each region would suggest that we
c            include radial points at the edges twice. however, this
c            means more function evaluations, potential problems in the
c            voronoi weights, etc.  it is simple to bias the pointer to
c            the radial points by one so as to avoid the overcounting
c            and do similar tricks in the actual quadrature to use
c            the common endpoints. however, the weights are a little
c            awkward. since each ring may have a different number of points
c            with a different order newton-cotes rule, the weights for
c            point b in ring 1 may be different from the weights used
c            for point b in ring 2.  for now, we get around this by 
c            duplicating the weights.
         ptr=ptr+nr(ring)-1
         ptwt=ptwt+nr(ring)*(nr(ring)-1)
         nrtot=nrtot+nr(ring)-1
         nwtot=nwtot+nr(ring)*(nr(ring)-1)
   20 continue
c     --- last point is not shared, account for this.
      nrtot=nrtot+1
c
c
      return
      end
