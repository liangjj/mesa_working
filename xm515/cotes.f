*deck @(#)cotes.f	1.1  4/25/95
      subroutine cotes(alpha,npts,rule,nrings,rpts,rwts,jacob,
     $                  nr,mxrings,nrtot,nwtot)
c***begin prologue     cotes.f
c***date written       940304   (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard(lanl)
c***source             @(#)cotes.f	1.1   4/25/95
c***purpose            
c***description
c   the plan is as follows:
c      transform the radial variable from the range r(0,Inf) to x(0,1)
c         via r=alpha*x/1-x            if mu=1
c             r=alpha*x**2/(1-x)**2    if mu=2
c             r=alpha*x**3/(1-x)**3    if mu=3
c      where alpha is an atomic size parameter -- see radial in m515.
c   generate pts in the new variable x which are evenly spaced at xi=i/n.
c   given these equally spaced points in x, generate the weights
c   for an extended newton-cotes quadrature which covers the region.
c   for this purpose, the grid is divided into rings which contain
c   m+1 points, and the mth order newton-cotes formulas generated.
c   note that the endpoints of one ring may be shared with the
c   next ring.
c
c   finally, include the term dr/dx in the weights returned. 
c            dr/dx=alpha/(1-x)**2         if mu=1
c            dr/dx=2*alpha*x/(1-x)**3     if mu=2
c            dr/dx=3*alpha*x**2/(1-x)**4  if mu=3
c     
c
c***references
c
c***routines called
c
c***end prologue       cotes.f
      implicit none
c     --- input variables -----
      integer mxrings,npts,rule
      real*8 alpha
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer nr(mxrings)
c     --- output variables ---
      integer nrings,nrtot,nwtot
      real*8 rpts(0:npts),rwts(*),jacob(0:npts)
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,nleft
      integer ptwt,pts,mu
      integer inp,iout
      real*8 one,two,three,step
c
      parameter (one=1.0d+00,two=2.0d+00,three=3.0d+00)
      parameter (mu=2)
c
      common/io/inp,iout
c
c     --- generate the points and dr/dx
      if(mu.eq.1) then
         do 10 i=0,npts-1
            rpts(i)=alpha*float(i)/float(npts-i)
            jacob(i)=alpha*float(npts*npts)
     $              /float((npts-i)*(npts-i))
   10    continue
      else if (mu.eq.2) then
         do 20 i=0,npts-1
            rpts(i)=alpha*float(i*i)/float((npts-i)*(npts-i))
            jacob(i)=two*alpha*float(npts*npts*i)
     $              /float((npts-i)*(npts-i)*(npts-i))
   20    continue
      else if(mu.eq.3) then
         do 30 i=0,npts-1
            rpts(i)=alpha*float(i*i*i)
     $              /float((npts-i)*(npts-i)*(npts-i))
            jacob(i)=three*alpha*float(npts*npts*i*i)
     $              /float((npts-i)*(npts-i)*(npts-i)*(npts-i))
   30    continue
      endif
c
      nrtot=npts
c
c     --- determine the number of rings to use for this atom.
c         divide the range into a number of rings containing 'rule' points,
c         and use a 'rule-1'-order newton-cotes quadrature in each ring.
c         the last ring uses whatever is left over.
c         note that the endpoint of ring1 is the beginning point of ring2, etc.
c
c     --- the biasing by -1 in the loop over points is a little sleaze.
c         in order to do the newton-cotes quadrature over the intervals
c         it is helpful to have points duplicated -- i.e. [a,b], [b,c],
c         the common endpoints in each region would suggest that we
c         include radial points at the edges twice. however, this
c         means more function evaluations, potential problems in the
c         voronoi weights, etc.  it is simple to bias the pointer to
c         the radial points by one so as to avoid the overcounting
c         and do similar tricks in the actual quadrature to use
c         the common endpoints. however, the weights are a little
c         awkward. since each ring may have a different number of points
c         with a different order newton-cotes rule, the weights for
c         point b in ring 1 may be different from the weights used
c         for point b in ring 2.  for now, we get around this by 
c         duplicating the weights.
c
c         note that the weights array is really 2-dimensional.
c         e.g. for a 3-point quadrature, it is a (3,2) array
c         whose rows denote the points and columns denote the
c         integration region.  that is
c         Int(0,1) = wt(1,1)f(1) +wt(2,1)f(2) +wt(3,1)f(3)
c         Int(1,2) = wt(1,2)f(1) +wt(2,2)f(2) +wt(3,2)f(3)
c
c
c     this code is left over from days when we didn't know how many
c     rings there were -- can change it someday.
      step=one/float(npts)
      nrings=0
      ptwt=1
      pts=0
  100 continue
         if((pts+rule-1).lt.nrtot) then
            nrings=nrings+1
            nr(nrings)=rule
            call ncwts(rpts(pts),rwts(ptwt),rule,step,.false.)
            do 40 j=1,rule-1
               call vmul(rwts(ptwt+(j-1)*rule),rwts(ptwt+(j-1)*rule),
     $                   jacob(pts),rule)
   40       continue
            ptwt=ptwt+nr(nrings)*(nr(nrings)-1)
            pts=pts+nr(nrings)-1
         else
c           --- last ring
            nrings=nrings+1
            nleft=nrtot-pts
            nr(nrings)=nleft
            call ncwts(rpts(pts),rwts(ptwt),nleft,step,.false.)
            do 50 j=1,nleft-1
               call vmul(rwts(ptwt+(j-1)*nleft),rwts(ptwt+(j-1)*nleft),
     $                   jacob(pts),nleft)
   50       continue
            ptwt=ptwt+nr(nrings)*(nr(nrings)-1)
            pts=pts+nr(nrings)-1
         endif
      if(pts.lt.(nrtot-1)) goto 100
      nwtot=ptwt-1
c
c     --- test for an error
      if(nrings.ne.mxrings) then
         call lnkerr('problem with nrings in cotes')
      endif
c
c
      return
      end
