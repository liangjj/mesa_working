*deck @(#)fmmult.f	5.1  11/6/94
      subroutine fmmult(xyzprp,xyz,npint,imax,jmax,
     $                  lmult,powx,powy,powz,
     $                  jcx,jcy,jcz)
c***begin prologue     fmmult.f
c***date written       860826  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl)
c***source             @(#)fmmult.f	5.1   11/6/94
c***purpose            assembles multipole integrlas 
c***description
c   module to assemble the one-dimensional integrals in xyz into
c   three-dimensional primitive multipole integrals in xyzprp.
c***references
c
c***routines called
c
c***end prologue       fmmult.f
      implicit none
c     --- input variables -----
      integer npint,imax,jmax,lmult
      integer powx,powy,powz
      real*8 jcx,jcy,jcz
c     --- input arrays (unmodified) ---
      real*8 xyz(npint,0:imax,0:jmax+lmult,3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 xyzprp(npint,0:imax,0:jmax,3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer len,j,k
      integer inp,iout
      real*8 bico,trans,one
      parameter (one=1.0d+00)
c
      common/io/inp,iout
c
c     --- the multipole integrals in xyz have been computed with 
c     respect to the coordinates of center j.  we are actually 
c     interested in the multipoles with respect to a third center c;
c     i.e., (x-Cx)**p = (x-Jx +Jx-Cx)**p.
c     this routine translates the origin of the multipole operator
c     to xc=(x-Cx). the translation component is jcx.
c
      len=npint*(imax+1)
      call rzero(xyzprp(1,0,0,1),len*(jmax+1)*3)
c
      do 50 j=0,jmax
c        --- handle the k=0 case individually. 
c            it's done separately to avoid the possibility of (0.0)**0
         trans=one
         call vwxs(xyzprp(1,0,j,1),xyzprp(1,0,j,1),
     $                    xyz(1,0,j+powx,1),trans,+1,len)
         call vwxs(xyzprp(1,0,j,2),xyzprp(1,0,j,2),
     $                    xyz(1,0,j+powy,2),trans,+1,len)
         call vwxs(xyzprp(1,0,j,3),xyzprp(1,0,j,3),
     $                    xyz(1,0,j+powz,3),trans,+1,len)
         do 40 k=1,powx
            trans=bico(powx,k)*(jcx**k)
            call vwxs(xyzprp(1,0,j,1),xyzprp(1,0,j,1),
     $                       xyz(1,0,j+powx-k,1),trans,+1,len)
   40    continue
         do 30 k=1,powy
            trans=bico(powy,k)*(jcy**k)
            call vwxs(xyzprp(1,0,j,2),xyzprp(1,0,j,2),
     $                       xyz(1,0,j+powy-k,2),trans,+1,len)
   30    continue
         do 20 k=1,powz
            trans=bico(powz,k)*(jcz**k)
            call vwxs(xyzprp(1,0,j,3),xyzprp(1,0,j,3),
     $                       xyz(1,0,j+powz-k,3),trans,+1,len)
   20    continue
   50 continue
c
c
      return
      end
