*deck @(#)errmat.f	5.1  11/28/95
      subroutine errmat(f,t1,t2,t3,t4,error,densty,s,diiser,
     $                  smhalf,nnp,nbf,ops)
c
c***begin prologue     errmat.f
c***date written       930612  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard (lanl)
c***source             @(#)errmat.f	5.1   11/28/95
c***purpose            computes the error matrix (fds-sdf)
c***description
c     
c    
c
c***references         p.pulay, j.comp.chem. 3, 556(1982).
c
c***routines called
c
c***end prologue       errmat.f
      implicit none
c     --- input variables -----
      integer nbf,nnp
c     --- input arrays (unmodified) ---
      character*(*) ops
      real*8 f(nnp),s(nnp),densty(nnp),smhalf(nnp)
c     --- input arrays (scratch) ---
      real*8 t1(nbf,nbf),t2(nbf,nbf),t3(nbf,nbf),t4(nbf,nbf)
c     --- output arrays ---
      real*8 error(nbf,nbf)
c     --- output variables ---
      real*8 diiser
c     --- scratch arrays ---
c     --- local variables ---
      logical logkey
      integer i,j,inp,iout
      real*8 zero
c
      common /io/     inp,iout
c
      parameter (zero=0.0d+00)
c
 1000 format (' the diis vector in the orthogonalized a. o. basis')
c
c     --- compute the diis error vector fds-sdf ---
      call trtosq(t1,f,nbf,nnp)
      call trtosq(t2,densty,nbf,nnp)
      call trtosq(t3,s,nbf,nnp)
      call ebc(t4,t1,t2,nbf,nbf,nbf)
      call ebc(error,t4,t3,nbf,nbf,nbf)
      call ebc(t4,t3,t2,nbf,nbf,nbf)
      call ambc(error,t4,t1,nbf,nbf,nbf)
c
c     --- transform the error matrix to an orthogonal basis using
c         s**(-1/2).
      call trtosq(t1,smhalf,nbf,nnp)
      call ebc(t2,error,t1,nbf,nbf,nbf)
      call ebtc(error,t1,t2,nbf,nbf,nbf)
c
      if (logkey(ops,'scf=print=diis-error',.false.,' ')) then
            write (iout,1000)
            call matout(error,nbf,nbf,nbf,nbf,iout)
      endif
c
c     --- find maximum element of error vector ---
      diiser=zero
      do 20 i=1,nbf
         do 10 j=1,i
            diiser=max(diiser,abs(error(i,j)))
   10    continue
   20 continue
c
c
      return
      end
