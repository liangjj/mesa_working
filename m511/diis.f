*deck @(#)diis.f	5.1  11/6/94
      subroutine diis(f,t1,t2,t3,t4,t5,densty,s,smhalf,error,
     $                diissv,bmatrx,focksv,c,diiser,stdiis,
     $                nnp,nbf,mxds,mxiter,iter,ndiis,mndiis,
     $                ptr,dsptr,ops)
c
c***begin prologue     diis.f
c***date written       840527  (yymmdd)
c***revision date      11/6/94
c   12 june, 1993      rlm at lanl
c      modifying to call errmat to compute the error vector over
c      an orthogonalized ao basis. also changing to represent the error
c      vector in this basis as opposed to the initial mo vectors when
c      diis is turned on. finally, linear dependence defined as the
c      determinant less than 10**-12, as opposed to the previous 10**-15.
c   15 may, 1993       rlm at lanl
c      updated to delete dead code.
c***keywords           pulay,diis,extrapolate
c***author             saxe, paul (lanl)
c***source             @(#)diis.f	5.1   11/6/94
c***purpose            perform pulay's diis extrapolation
c***description
c     
c    
c
c***references
c   p.pulay, j.comp. chem. 3, 556(1982).
c   t.p.hamilton and p.pulay, jcp 84, 5728(1986).
c
c***routines called
c                      sdot,trtosq,flin,ebc,ebtc,ambc,saxpy
c                      matout,vmove,rzero,print,error
c***end prologue       diis.f
      implicit none
c     --- input variables -----
      integer stdiis,nnp,nbf,mxds,mxiter,iter,ndiis,mndiis
c     --- input arrays (unmodified) ---
      character*(*) ops
      real*8 c(nbf,nbf),densty(nnp),s(nnp),smhalf(nnp)
c     --- input arrays (scratch) ---
      real*8 t1(nbf,nbf),t2(nbf,nbf),t3(nbf,nbf),t4(nbf,nbf)
c     note that t1 and bmatrx are implicitly equivalenced.
      real*8 t5(nbf,nbf),bmatrx(*)
c     --- output arrays ---
      integer ptr(mxiter),dsptr(mxds)
      real*8 f(nnp),error(nbf,nbf,mxds),diissv(mxiter)
      real*8 focksv(nnp,mxds)
c     --- output variables ---
      real*8 diiser
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,n,iq
      logical logkey,debug
      real*8 scalar
      real*8 sdot
      real*8 zero,one,ten,tenm12
c
      common /io/     inp,iout
c
      parameter (zero=0.0d+00,one=1.0d+00,ten=10.0d+00)
      parameter (debug=.false.)
      parameter (tenm12=1.0d-12)
c
 1010 format (' the diis vector in the a. o. basis')
 1020 format (' the diis error vector in the original mo basis')
 1030 format (' small matrix for diis ')
 1040 format (9x,' diis iteration',i4,', determinant: ',e9.1)
 1050 format (' solution vector for diis ',/,
     $                (1x,2i4,f15.9))
 1060 format (' fock matrix after diis ')
c
c     --- compute the diis error vector fds-sdf ---
c     note that pulay recommends that this be transformed to an
c     orthogonalized ao  basis to generate a more "balanced" error vector.
c     the maximum element of this orthogonalized error matrix is
c     returned in diiser, and is used to test convergence.
c     the orthogonalized error matrix, s**-1/2*(fds-sdf)*s**-1/2
c     is returned in t5.
      call errmat(f,t1,t2,t3,t4,t5,densty,s,diiser,smhalf,nnp,nbf,ops)
      diissv(iter)=diiser
c
c     --- check if error small enough for diis procedure ---
      if (diiser.lt.ten**(-stdiis)) then
c
         if (logkey(ops,'scf=print=ao-diis-error',.false.,' ')) then
            write (iout,1010)
            call matout(t5,nbf,nbf,nbf,nbf,iout)
         endif
c
         ndiis=ndiis+1
c
c        --- find an available diis storage spot ---
         do 30 i=1,mxds
            if (dsptr(i).le.0) then
               dsptr(i)=ndiis
               ptr(ndiis)=i
               go to 40
            endif
  30    continue
c
c        --- none available, so use oldest one ---
         i=ptr(mndiis)
         ptr(mndiis)=-9999998
         dsptr(i)=ndiis
         ptr(ndiis)=i
         mndiis=mndiis+1
c
   40    continue
c
         call vmove(error(1,1,ptr(ndiis)),t5,nbf**2)
         call vmove(focksv(1,ptr(ndiis)),f,nnp)
         if (ndiis.ge.2) then
c
c           --- form small matrix ---
  100       continue
            bmatrx(1)=zero
            n=1
            do 110 i=mndiis,ndiis
               n=n+1
               bmatrx(n)=-one
  110       continue
            scalar=sdot(nbf**2,error(1,1,ptr(mndiis)),1,
     $                 error(1,1,ptr(mndiis)),1)
            do 130 i=mndiis,ndiis
               n=n+1
               bmatrx(n)=-one
               do 120 j=mndiis,ndiis
                  n=n+1
                  bmatrx(n)=sdot(nbf**2,error(1,1,ptr(i)),1,
     $                           error(1,1,ptr(j)),1)/scalar
  120          continue
  130       continue
            n=n+1
            bmatrx(n)=-one
            do 140 i=mndiis,ndiis
               n=n+1
               bmatrx(n)=zero
  140       continue
c
c           --- solve small set of equations ---
            if (logkey(ops,'scf=print=diis-matrix',.false.,' ')) then
               write (iout,1030)
               call matout(bmatrx,ndiis-mndiis+2,ndiis-mndiis+3,
     $                     ndiis-mndiis+2,ndiis-mndiis+3,iout)
            endif
c           note that hamilton and pulay scale the b-matrix so that
c           b(ndiis,ndiis)=1 in order to avoid underflow.  so far we haven't
c           had a problem with that.
            call flin(bmatrx,ndiis-mndiis+2,ndiis-mndiis+2,1,scalar)
c
            if(debug) then
               write (iout,1040) ndiis,scalar
            endif
c
c           --- check if problems solving the equations ---
c               if there is, throw out earlier iterates until they behave.
            if (abs(scalar).lt.tenm12) then
               i=ptr(mndiis)
               dsptr(i)=-1
               ptr(mndiis)=-9999998
               mndiis=mndiis+1
               if (mndiis.eq.ndiis-1) go to 200
               go to 100
            endif
c
c           --- form new fock matrix ---
            call rzero(f,nnp)
            n=(ndiis-mndiis+2)**2+1
            if(debug) then
               write (iout,1050) (iq,ptr(iq),bmatrx(n+ptr(iq)),
     $                           iq=mndiis,ndiis)
            endif
c
            do 150 i=mndiis,ndiis
               n=n+1
               call saxpy(nnp,bmatrx(n),focksv(1,ptr(i)),1,f,1)
  150       continue
            if(debug) then
               write (iout,1060)
               call print(f,nnp,nbf,iout)
            endif
         endif
      endif
c
  200 continue
c
c
      return
      end
