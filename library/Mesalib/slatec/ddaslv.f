*deck ddaslv
      subroutine ddaslv (neq, delta, wm, iwm)
c***begin prologue  ddaslv
c***subsidiary
c***purpose  linear system solver for ddassl.
c***library   slatec (dassl)
c***type      double precision (sdaslv-s, ddaslv-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------------
c     this routine manages the solution of the linear
c     system arising in the newton iteration.
c     matrices and real temporary storage and
c     real information are stored in the array wm.
c     integer matrix information is stored in
c     the array iwm.
c     for a dense matrix, the linpack routine
c     dgesl is called.
c     for a banded matrix,the linpack routine
c     dgbsl is called.
c-----------------------------------------------------------------------
c***routines called  dgbsl, dgesl
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c***end prologue  ddaslv
c
      integer  neq, iwm(*)
      double precision  delta(*), wm(*)
c
      external  dgbsl, dgesl
c
      integer  lipvt, lml, lmu, lmtype, meband, mtype, npd
      parameter (npd=1)
      parameter (lml=1)
      parameter (lmu=2)
      parameter (lmtype=4)
      parameter (lipvt=21)
c
c***first executable statement  ddaslv
      mtype=iwm(lmtype)
      go to(100,100,300,400,400),mtype
c
c     dense matrix
100   call dgesl(wm(npd),neq,neq,iwm(lipvt),delta,0)
      return
c
c     dummy section for mtype=3
300   continue
      return
c
c     banded matrix
400   meband=2*iwm(lml)+iwm(lmu)+1
      call dgbsl(wm(npd),meband,neq,iwm(lml),
     *  iwm(lmu),iwm(lipvt),delta,0)
      return
c------end of subroutine ddaslv------
      end
