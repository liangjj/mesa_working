*deck @(#)hrdwrt.f	5.1  11/6/94
      subroutine hrdwrt(row)
c
c  1. read the hamiltonian matrix identifiers from iunt2a.
c     if ihamrd.ge.1 the identifiers are written out.
c
c  2. print the hamiltonian matrix if iciwrt.ge.2
c
      implicit real*8 (a-h,o-z)
      character*80 rtitle, ctitle, blabel
      integer and
c
      common /headng/ rtitle, ctitle, blabel
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
c
*mdc*if cray
*      dimension row(*), hij(4), i(4), j(4)
*      equivalence (next,rnext)
*mdc*else
      dimension row(*), hij(4), i(4), j(4), nexta(2)
      equivalence (nexta(1),rnext), (nexta(2),next)
*mdc*endif
c
      rewind iunt2a
      read (iunt2a) rtitle, blabel, ecore, thresh
      if(ihamrd.ge.1) write (iw,1000) rtitle, blabel, ecore, thresh
c
      if(iciwrt.lt.2) return
c
      write (iw,1400)
c
      next = 2
      icount = 0
c
      do 30 ii=1,nsef
c
      call hin(row,next)
c
      last = next - 1
      do 20 iii=1,last
      if(icount.eq.4) then
        write (iw,1500) (i(ij), j(ij), hij(ij), ij=1,4)
        icount = 0
      endif
      icount = icount + 1
      i(icount) = ii
      rnext = row(iii)
*mdc*if sun cray
      j(icount) =  and(next,mask)
*mdc*else
*      j(icount) = iand(next,mask)
*mdc*endif
   20 hij(icount) = row(iii)
c
   30 rnext = row(last+1)
c
      write (iw,1500) (i(ij), j(ij), hij(ij), ij=1,icount)
c
      return
c
 1000 format(///5x,'hamiltonian matrix label - ',a80/
     1       /5x,'integral tape label - ',a80/
     2       /10x,'core energy =',f20.8////10x,
     3   'threshold for the hamiltonian matrix elements =',1p,d10.1)
 1400 format(1h1/5x,41hthe hamiltonian matrix (lower triangular)///)
 1500 format( 4(4x,2i5,f15.8) )
c
      end
