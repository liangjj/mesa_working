*deck @(#)zget.f	5.2  2/5/95
      subroutine zget(nz,nvar,mxatm,toang,inau,inrad,ianz,z,iz,
     $                bl,alpha,beta,lbl,lalpha,lbeta,
     $                values,fpvec,intvec,anames,symbls,namcnt,bastyp,
     $                grdtyp,radius)
c***begin prologue     zget.f
c***date written       850601  yymmdd
c***revision date      2/5/95
c***keywords           z-matrix
c***author             binkley, et.al. guassian 82.
c***                   martin, richard (lanl)
c***source             @(#)zget.f	5.2   2/5/95
c***purpose            parses the z-matrix section.
c***description
c     call zget(nz,nvar,mxatm,toang,inau,inrad,ianz,z,iz,
c               bl,alpha,beta,lbl,lalpha,lbeta,
c               values,fpvec,intvec,anames,symbls,namcnt,bastyp,grdtyp,radius)
c***references
c***routines called    izero(math), rzero(math), zsymbl(m101), zvar(m101),
c                      iosys(io)
c
c***end prologue       zget.f
      implicit none
c     --- input variables -----
      integer nz,nvar,mxatm
      logical inau,inrad
      real*8 toang
c     --- input arrays (unmodified) ---
      real*8 z(mxatm)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ianz(mxatm),iz(4,mxatm),lbl(mxatm),lalpha(mxatm)
      integer lbeta(mxatm),intvec(mxatm)
      character*(*) anames(mxatm),namcnt(mxatm),bastyp(mxatm)
      character*(*) grdtyp(mxatm)
      character*80 symbls(mxatm)
      real*8 bl(mxatm),alpha(mxatm),beta(mxatm),values(mxatm)
      real*8 fpvec(mxatm),radius(mxatm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nsymbl,i,j
      logical ok
      real*8 one,tobohr,con,f180,torad,pi
c
      data one/1.0d+00/, f180/180.0d+00/
      save one,f180
c
      common/io/inp,iout
c
 1000 format(' integer parameters encountered on z-matrix card.',i3)
 1010 format(1x,'symbolic values:')
c
c     module to read the z-matrix and variables from the input stream.
c
c     initialize arrays.
      nsymbl=0
      call izero(lbl,mxatm)
      call izero(lalpha,mxatm)
      call izero(lbeta,mxatm)
      call izero(intvec,mxatm)
      call rzero(z,mxatm)
      call rzero(fpvec,mxatm)
      call rzero(bl,mxatm)
      call rzero(alpha,mxatm)
      call rzero(beta,mxatm)
      do 10 i=1,mxatm
         anames(i)=' '
         bastyp(i)=' '
         grdtyp(i)=' '
         namcnt(i)=' '
         symbls(i)=' '
   10 continue
c
c     get the symbolic z-matrix. upon return from zsymbl:
c        nsymbl ... total number of symbols (non-numbers) in z-matrix.
c        lbl   ... bond length flag.
c        lalpha... alpha angle flag.
c        lbeta ... beta angle flag.
c           the flags are:
c              1 = integer
c              2 = floating point
c              3 = symbol
c             -3 = -symbol
c
      call zsymbl(ianz,z,iz,bl,alpha,beta,lbl,lalpha,lbeta,
     $           nsymbl,nz,bastyp,grdtyp,radius,symbls,namcnt)
c
c     scan for integers as z-matrix parameters, and convert the
c     code used for constants.
      ok=.true.
      do 40 i=1,nz
         if(lbl(i)    .eq.2) lbl(i)    =0
         if(lalpha(i) .eq.2) lalpha(i) =0
         if(lbeta(i)  .eq.2) lbeta(i)  =0
         if(lbl(i).eq.1.and.lalpha(i).eq.1.and.lbeta(i).eq.1) then
            ok=.false.
            write(iout,1000) i
         endif
   40 continue
c
c
      nvar=0
      if(nsymbl.ne.0) then
c        only 0 and +/- 3 remains in lbl,lalpha, and lbeta.  convert
c        these to +/- 3000 so that they won't be confused with
c        pointers to variable 3.
         do 60 i=1,nz
            lbl(i)    = lbl(i)    * 1000
            lalpha(i) = lalpha(i) * 1000
            lbeta(i)  = lbeta(i)  * 1000
   60    continue
c
c        process variables section.
c        write(iout,1010)
         call zvar(bl,alpha,beta,lbl,lalpha,lbeta,values,
     $          fpvec,intvec,nsymbl,nz,nvar,anames,symbls)
      endif
c
c     convert to atomic units.
      tobohr=one/toang
      call iosys('read real pi from rwf',1,pi,0,' ')
      if(inau) tobohr=one
      torad=pi/f180
      if(inrad) torad=one
      do 80 i=1,nz
         bl(i)    =bl(i)    *tobohr
         alpha(i) =alpha(i) *torad
         beta(i)  =beta(i)  *torad
         radius(i)=radius(i)*tobohr
   80 continue
      if(nvar.ne.0) then
         do 90 i=1,nvar
            con=torad
            do 85 j=1,nz
               if(iabs(lbl(j)).eq.i) con=tobohr
   85       continue
         values(i)=values(i)*con
   90    continue
      endif
c
c
      return
      end
