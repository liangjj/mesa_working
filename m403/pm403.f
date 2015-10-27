*deck @(#)pm403.f	5.1  11/6/94
      subroutine pm403(z,a)
c***begin prologue     m552
c***date written       910110   (yymmdd)
c
c
c***keywords           m403, link 403, overlap
c                     
c***author             braunstein,matt ; martin, richard 
c***source             @(#)pm403.f	5.1   11/6/94
c***purpose            checks orbital overlaps between two geometries
c***description
c     after determining orbitals check overlap of these orbitals
c     with those of a previous run and reorder them
c     m403 currently recognizes:
c       overlap=(order=(xx,xx))
c       where xx is an integer
c
c          like merge documentation
c       timing        print timing statistics for this link.
c
c***references
c
c***routines called
c     m403
c       drum(mdutil)
c       iosys(io)
c       traktm(mdutil)
c       getscm(mdutil)
c       sinv(util)
c       prjges(local)
c         sizmb(local)
c         rdbas(local)
c         rdmo(local)
c         minbas(local)
c           replic(util)
c           sto3g(local)
c             coefs(local)
c             scalef(local)
c             defsc(local)
c           togenc(util)
c           fmstrt(util)
c           ehuckl(local)
c             eneg(local)
c           normal(util)
c         ovrlap(local)
c           stwod(util)
c           fmonel(util)
c           trans1(util)
c           putrec(local)
c         hukmat(local)
c           getmo(local)
c         alter(local)
c       rwfges(local)
c         rdmo(local)
c         alter(local)
c       corges(local)
c         getmo(local)
c         alter(local)
c       wvec(util)
c       chainx(mdutil)
c
c***end prologue       m403
      implicit integer(a-z)
      integer a(*)
      real*8 z(*)
      real*8 maxerr
      character ops*4096, gestyp*8
      logical altges,prnt
      logical logkey
c
      common/io/inp,iout
c
      data maxcor/1/, maxerr/1.0d-06/, altges/.false./, prnt/.true./
      save maxcor,maxerr,altges,prnt
c
 1010 format(5x,'orbital overlaps')
c
c
c
c
c
c     retrieve the options string.
      call iosys('read character options from rwf',-1,0,0,ops)
      
      if(logkey(ops,'overlap=order',.false.,' ')) gestyp='order'
       write(iout,*)gestyp
c
c     start timing routines if enabled.
c
c     signal our presence.
      if(gestyp.eq.'order') write(iout,1010)
c
c     retrieve the basis set information.
      call iosys('read integer "number of basis functions" from rwf',
     $     1,nbasis,0,' ')
      nnp=nbasis*(nbasis+1)/2
      nbsq=nbasis*nbasis
c
c     allocate core needed by all options.
      c=1
      eigval=c+nbsq
      s=eigval+nbasis
      smhalf=s+nnp
      u=smhalf+nnp
      iu=wpadti(u)
      top=iu
c
      call getscm(top+30000,z,maxi,'m403',0)
      maxcor=iadtwp(maxi)
c
c     find orbital overlaps
      if(gestyp.eq.'order') then
         navl=maxcor-u+1
         call prjovr(gestyp,nbasis,nnp,z(c),z(eigval),z(smhalf),a(iu),
     $               z(u),navl,altges,ops)
      else
         call lnkerr('unrecognized guess type in m403:'//gestyp)
      endif
c
c
      return
      end
