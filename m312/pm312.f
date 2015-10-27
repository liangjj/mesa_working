*deck @(#)pm312.f	5.2  2/5/95
      subroutine pm312
c***begin prologue     pm312.f
c***date written       840710   (yymmdd)
c***revision date      2/5/95
c
c   20 september 1994  rlm at lanl
c        when doing intergrals, the int file (produced by m330)
c        is now destroyed if it exists before the raw integrals
c        are written.
c   11 june 1991  bhl at llnl, rlm at lanl
c        putting raw unsorted ao integrals on their own file (rints).
c        this filename is specified in m1/inirwf as "raw integral filename"    
c
c   27 june 1989  bhl at llnl, rlm at lanl.
c        added an option to skip the integral calculation
c        without using a nonstandard route .. the code
c        scans for the option  noints and if found calls
c        chainx.  in addition, one can use int=reuse, or
c        int=reuse1, int=reuse2 to recompute either 1e- or
c        2e- integrals only.
c
c   16 february 1987  pws at lanl
c        commenting out test on 'flbls', since factor of two on
c        32 bit machines causes problems. subsequent shuffle should
c        stop any problems
c
c   15 february 1987 pws at lanl
c        changing so that integrals unit is named with character
c        variable namint
c
c   10 june 1985   pws at lanl
c        modified fairly drastically to reduce memory
c        requirements as well as to test integrals for
c        negligible pre-exponential factors (pws,lanl).
c
c***keywords           m312, link 312, integrals, cutoff, preexponential
c***author             saxe, paul    (lanl)
c***source             @(#)pm312.f	5.2   2/5/95
c***purpose            computes symmetry-orbital 2-e integrals over a
c                      general contraction scheme.
c***description
c
c***references
c
c***routines called
c
c***end prologue       pm312.f
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a
      real*8 z
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer lenbuf
      integer n,intkey,nat,nbasis,nprim,ncont,ntypes,nbtype,lenxyz
      integer ndrop,newnbf,len,nintgr,zan,ptprim,sints,cont,c
      integer noprim,nocont,ptcont,nocart,nobf,minmom,maxmom
      integer mintyp,start,nx,ny,nz,pkindx,labels,nnp,mnnp
      integer ex,need
      integer core,maxcor,iadtwp,wpadti,wptoin, ngot, idum
      character*4096 ops
      character*8 prtflg
      character*128 rints,ints
      real*8 cutoff,cutexp
      logical logkey,prnt
      logical drop
      pointer (p,z(1)), (p,a(1))
c
      common /io/     inp,iout
      common /toler/  cutoff
c
      data lenbuf /524288/
      data prnt/.true./
      save lenbuf,prnt
c
c     --- lenbuf is number of integrals to write in a buffer ---
      lenbuf=lenbuf/512*512
c
c     --- recover the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- see if we are to compute 2e- integrals.
      if(logkey(ops,'noints',.false.,' ').or.
     $   logkey(ops,'int=reuse',.false.,' ').or.
     $   logkey(ops,'int=reuse2',.false.,' ').or.
     $   logkey(ops,'int=reuse=m312',.false.,' ')) then
         write(iout,*)'m312: skip two-electron integrals '
      else
c        --- has printing been turned off externally?
         call iosys('read character "print flag" from rwf',
     $               -1,0,0,prtflg)
         if(prtflg.eq.'minimum') prnt=.false.
c
c        --- set cutoffs for preexponential and final integral ---
c        --- cutoff is the threshhold for writing out integrals ---
         n=intkey(ops,'int=preex',15,' ')
         cutexp=float(n)
         cutexp=cutexp*log(10.0d+00)
         n=intkey(ops,'int=cutoff',10,' ')
         cutoff=10.0d+00**(-n)
c
         if(prnt) then
            write (iout,1) exp(-cutexp),cutoff
    1       format(1x,'m312:'
     $             /,5x,'pre-exponential integral cutoff ',1pe8.1,
     $             /,5x,'integral cutoff                 ',1pe8.1)
         endif
c
c        --- get lengths, dimensions, etc. needed for core allocation --
         call iosys('read integer "number of atoms" from rwf',1,nat,
     $               0,' ')
         call iosys('read integer "number of basis functions" from rwf',
     $              1,nbasis,0,' ')
         call iosys('length of exponents on rwf',nprim,0,0,' ')
         call iosys('length of "contraction coefficients" on rwf',
     $               ncont,0,0,' ')
         call iosys('length of "number of pure functions" on rwf',
     $               ntypes,0,0,' ')
         call iosys('read integer "number basis types" from rwf',
     $               1,nbtype,0,' ')
c
c        --- ntypes is the total number of types used to define the lengths
c            in m301.  this total is composed of two classes:  the first
c            nbtype refer to the basis function types, the remainder refer to
c            ecp types, and are not used here.
         call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
         nnp=(nbasis+1)*nbasis/2
c        
c        --- before opening the raw integral file, destroy the sorted integral 
c            file if it exists.
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,ints)
         call iorm(ints)
c
c        --- check for option to drop integral components ---
         drop=logkey(ops,'drop',.false.,' ')
         ndrop=intkey(ops,'drop',0,' ')
         newnbf=nbasis-ndrop
         mnnp=newnbf*(newnbf+1)/2
c
c        --- open the integral unit ---
c            len is a length of the files used...the fancy footwork
c            is a guess at how long small files will be.
c            n**4 integrals and labels plus few extra words
         len=3*(mnnp**2)/4 +5*nnp +20000
         len=max(len,100000)
         len=min(len,5000000)
         call iosys('read character "raw integral filename" from rwf',
     $               0,0,0,rints)
c
         if(logkey(ops,'unit=ssd=rints',.false.,' ')) then
            call iosys('open rints as new on ssd',len,0,0,rints)
         else
            call iosys('open rints as new',len,0,0,rints)
         end if
c
         call iosys('write integer lenbuf to rints',1,lenbuf,0,' ')
         call iosys('create real file "unsorted ao integrals" on rints',
     $               -1,0,0,' ')

c        --- the labels require 64 bits for packing. the following
c            variable is '1' on 64 bit machines, '2' on 32 bit
c            machines. it gives the number of integer words for
c            the integral labels
         nintgr=wptoin(1)
c
c        --- divide core for basis set information ---
         zan=1
         c=zan+nat
         cont=c+3*nat
         ex=cont+ncont
         sints=ex+nprim
         ptprim=wpadti(sints+lenbuf)
         noprim=ptprim+nat*ntypes
         nocont=noprim+nat*ntypes
         ptcont=nocont+nat*ntypes
         start=ptcont+nat*ntypes
         nocart=start+nat*ntypes
         nobf=nocart+ntypes
         minmom=nobf+ntypes
         maxmom=minmom+ntypes
         mintyp=maxmom+ntypes
         nx=mintyp+ntypes
         ny=nx+lenxyz
         nz=ny+lenxyz
         labels=nz+lenxyz
         pkindx=labels+nintgr*lenbuf
         core=iadtwp(pkindx+nbasis)
         need=wpadti(core)
c
c         call getscm(0,z(1),maxcor,'m312:main',0)
         call getmem(need,p,ngot,'m312:main',0)
c         if (need+20000.gt.maxcor) then
c            write(iout,*) 'real*8 available:',iadtwp(maxcor)
c            write(iout,*) 'real*8 needed:',iadtwp(need+20000)
c            call lnkerr('m312 needs more memory')
c         end if
c
c        --- read in basis set information from the read-write file ---
         call iosys('read real exponents from rwf',-1,z(ex),0,' ')
         call iosys('read real "contraction coefficients" from rwf',
     $              -1,z(cont),0,' ')
         call iosys('read real "nuclear charges" from rwf',
     $              -1,z(zan),0,' ')
         call iosys('read real coordinates from rwf',-1,z(c),0,' ')
         call iosys('read integer "pointer to primitives" from rwf',
     $              -1,a(ptprim),0,' ')
         call iosys('read integer "number of primitives" from rwf',
     $              -1,a(noprim),0,' ')
         call iosys('read integer "pointer to contraction coefficients"'
     $              //' from rwf',-1,a(ptcont),0,' ')
         call iosys('read integer "number of contraction coefficients"'
     $              //' from rwf',-1,a(nocont),0,' ')
         call iosys('read integer "number of cartesians" from rwf',
     $              -1,a(nocart),0,' ')
         call iosys('read integer "number of pure functions" from rwf',
     $              -1,a(nobf),0,' ')
         call iosys('read integer "minimum momentum" from rwf',
     $              -1,a(minmom),0,' ')
         call iosys('read integer "maximum momentum" from rwf',
     $              -1,a(maxmom),0,' ')
         call iosys('read integer "pointer to cartesians"from rwf',
     $              -1,a(mintyp),0,' ')
         call iosys('read integer "pointer to first function" from rwf',
     $              -1,a(start),0,' ')
         call iosys('read integer "power of x" from rwf',
     $              -1,a(nx),0,' ')
         call iosys('read integer "power of y" from rwf',
     $              -1,a(ny),0,' ')
         call iosys('read integer "power of z" from rwf',
     $              -1,a(nz),0,' ')
         if (drop) then
            call iosys('read integer "packing index vector" from rwf',
     $              -1,a(pkindx),0,' ')
         end if
c
c        --- calculate two electron integrals ---
         call driver(a(ptprim),a(noprim),nbtype,a(nocont),a(ptcont),
     $               z(ex),z(cont),ncont,z(c),z(sints),a(labels),lenbuf,
     $               a(start),a(nx),a(ny),a(nz),lenxyz,a(nobf),
     $               a(nocart),a(mintyp),a(minmom),a(maxmom),
     $               nat,nbasis,nnp,nprim,
     $               ops,cutexp,prnt,nintgr,drop,a(pkindx))
         call iosys('close rints',0,0,0,' ')
         call getmem(-ngot,p,idum,'m312:main',idum)
      endif
c
c     --- and exit gracefully ---
      call chainx(0)
c
c
      stop
      end
