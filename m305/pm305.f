*deck @(#)pm302.f	5.1  11/6/94
      subroutine pm302(z,a)
c***begin prologue     pm302.f
c***date written       840710   (yymmdd)
c***revision date      11/6/94
c   january 8, 1994    rlm at lanl
c      adding capability to redo nuclear attraction integrals in
c      the event of a solvent calculation.
c   june   10, 1985    rlm at lanl
c      modified fairly drastically to compute
c      effective core potential integrals.
c***keywords           m302, link 302, one-electron, integrals, ecp
c***author             saxe, paul and martin, richard    (lanl)
c***source             @(#)pm302.f	5.1   11/6/94
c***purpose            computes symmetry-orbital 1-e integrals over a
c                      general contraction scheme.
c***description
c***references
c
c***routines called
c***end prologue       m302
      implicit none
c     --- input variables -----
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      integer a(*)
      real*8 z(*)
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer maxnbf
      parameter (maxnbf=2000)
c
      integer inp,iout
      integer maxcor,bufsiz,ngot
      integer nat,nbasis,nprim,ncont,ntypes,nbtype,lenxyz
      integer ncharge,zan,ptprim,s,nnp
      integer maxprm,mxcont,maxl,maxblk
      integer ex,cont,c,noprim,ptcont,nocont,nocart,nobf,maxmom,minmom
      integer mintyp,start,nx,ny,nz
      integer wpadti,iadtwp
      integer i,npint,prmint
      integer top,top1,top2,need,rneed
      character*4096 ops
      character*128 namint
      character*16 bflabl(maxnbf)
      character*3 answer
      logical dolp,logkey,debug
      real*8 half
c
      common /io/     inp,iout
c
      parameter (debug=.false.)
      parameter (half=0.5d+00)
c
      data bufsiz/4096/
      save bufsiz
c
 1000 format(1x,'m302: skip one-electron integrals')
c
c     --- collect the options string ---
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     --- see if we are to compute integrals.
      if(logkey(ops,'int=reuse',.false.,' ').or.
     $   logkey(ops,'int=reuse1',.false.,' ').or.
     $   logkey(ops,'noints',.false.,' ')) then
         write(iout,1000)
c        --- restore the rwf file.
c            get some buffer space.
         call getscm(bufsiz,z(1),ngot,'m302: buffer',0)
         call iosys('read character "integral filename" from rwf',
     $               0,0,0,namint)
         call iosys('open ints as old',0,0,0,namint)
         call iosys('copy "nuclear repulsion energy" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "overlap integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "kinetic integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
         call iosys('copy "potential integrals" from ints to rwf',
     $               bufsiz,a,0,'noerror')
      else
c
c        --- get the lengths of arrays needed for core allocation ---
c        ntypes is the total number of types used to define the lengths
c        in m102.  this total is composed of two sets: the first nbtype
c        types are used to mark the basis function types.  the remainder
c        refer to the ecp types.
c
         call iosys('read integer "number of atoms" from rwf',1,nat,
     $               0,' ')
         call iosys('read integer "number of basis functions" from rwf',
     $                  1,nbasis,0,' ')
         call iosys('length of exponents on rwf',nprim,0,0,' ')
         call iosys('length of "contraction coefficients" on rwf',
     $               ncont,0,0,' ')
         call iosys('length of "number of pure functions" on rwf',
     $               ntypes,0,0,' ')
         call iosys('read integer "number basis types" from rwf',
     $               1,nbtype,0,' ')
         call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
         nnp=(nbasis+1)*nbasis/2
c
c        --- see if solvent charges exist ---
c            if so, they will be added to the list of centers
c            for nuclear attraction integrals
         call iosys('does "number of solvent surface points"'
     $              //' exist on rwf',0,0,0,answer)
         if(answer.eq.'yes') then
            call iosys('read integer "number of solvent surface'
     $                 //' points" from rwf',1,ncharge,0,' ')
         else
            ncharge=0
         endif
c
c        --- divide core for basis set information ---
         zan=1
         c=zan+nat+ncharge
         cont=c+3*(nat+ncharge)
         ex=cont+ncont
         s=ex+nprim
         ptprim=wpadti(s+nnp)
         noprim=ptprim+ntypes*nat
         nocont=noprim+ntypes*nat
         ptcont=nocont+ntypes*nat
         start=ptcont+ntypes*nat
         nocart=start+ntypes*nat
         nobf=nocart+ntypes
         minmom=nobf+ntypes
         maxmom=minmom+ntypes
         mintyp=maxmom+ntypes
         nx=mintyp+ntypes
         ny=nx+lenxyz
         nz=ny+lenxyz
         need=nz+lenxyz
         rneed=iadtwp(need)
c
c        --- allocate core for integral computation.
c            retrieve information about the most demanding shell block.
         call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
         call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
         call iosys('read integer maxl from rwf',1,maxl,0,' ')
         call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
         call iosys('read integer dolp from rwf',1,dolp,0,' ')
c
c        --- potential-energy integrals are most-demanding in a normal run.
c            they need:
         npint=maxprm*maxprm
         prmint=1+16*npint+2*npint*(maxl+1)
         top1=prmint+npint*(maxblk+6*(maxl+1)*(maxl+1)+1)
         top1=max(top1,prmint+nbasis*nbasis)
         top2=prmint+maxblk*(npint+mxcont*mxcont) +mxcont*maxprm
c        --- effective core potentials?
         if(dolp) then
            top1=max(top1,npint*(11*9*9+maxblk)+310)
         endif
         top=need+wpadti(max(top1,top2))
c
         call getscm(top,z(1),maxcor,'m302: main',0)
         maxcor=iadtwp(maxcor)
c
c        --- read in basis set information from read-write file ---
         call iosys('read real exponents from rwf',-1,z(ex),0,' ')
         call iosys('read real "contraction coefficients" from rwf',
     $               -1,z(cont),0,' ')
         call iosys('read real "nuclear charges" from rwf',
     $               -1,z(zan),0,' ')
         call iosys('read real coordinates from rwf',-1,z(c),0,' ')
         call iosys('read integer "pointer to primitives" from rwf',
     $              -1,a(ptprim),0,' ')
         call iosys('read integer "number of primitives" from rwf',
     $              -1,a(noprim),0,' ')
         call iosys('read integer "pointer to contraction '//
     $              'coefficients" from rwf',-1,a(ptcont),0,' ')
         call iosys('read integer "number of contraction coefficients"'
     $              //' from rwf',-1,a(nocont),0,' ')
         call iosys('read integer "number of cartesians" from rwf',
     $              -1,a(nocart),0,' ')
         call iosys('read integer "number of pure functions" from rwf',
     $              -1,a(nobf),0,' ')
         call iosys('read integer "minimum momentum" from rwf',
     $               -1,a(minmom),0,' ')
         call iosys('read integer "maximum momentum" from rwf',
     $              -1,a(maxmom),0,' ')
         call iosys('read integer "pointer to cartesians" from rwf',
     $               -1,a(mintyp),0,' ')
         call iosys('read integer "pointer to first function" from rwf',
     $               -1,a(start),0,' ')
         call iosys('read integer "power of x" from rwf',-1,a(nx),0,' ')
         call iosys('read integer "power of y" from rwf',-1,a(ny),0,' ')
         call iosys('read integer "power of z" from rwf',-1,a(nz),0,' ')
         call iosys('read character "basis function labels" from rwf',
     $               -1,0,0,bflabl)
c        --- get solvent charges and coordinates ---
         if(ncharge.ne.0) then
            call iosys('read real "solvent surface charges" from rwf',
     $                  ncharge,z(zan+nat),0,' ')
            call iosys('read real "solvent surface coordinates"'
     $                  //' from rwf',3*ncharge,z(c+3*nat),0,' ')
c
c           --- note this little sleight of hand.
c               the point charges which come from m611 interact with
c               the solute in such a way that the potential is 
c               q(i)*q(j)/2*rij if i denotes solute and j denotes solvent
c               coordinates. the routines which compute the one-electron
c               integrals and the nuclear repulsion energy assume the
c               standard form q(i)q(j)/rij. in order to trick the code
c               into doing the different potential for all interactions
c               involving the solvent induced point charges, multiply
c               the solvent charges by 1/2 and proceed.
            call smul(z(zan+nat),z(zan+nat),half,ncharge)
            if(debug) then
               write(iout,*) 'ncharge',ncharge
               write(iout,*) 'charges',(z(zan+nat+i),i=0,ncharge-1)
               write(iout,*) 'coordinates',
     $                       (z(c+3*nat+i),i=0,3*ncharge-1)
            endif
         endif
c
c        --- calculate the nuclear repulsion energy ---
         call nucrep(z(zan),z(c),nat,ncharge)
c
c        --- calculate one-electron integrals ---
         call oneint(z(c),z(ex),z(rneed),a(need),z(cont),z(s),a(ptprim),
     $               a(noprim),a(nocont),a(ptcont),nat,nprim,
     $               maxcor-rneed+1,ntypes,nbtype,nnp,ncont,a(start),
     $               nbasis,z(zan),a(nocart),a(nobf),a(maxmom),
     $               a(mintyp),a(nx),a(ny),a(nz),a(minmom),dolp,
     $               bflabl,ops,ncharge)
      endif
c
c     --- and exit gracefully ---
      call chainx(0)
c
c
      stop
      end
