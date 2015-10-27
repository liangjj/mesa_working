*deck @(#)dirsrt.f	5.1  11/6/94
      subroutine dirsrt
c
c***begin prologue     dirsrt
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)dirsrt.f	5.1   11/6/94
c***purpose            to move file lists to equalize the space left
c                      among various units.
c***description        #
c
c
c***references
c
c***routines called    lnkerr (mdutil)
c
c   common blocks:     (none)
c
c***end prologue       dirsrt
c
      implicit integer (a-z)
c
      character*80 tmplin
      character*32 keylis,key
      character*16 cconst,unlist,unityp,type,unit
      character*8  filtyp
      integer base,readpt,writpt,eof,end,iconst,unitpt,nfile
      integer un,file,nunits,unitno
      real*8 rconst
      logical align,locked
c
      parameter (maxfil=2000,maxun=20)
c
      common /ioqqq1/ keylis(maxfil),filtyp(maxfil),cconst(maxfil),
     #                unlist(maxun),type,key,unit,tmplin,unityp(maxun)
      common /ioqqq2/ base(maxfil),readpt(maxfil),writpt(maxfil),
     #                eof(maxfil),end(maxfil),iconst(maxfil),
     #                rconst(maxfil),unitpt(maxun),nfile(maxun),
     #                align(maxun),locked(maxun),unitno(maxun),
     #                un,file,nunits
c
c
c     ----- first pack the files down in the arrays -----
c
      pt=0
      do 2 unt=1,nunits
         ptold=unitpt(unt)
         do 1 i=1,nfile(unt)
            keylis(pt+i)=keylis(ptold+i)
            base(pt+i)=base(ptold+i)
            readpt(pt+i)=readpt(ptold+i)
            writpt(pt+i)=writpt(ptold+i)
            end(pt+i)=end(ptold+i)
            eof(pt+i)=eof(ptold+i)
            filtyp(pt+i)=filtyp(ptold+i)
            cconst(pt+i)=cconst(ptold+i)
            iconst(pt+i)=iconst(ptold+i)
            rconst(pt+i)=rconst(ptold+i)
c
            if (ptold.eq.pt) go to 1
            keylis(ptold+i)=' '
            base(ptold+i)=0
            readpt(ptold+i)=0
            writpt(ptold+i)=0
            end(ptold+i)=0
            eof(ptold+i)=0
            filtyp(ptold+i)=' '
            cconst(ptold+i)=' '
            iconst(ptold+i)=0
            rconst(ptold+i)=0.0d+00
    1    continue
         unitpt(unt)=pt
         pt=pt+nfile(unt)
    2 continue
c
c     ----- check that there is enough space for all the files -----
c
      if (pt.ge.maxfil) then
         call lnkerr('io system: exceeded the maximum number of files')
      end if
c
c     ----- now expand the list again, giving each unit the same
c           number of 'open' files. the funny business with pt
c           is to arrange that the last file gets the extra files
c           due to roundoff. thus, if adding a new unit and the
c           directory is almost full, the added unit will at least
c           have a few entries available.
c
      free=(maxfil-pt)/nunits
      pt=pt-nfile(nunits)+(nunits-1)*free
      do 4 unt=nunits,2,-1
         ptold=unitpt(unt)
         do 3 i=nfile(unt),1,-1
            keylis(pt+i)=keylis(ptold+i)
            base(pt+i)=base(ptold+i)
            readpt(pt+i)=readpt(ptold+i)
            writpt(pt+i)=writpt(ptold+i)
            end(pt+i)=end(ptold+i)
            eof(pt+i)=eof(ptold+i)
            filtyp(pt+i)=filtyp(ptold+i)
            cconst(pt+i)=cconst(ptold+i)
            iconst(pt+i)=iconst(ptold+i)
            rconst(pt+i)=rconst(ptold+i)
c
            if(ptold.eq.pt) go to 3
            keylis(ptold+i)=' '
            base(ptold+i)=0
            readpt(ptold+i)=0
            writpt(ptold+i)=0
            end(ptold+i)=0
            eof(ptold+i)=0
            filtyp(ptold+i)=' '
            cconst(ptold+i)=' '
            iconst(ptold+i)=0
            rconst(ptold+i)=0.0d+00
    3    continue
         unitpt(unt)=pt
         pt=pt-free-nfile(unt-1)
    4 continue
c
c
      return
      end
