*deck @(#)setind.f	1.2  7/30/91
      subroutine setind(a,pkndx,nbf,newnbf,bflabl,nulabl,
     #                 nstart,bstart,blen,nocart,nocont,
     $                 nat,ntypes,prnt,lenbf)
      implicit integer(a-z)
      integer a(*),pkndx(*),zero(10)
      character*80 line,cjunk
      character*16 bflabl(*),nulabl(*)
      character*8 aname,blanc
      character*8 chrkey
      character*3 center,blank
      character*3 ctype(20)
      character*1 stype,ptype,dtype,ftype,type
      character*3 drop(10),keep(10)
      logical posinp,find
      logical prnt
      logical logkey
      integer lenbf(nat)
      integer nstart(nat),bstart(nat,ntypes),blen(nat,ntypes)
      integer nocart(ntypes),nocont(nat,ntypes),ptr(4)
      common /io/ inp,iout
      data blank/'   '/,blanc/'        '/
      data ctype/'s','x','y','z',
     #           'xx','yy','zz','xy','xz','yz',
     #           'xxx','yyy','zzz','xyy','xxy','xxz','xzz',
     #           'yzz','yyz','xyz'/
      data stype/'s'/,ptype/'p'/,dtype/'d'/,ftype/'f'/
      data ptr/0,1,4,10/
c
c
      nztot=0
c
      do 100 i=1,nbf
         pkndx(i)=i
 100  continue
c
      find=posinp('$drop',cjunk)
      if(.not.find) then
         write(iout,*)' cant find $drop in the input'
         call lnkerr(' m302: cant find $drop')
      end if
c
   1  continue
         read (inp,845) line
 845     format (a80)
         call locase(line,line)
         ok=0
         if(.not.logkey(line,'$end',.false.,' ')) then
            ok=1
c
            do 2 i=1,10
               drop(i)=blank
               keep(i)=blank
   2        continue
c
c -- read the center number or the center name
c
            atom=intkey(line,'atom',-1,' ')
            center=chrkey(line,'center',' ',' ')
c
c -- read the type of function to be dropped
c
            type=chrkey(line,'type',' ',' ')
c
c -- determine which components are to be retained
c
            write(iout,*) '     '//line(1:40)
            idrop=0
            if(logkey(line,'drop',.false.,' ')) then
               idrop=1
               call getdrp(line,drop,nzero)
            else if(logkey(line,'keep',.false.,' ')) then
               call getkep(line,keep,nkeep)
            else
               write(iout,*)' error: drop or keep must be specified '
               call lnkerr(' m302: input error ')
            end if
c
c -- determine the number of these functions which are to be processed
c
            number=intkey(line,'number',-1,' ')
c
c -- match the center name
c
            if(atom.eq.-1) then
               do 10 i=1,nat
               ic=nstart(i)
               if(logkey(bflabl(ic),center,.false.,' ')) then
                  atom=i
                  go to 11
               end if
  10           continue
               write(iout,*)' center = ',center
               call lnkerr(' m302: cant find this atomic center')
  11           continue
            end if
c
            angmom=-1
c
            if(type.eq.stype) angmom=1
            if(type.eq.ptype) angmom=2
            if(type.eq.dtype) angmom=3
            if(type.eq.ftype) angmom=4
c
            if(angmom.eq.-1) then
               write(iout,*)' cant match this orbital type ',type
               call lnkerr(' m302: input error ')
            end if
c
            ncart=nocart(angmom)
            ncont=nocont(atom,angmom)
c
            write(iout,*)'     ncart ncont ',ncart,ncont
c
            if(idrop.eq.1) then
c
               a0=ptr(angmom)
               do 40 j=1,nzero
                  id=-1
                  do 41 i=1,ncart
                     if(drop(j).eq.ctype(a0+i)) then
                        id=i
                        go to 42
                     end if
   41             continue
   42             continue
                  if(id.eq.-1) then
                     write(iout,*)' error in 40 loop atom angmom ',
     $                            atom,angmom
                     call lnkerr(' m302: error ')
                  end if
                  zero(j)=id
   40          continue
            else
               a0=ptr(angmom)
               do 49 i=1,ncart
                  zero(i)=i
  49           continue
               do 50 j=1,nkeep
                  id=-1
                  do  51 i=1,ncart
                     if(keep(j).eq.ctype(a0+i)) then
                        id=i
                        go to  52
                     end if
   51             continue
   52             continue
                  if(id.eq.-1) then
                     write(iout,*)' error in 50 loop atom angmom ',
     $                            atom,angmom
                     call lnkerr(' m302: error ')
                  end if
                  zero(id)=0
   50          continue
c23456
               i0=1
               do 53 i=1,ncart
                  if(zero(i).ne.0) then
                     zero(i0)=i
                     i0=i0+1
                  end if
c23456
   53          continue
               nzero=i0-1
               nzt=ncart-nkeep
               if(nzero.eq.0 .or. nzero.ne.nzt) then
                  write(iout,*)' error in keep nzero nkeep nzt ',
     #                           nzero,nkeep,nzt
                  call lnkerr(' m302 keep section')
               end if
 
            end if
c
            if(number.eq.-1) then
               ns=bstart(atom,angmom)
               ne=ns+blen(atom,angmom)-1
               nztot=nztot+nzero*ncont
            else
               ns=bstart(atom,angmom)+(number-1)*ncart
               ne=ns+ncart-1
               nztot=nztot+nzero
            end if
c
            do 3 i=ns,ne,ncart
               do 4 j=1,nzero
                  iz=zero(j)
                  pkndx(i-1+iz)=0
    4          continue
    3       continue
c
c
         end if
c
         if(ok.eq.1) go to 1
c
      newnbf=nbf-nztot
c
      i0=0
      do 99 i=1,nbf
         if(pkndx(i).ne.0) then
            i0=i0+1
            pkndx(i)=i0
         end if
  99  continue
c
      is=1
      i0=0
      do 96 i=2,nbf
         if(bflabl(i)(1:8).ne.blanc) then
            i0=i0+1
            lenbf(i0)=i-is
            is=i
         end if
  96  continue
      lenbf(i0+1)=nbf-is+1
c
      j0=0
      i0=0
      do 98 i=1,nat
         aname=bflabl(j0+1)(1:8)
         iat=0
         do 97 j=1,lenbf(i)
            j0=j0+1
            if(pkndx(j0).ne.0) then
               i0=i0+1
               nulabl(i0)=bflabl(j0)
               if(iat.eq.0) then
                  nulabl(i0)(1:8)=aname
                  iat=1
               end if
            end if
  97     continue
  98  continue
c23456
      if(prnt) then
         write(iout,*)' '
         write(iout,*)' orbital   index   label '
         do 101 i=1,nbf
            write(iout,102) i,pkndx(i),bflabl(i)
 101     continue
 102    format(2x,i5,2x,i5,2x,a16)
c
         write(iout,*)' '
         write(iout,*)' orbital   new label '
         do 103 i=1,newnbf
            write(iout,104) i,nulabl(i)
 103     continue
 104    format(2x,i5,2x,a16)
      end if
c
c
      return
      end
