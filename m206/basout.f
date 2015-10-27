*deck @(#)basout.f	1.1 9/7/91
c***begin prologue     basout
c***date written       901214   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6010, link 6010, kohn data
c***author             schneider, barry (lanl)
c***source             m6010
c***purpose            basis from mesa to compatible
c***                   file for kohn codes.
c
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       basout
      subroutine basout(coords,ex,cont,ptprim,noprim,nocont,ptcont,nat,
     1                 nprim,ntypes,nbtype,nnp,ncont,start,nbasis,zan,
     2                 nocart,nobf,maxmom,mintyp,nx,ny,nz,minmom,index,
     3                 idrop,oldnbf,ncon,prbaso,prbasn)
      implicit integer (a-z)
      real *8 ex, cont, zan, coords, eta
      character*3 ctype
      dimension ex(nprim), cont(ncont), zan(nat), ptprim(nat,ntypes)
      dimension noprim(nat,ntypes), nocont(nat,ntypes)
      dimension ptcont(nat,ntypes), start(nat,ntypes), nocart(ntypes)
      dimension nobf(ntypes), maxmom(ntypes), mintyp(ntypes)
      dimension minmom(ntypes), nx(*), ny(*), nz(*), index(oldnbf)
      dimension coords(3,nat), eta(300,5), nstart(300), nstop(300)
      dimension lnew(300), mnew(300), nnew(300), cnt(300)
      dimension ctype(0:3,0:3,0:3)
      common/io/inp,iout
      logical prbaso, prbasn
      ctype(0,0,0)='s'
      ctype(1,0,0)='x'
      ctype(0,1,0)='y'
      ctype(0,0,1)='z'
      ctype(2,0,0)='xx'
      ctype(0,2,0)='yy'
      ctype(0,0,2)='zz'
      ctype(1,1,0)='xy'
      ctype(1,0,1)='xz'
      ctype(0,1,1)='yz'
      ctype(3,0,0)='xxx'
      ctype(0,3,0)='yyy'
      ctype(0,0,3)='zzz'
      ctype(2,1,0)='xxy'
      ctype(2,0,1)='xxz'
      ctype(1,2,0)='xyy'
      ctype(0,2,1)='yyz'
      ctype(1,0,2)='xzz'
      ctype(0,1,2)='yzz'
      ctype(1,1,1)='xyz'
c----------------------------------------------------------------------c
c     retrieve information about the most demanding shell block        c
c----------------------------------------------------------------------c
      call iosys('read integer maxprm from rwf',1,maxprm,0,' ')
      call iosys('read integer maxcont from rwf',1,mxcont,0,' ')
      call iosys('read integer maxl from rwf',1,maxl,0,' ')
      call iosys('read integer maxblk from rwf',1,maxblk,0,' ')
      call iosys('read integer d1maxblk from rwf',1,dlen,0,' ')
      call iosys('read integer dolp from rwf',1,dolp,0,' ')
c----------------------------------------------------------------------c
c     read in basis set information from read-write file               c
c----------------------------------------------------------------------c
      call iosys('read real exponents from rwf',-1,ex,0,' ')
      call iosys('read real "contraction coefficients" from rwf',
     $     -1,cont,0,' ')
      call iosys('read real "nuclear charges" from rwf',-1,zan,0,' ')
      call iosys('read real coordinates from rwf',-1,coords,0,' ')
      call iosys('read integer "pointer to primitives" from rwf',
     $     -1,ptprim,0,' ')
      call iosys('read integer "number of primitives" from rwf',
     $     -1,noprim,0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $     ' from rwf',-1,ptcont,0,' ')
      call iosys('read integer "number of contraction coefficients" '//
     $     'from rwf',-1,nocont,0,' ')
      call iosys('read integer "number of cartesians" from rwf',
     $     -1,nocart,0,' ')
      call iosys('read integer "number of pure functions" from rwf',
     $     -1,nobf,0,' ')
      call iosys('read integer "minimum momentum" from rwf',
     $     -1,minmom,0,' ')
      call iosys('read integer "maximum momentum" from rwf',
     $     -1,maxmom,0,' ')
      call iosys('read integer "pointer to cartesians" from rwf',
     $     -1,mintyp,0,' ')
      call iosys('read integer "pointer to first function" from rwf',
     $     -1,start,0,' ')
      call iosys('read integer "power of x" from rwf',-1,nx,0,' ')
      call iosys('read integer "power of y" from rwf',-1,ny,0,' ')
      call iosys('read integer "power of z" from rwf',-1,nz,0,' ')
c----------------------------------------------------------------------c
c             write out the location and charge of the atoms           c
c----------------------------------------------------------------------c
c
      nstart(1) = 1
      ikount = 0
      iloc = 0
c
c.. bhl 9/13/89         using index array to drop functions
c
      if(idrop.ne.0) then
          call iosys('read integer "packing index vector" from rwf',
     $                oldnbf,index,0,' ')
      else
          do 999 i=1,oldnbf
             index(i)=i
 999      continue
      endif
c
      kbf=0
c
c.. bhl 9/13/89             kbf counter
c
      do 9 iatom=1,nat
         write(iout,101) iatom, zan(iatom)
         do 7 itype=1,nbtype
            if (noprim(iatom,itype).gt.0) then
                mini=mintyp(itype)
                maxi=mintyp(itype)+nocart(itype)-1
                imax=maxmom(itype)
                write(iout,102) itype,imax
                nprimi=noprim(iatom,itype)
                nconti=nocont(iatom,itype)
                write(iout,103) nprimi, nconti
c
c  build arrays for quad codes: nstart, nstop, lnew, mnew, nnew, eta
c
c
c  the way the  indexing for quad codes works is that there are
c  ncontra contracted functions.  for the ith contracted function nstart(i)
c  and nstop(i) give the first and last locations of the relevant quantities
c  in the arrays lnew(i), mnew(i), nnew(i), and eta(i,j=1,5)
c   lnew, mnew, and new are the powers of x, y, and z of the primitive
c   and eta(i,j=1,5) contains the center, exponent and contraction coefficient
c
                do 6 icont=1,nconti
                   do 3 m=mini,maxi
c..bhl 9/13/89                  kbf counter
                      kbf=kbf+1
                      if(index(kbf).ne.0) then
c..bhl 9/13/89                  if statement
                         ikount = ikount + 1
                         do 4 iii=1,nprimi
                            iiirel = iii-1+ (icont-1)*nprimi
                            if (cont(ptcont(iatom,itype)+iiirel).gt.
     1                                           1.e-20) then
                                iloc = iloc+1
                                cnt(iloc)=iatom
                                eta(iloc,1) = coords(1,iatom)
                                eta(iloc,2) = coords(2,iatom)
                                eta(iloc,3) = coords(3,iatom)
                                eta(iloc,4) = ex(ptprim(iatom,itype)
     1                                           +iii - 1)
                                eta(iloc,5) = cont(ptcont(iatom,itype)
     1                                             +iiirel)
                                lnew(iloc) = nx(m)
                                mnew(iloc) = ny(m)
                                nnew(iloc) = nz(m)
                            endif
    4                   continue
                        nstop(ikount) = iloc
                        nstart(ikount+1) = nstop(ikount)+ 1
                      endif
c..bhl 9/13/89                  end if
    3              continue
    6           continue
c
      if(prbaso) then
         write(iout,104)
         call matout(ex(ptprim(iatom,itype)),nprimi,1,
     $               nprimi,1,iout)
      endif
      if(prbaso) then
         write(iout,105)
         call matout(cont(ptcont(iatom,itype)),nprimi,nconti,
     $               nprimi,nconti,iout)
      endif
      if (prbaso) then
          write(iout,106)
          do 5 m=mini,maxi
             write(iout,107) nx(m), ny(m), nz(m)
    5      continue
      endif
           endif
    7    continue
    9 continue
      npr = iloc
      ncon = ikount
c
c write the newly indexed basis set information to kohndt.
c
c
      write(iout,*) '   number primitives ',i4,
     1              '   number contracted ',i4
      write(iout,*) '   beginning and ending primitive for each '//
     1                  'contracted function'
      do 20 i=1,ncon
         write(iout,*) nstart(i), nstop(i)
   20 continue
      if(prbasn) then
         write (iout,109)
         write (iout,110)
         do 12 i=1,ncon
            nlow = nstart(i)
            nhi = nstop(i)
            do 11 j=nlow,nhi
               write(iout,111) lnew(j),mnew(j),nnew(j),
     1                         eta(j,4),eta(j,5),eta(j,1),
     2                         eta(j,2),eta(j,3)
   11       continue
   12    continue
      endif
c
      return
  101 format (/,5x,'atom number',1x,i4,2x,'charge',f15.8)
  102 format (/,5x,'type primitive',1x,i2,2x,'maximum angular momentum',
     1        1x,i2)
  103 format (/,5x,'no. primitives',1x,i3,2x,'no. contracted',1x,i3)
  104 format(/,5x,'primitive exponents')
  105 format(/,5x,'contraction coefficients')
  106 format (/,5x,'xyz convention')
  107 format (/,5x,'total no. primitive functions',1x,i4,1x,'total no. c
     1ontracted functions',1x,i4)
  109 format (/,5x,'basis set in new format')
  110 format(/,6x,' l ',' m ',' n ','      exp       ',
     1     '  cont. coef.  ','    Rx    ','    Ry    ','    Rz    ',/)
  111 format(5x,3i3,2e15.8,3f10.5)
      end
