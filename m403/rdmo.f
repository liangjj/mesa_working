*deck @(#)rdmo.f	5.1  11/6/94
      subroutine rdmo(file,num,nummo,c,eigval,ops)
c***begin prologue     rdmo
c***date written       850601  yymmdd
c***revision date      890728  yymmdd
c
c  28 july 1989        bhl at llnl
c                      changing the format of rdinp so a blank
c                      read occurs before each vector. this is
c                      consistent with the scf=punch output in m501
c
c***keywords           molecular orbitals, i/o, rwf, chk
c***author             martin, richard (lanl)
c***source             @(#)rdmo.f	5.1   11/6/94
c***purpose            reads molecular orbitals/eigenvalues from a specific
c                      file.
c***description
c     call rdmo(file,num,nummo,c,eigval,ops)
c       file    the file to access.  may be 'rwf','chk',or 'inp'.
c       num     the basis set dimension.
c       nummo   the number of functions returned.
c       c       the eigenvector array(num,num).
c       eigval  the eigenvalue vector(num).
c       ops     the option string
c***references         (none)
c***routines called    iosys(io), lnkerr(mdutil), positn(chr)
c***end prologue       rdmo
      implicit integer(a-z)
      character*(*) file,ops
      character filnam*8,card*80
      real*8 c(num,num),eigval(num)
      logical positn
c
      common/io/inp,iout
c
c     nummo is identical to num except when the guess is read from the
c     input file, in which case you may read only a subset of the
c     full space (just the occupied orbitals, for instance).
c
      filnam=file
      if(filnam.eq.'rwf'.or.filnam.eq.'chk') then
         call iosys('read real "scf vector" from '//filnam,-1,c,0,' ')
         call iosys('read real "orbital energies" from '//filnam,
     $              -1,eigval,0,' ')
         nummo=num
      else if(filnam.eq.'inp') then
         if(.not.positn('$vectors',card,inp)) then
            call lnkerr(' no $vectors section found on '//filnam)
         endif
         nummo=intkey(ops,'guess=nocc',nummo,' ')
         do 20 j=1,nummo
            read(inp,*)
            read(inp,*) (c(i,j),i=1,num)
   20    continue
         write(iout,21)
   21    format(/,' input vectors from cards ')
         call matout(c,num,nummo,num,nummo,iout)
      else
         call lnkerr(' unrecognized file name in rdmo.')
      endif
c
c
      return
      end
