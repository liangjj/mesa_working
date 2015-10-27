*deck @(#)natorb.f	5.1  11/6/94
      subroutine natorb(d1,nnp,natocc,eigvec,t1,t2,norbs,iout,bftorb,
     #                  aotono,occ,aotomo,motono,nbf,ops,bflabl,mcscf,
     #                  s,smhalf,tocc,u,eigval,nnpbf,triang)
c
c***begin prologue natorb
c***date written   851003   (yymmdd)
c***revision date  891128   (yymmdd)
c
c   28 november 1989   rlm at lanl
c      removing hubbard options.
c
c   02 february 1988   bhl at brl
c      return if mcscf failed to converge
c
c   11 january 1988    bhl at brl
c      changed to print nos
c
c   28 august 1986     modified by pws at lanl to not rotate mcscf
c      orbitals since that has been done in mcscf
c
c***keywords
c
c***author  saxe, paul,    (lanl)
c***purpose  to compute the natural orbitals for a ci wavefunction.
c
c***description  natorb diagonalizes the mo one-particle density
c         matrix, then forms a transformation matrix which includes
c         frozen core and frozen virtual orbitals. the original scf
c         vector is then transformed using this matrix to give the
c         ao to no transformation matrix. in this matrix, the no's
c         are ordered by descending occupation number.
c
c         on input:
c
c            d1       real(nnp)
c                     the mo one-particle density matrix.
c
c            nnp      integer
c                     the size of a triangular matrix 'norbs' on a side
c                     (norbs+1)*norbs/2
c
c            norbs    integer
c                     the number of orbitals in the ci.
c
c            iout     integer
c                     the name or number of the fortran output unit.
c
c            bftorb   integer (nbf)
c                     the correspondance between scf orbitals and the
c                     ci orbitals. -1 indicates a frozen core orbital;
c                     0, a frozen virtual.
c
c            aotomo   real (nbf,nbf)
c                     the scf vector.
c
c            nbf      integer
c                     the number of orbitals in the scf.
c
c            ops      character*(*)
c                     option string which may contain print_nos or
c                     print_occ
c
c            bflabl   character*(*)
c                     character string describing basis function labels.
c
c
c         on output:
c
c            natocc   real (norbs)
c                     the occupations of the natural orbitals in the
c                     ci basis, ordered in ascending order.
c
c            eigvec   real (norbs,norbs)
c                     the eigenvectors of the mo one-particle density
c                     matrix.
c
c            aotono   real (nbf,nbf)
c                     the ao to no transformation matrix.
c
c            occ      real (nbf)
c                     the occupancies of the no's in the scf basis,
c                     in descending order, corresponding to 'aotono'.
c
c            motono   real (nbf,nbf)
c                     the transformation matrix to take the scf
c                     orbitals into natural orbitals.
c
c
c         scratch:
c
c            t1       real (norbs)
c            t2       real (norbs)
c
c***references
c
c***routines called  rsp (clams), lnkerr (mdutil), rzero (util),
c                    ebc (math), wvec (util), itoc (char)
c***end prologue natorb
c
      implicit integer (a-z)
c
      real*8 d1(nnp),natocc(norbs),eigvec(norbs,norbs),aotono(nbf,nbf)
      real*8 occ(nbf),aotomo(nbf,nbf),motono(nbf,nbf)
      real*8 t1(norbs),t2(norbs)
      real*8 s(nnpbf),smhalf(nnpbf),tocc(nbf),eigval(nbf)
      real*8 u(nbf,nbf),triang(nnpbf)
      integer bftorb(nbf)
      character*(*) ops, bflabl(*)
      real*8 sum
      character*4 itoc
      character*32 cjunk,xform
      logical mcscf,logkey
c
c
    2 format(5x,'number of electrons:',f18.12)
   40 format(5x,'natural orbitals:')
   45 format(5x,'natural orbital occupation:',
     $          (/5x,5f12.6))
c
c     ----- diagonalize the guga-mo one-particle density matrix -----
c
      call rsp(norbs,norbs,nnp,d1,natocc,1,eigvec,t1,t2,ierr)
c
      if (ierr.ne.0) call lnkerr('error in rsp diagonalizing the '//
     #        'guga-mo one-particle density matrix:'//itoc(ierr))
c
      sum=0.0
      do 1 i=1,norbs
         sum=sum+natocc(i)
    1 continue
c
      write (iout,2) sum
c
c     ----- prepare a transformation matrix from scf mo's to no's -----
c
      call rzero(motono,nbf**2)
c
      nfzc=0
      do 10 i=1,nbf
         if (bftorb(i).eq.-1) then
            nfzc=nfzc+1
            motono(i,nfzc)=1.0d+00
            occ(nfzc)=2.0d+00
         end if
   10 continue
c
      do 20 bf=1,nbf
         mo=bftorb(bf)
         if (mo.gt.0) then
            do 15 i=norbs,1,-1
               motono(bf,norbs+1-i+nfzc)=eigvec(mo,i)
   15       continue
         end if
   20 continue
c
      do 25 i=norbs,1,-1
         occ(nfzc+norbs-i+1)=natocc(i)
         tocc(nfzc+norbs-i+1)=natocc(i)
   25 continue
c
      n=nfzc+norbs
      do 30 i=1,nbf
         if (bftorb(i).eq.0) then
            n=n+1
            motono(i,n)=1.0d+00
            occ(n)=0.0d+00
         end if
   30 continue
c
      if (mcscf) then
         cjunk='"mcscf failed"'
         call iosys('read character "transformation vector" from rwf',
     $           -1,0,0,xform)
         if(xform.eq.cjunk) return
         call iosys('read real "mcscf vector" from rwf',nbf**2,
     #               aotono,0,' ')
         call iosys('read real "mcscf orbital energies" from rwf',nbf,
     #               occ,0,' ')
      else
         call ebc(aotono,aotomo,motono,nbf,nbf,nbf)
      end if
c
c     ----- rotate degenerate orbitals -----
c
      call iosys('read real "overlap integrals" from rwf',nnpbf,s,0,' ')
      call sinv(s,smhalf,u,eigval,t1,t2,nbf,nnpbf,triang,0)
c..bhl
      if(logkey(ops,'ci=no=enforce',.false.,' ')) then
         call enforc(nbf,nnpbf,s,smhalf,tocc,aotono,t1,t2)
      endif
c
c..bhl
      if(logkey(ops,'print=ci=no',.false.,' ')) then
         write(iout,40)
         call wvec(aotono,occ,nbf,nbf,bflabl,' ')
      else if(logkey(ops,'print=ci=no_occ',.false.,' ')) then
         write(iout,45) occ
      endif
c..bhl
c
c
      return
      end
