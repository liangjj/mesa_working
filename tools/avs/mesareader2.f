      subroutine avsinit_modules
      include 'avs/avs.inc'
      external mesa_r
      integer mesa_r,param
      integer dummy

      call AVSset_module_name('Mesa Checkfile Reader','data')
      dummy = AVScreate_output_port('molecule','molecule')
      param=AVSadd_parameter('Check file','string',' ',' ',' ')
      call AVSconnect_widget(param,'browser')
      param=AVSadd_parameter('Use only active sites','boolean',0,0,0)
      call AVSconnect_widget(param,'toggle')
      call AVSset_compute_proc(mesa_r)
      call AVSload_user_data_types('/usr/avs/MSI/include/mol_udat.h')
      call AVSload_user_data_types('/usr/avs/MSI/include/atom_udat.h')
      return
      end

      function mesa_r(mol_out,file_name,usenact)
      integer mesa_r
      include 'avs/avs.inc'
      include 'chemistry/CHEMlong.inc'
      integer mol_out
      integer error
      integer TRUE,FALSE
      character*(*) file_name
      character*3 answer
      integer nat,usenact,nactiv
      integer mesa_read
      logical opened
      save opened
      data opened/.false./

      parameter (MAXATM=1000)
      parameter (MAXBONDS=6)
      parameter (TRUE=1)
      parameter (FALSE=0)

      mesa_r=FALSE

      if ( file_name .ne. ' ') then
         opened = .true.
         call iosys('open chk as old',1000000,0,0,file_name)
         call iosys('read integer "number of atoms" from chk',1,nat,0,
     $        ' ')
         call iosys('does "number of atoms with basis functions" '//
     $        'exist on chk',0,0,0,answer)
         if (answer .eq. 'yes') then
            call iosys('read integer '//
     $           '"number of atoms with basis functions"'//
     $           ' from chk',1,nactiv,0,' ')
         else
            nactiv=nat
         endif 

         call mesa_read(mol_out,nat,usenact,nactiv,error)
         call iosys('close chk',0,0,0,' ')
         mesa_r=error
      else
         mol_out=0
         mesa_r=TRUE
      endif
      
      close(unit=66)
      return
      end

      subroutine mesa_read(mol_out,nat,usenact,nactiv,success)
      implicit none
      include 'avs/avs.inc'
      include 'chemistry/CHEMlong.inc'
      integer mol_out,success
      integer nat
      integer nactiv,usenact,na
      real*8 c(3,nat)
      integer ian(nat)
      integer mol, atom_list, atom,quantum
      integer nbf,nprim,TRUE,FALSE,i,npf,ncont,nbtype,ntypes
      integer do_color,setmsimol,setmsiatom
      integer ANGSTROMS,BOHR
      parameter (ANGSTROMS=0)
      parameter (BOHR=1)
      save mol
      data mol/0/
      character*2 pertable(103)
      character*256 title
      data pertable/'H ','He',
     $ 'Li','Be','B ','C ','N ','O ','F ','Ne',
     $ 'Na','Mg',
     $                'Al','Si','P ','S ','Cl','Ar',
     $ 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     $                'Ga','Ge','As','Se','Br','Kr',
     $ 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     $                'In','Sn','Sb','Te','I ','Xe',
     $ 'Cs','Ba','La',
     $     'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     $     'Ho','Er','Tm','Yb','Lu',
     $                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     $                'Tl','Pb','Bi','Po','At','Rn',
     $ 'Fr','Ra','Ac',
     $     'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf',
     $     'Es','Fm','Md','No','Lr'/
      parameter (TRUE=1)
      parameter (FALSE=0)

      success=FALSE

      call iosys('read real coordinates from chk',-1,c,0,' ')
      call iosys('read integer "atomic numbers" from chk',-1,ian,0,' ')

      if (mol .ne. 0) then
         call  CHEMmolecule_free(mol)
      endif
      mol=CHEMmolecule_alloc()
      atom_list=CHEMatom_alloc()
      atom=atom_list
      na=nat
      if (usenact .eq. 1) na=nactiv
c
      do 10 i=1,na
         if (chematom_set_inumber(atom,i).ne. 0) return 
         if (chematom_set_name(atom,pertable(ian(i))).ne.0)return
         if (chematom_set_color(atom,do_color(ian(i))).ne.0)return
         if (chematom_set_xyz(atom,c(1,i),c(2,i),c(3,i)).ne.0)return
         if (setmsiatom(atom,ian(i)).ne.0) return
         if (chematom_add(atom_list,atom).ne.0) return
         atom=chematom_alloc()
 10   continue 
c
c free the last one, coz it wasn't used.
c
      call chematom_free(atom)
c
C the atom list is now allocated.  Now we have to do make the 
c quantum stuff.
c
      call iosys('read integer "number of basis functions" from chk',
     $     1,nbf,0,' ')
      call iosys('read integer "number of primitive functions"'//
     $     'from chk',1,npf,0,' ')
      call iosys('length of exponents on chk',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on chk',
     $     ncont,0,0,' ')
      call iosys('read integer "number basis types" from chk',
     $     1,nbtype,0,' ')
      call iosys('length of "number of pure functions" on chk',
     $     ntypes,0,0,' ')
      call read_quant(quantum,nat,nbf,npf,nprim,ncont,nbtype,ntypes,
     $     c,success)
      if (success.ne.TRUE) return
c         
      call iosys('read character "title" from chk',-1,0,0,title)
      if (chemmolecule_set_atom(mol,atom_list).ne.0) return
      if (chemmolecule_set_name(mol,title).ne.0) return
      if (chemmolecule_set_units(mol,BOHR).ne.0) return
      if (chemmolecule_set_natom(mol,nat).ne.0) return
      if (chemquantum_add(quantum,quantum).ne.0) return
      if (chemmolecule_set_quant(mol,quantum).ne.0) return
      if (chemmolecule_set_nquant(mol,1).ne.0) return
      if (setmsimol(mol).ne.0) return
      mol_out=mol
c      call myctofdeal(mol)
      success=TRUE
      return
      end

      subroutine read_quant(quantum,nat,nbf,npf,nprim,ncont,nbtype,
     $     ntypes,c,success)
      implicit none
      include 'avs/avs.inc'
      include 'chemistry/CHEMlong.inc'
      integer nat,nbf,npf,nprim,ncont,nbtype,ntypes
      integer ptprim(nat,ntypes),noprim(nat,ntypes),start(nat,ntypes),
     $     ptcont(nat,ntypes),nocart(0:ntypes),nobf(ntypes),
     $     minmom(ntypes),maxmom(ntypes),nocont(nat,ntypes)
      real*8 cont(ncont),ex(nprim),c(3,nat)
      integer quantum,shell,shell_list,gauss,gauss_list,icdegree,
     $     prim,icptr,thecoef,icontr,mystart
      integer success,TRUE,FALSE
      parameter (TRUE=1)
      parameter (FALSE=0)
      integer ichg,imul,ialp,ibet,bf,ngauss,nshell,nbasis,iatom,itype
      integer qpage,qpage_list
      success=FALSE

      call iosys('read integer "pointer to primitives" from chk',
     $     -1,ptprim,0,' ')
      call iosys('read integer "number of primitives" from chk',
     $     -1,noprim,0,' ')
      call iosys('read integer "pointer to contraction coefficients"'//
     $     ' from chk',-1,ptcont,0,' ')
      call iosys('read integer "number of contraction coefficients"'//
     $     ' from chk',-1,nocont,0,' ')
      call iosys('read integer "pointer to first function" from chk',
     $     -1,start,0,' ')
      call iosys('read integer "number of cartesians" from chk',
     $     -1,nocart,0,' ')
      call iosys('read integer "number of pure functions" from chk',
     $     -1,nobf,0,' ')
      call iosys('read integer "minimum momentum" from chk',
     $     -1,minmom,0,' ')
      call iosys('read real exponents from chk',-1,ex,0,' ')
      call iosys('read real "contraction coefficients" from chk',
     $     -1,cont,0,' ')
      call iosys('read integer charge from chk',1,ichg,0,' ')
      call iosys('read integer "spin multiplicity" from chk',1,
     $     imul,0,' ')
      call iosys('read integer "number of alpha electrons" from chk',
     $     1,ialp,0,' ')
      call iosys('read integer "number of beta electrons" from chk',
     $     1,ibet,0,' ')
      quantum=CHEMquantum_alloc()
      shell_list=chemshell_alloc()
      shell=shell_list

      if (chemquantum_set_nbasis(quantum,nbf).ne.0) return
      if (chemquantum_set_ich(quantum,ichg).ne.0)return
      if (chemquantum_set_mul(quantum,imul).ne.0)return
      if (chemquantum_set_ne(quantum,ialp+ibet).ne.0) return
      if (chemquantum_set_na(quantum,ialp).ne.0) return
      if (chemquantum_set_nb(quantum,ibet).ne.0) return
cKLUDGE
      if (chemquantum_set_scftype(quantum,1).ne.0)return
      if (chemquantum_set_name(quantum,'FOOBAR').ne.0)return
      bf=1
      ngauss=0
      nshell=0
      nbasis=0
      mystart=1
      do 1 iatom=1,nat
         do 2 itype=1,nbtype
            if (noprim(iatom,itype).gt.0) then
C
c there are basis functions of this type on this atom               
c     ChemViewer only groks up to f functions, and I don't feel like handling
c     sp[d[f]] blocks, so blow this off
c
               if (itype.gt.4) then
                  call AVSerror('Sorry, >F funs not supported yet')
                  return
               endif
c
c for each contraction of this type we make a "shell" for chemviewer.
c MESA uses one contraction degree for all types, but we'll fix that.
c we also need to duplicate primitives so that general contractions
c do not appear as such.
c
               do 40 icontr=0,nocont(iatom,itype)-1
                  nshell=nshell+1
                  gauss_list=chemgauss_alloc()
                  gauss=gauss_list
                  icdegree=0
                  do 50 prim=0,noprim(iatom,itype)-1
c pointer to contraction coefficient 
                     icptr=ptcont(iatom,itype)+
     $                    icontr*noprim(iatom,itype)+prim
                     if (cont(icptr).ne.0.0) then
                        ngauss=ngauss+1
                        thecoef=chemcoef_alloc()
                        if (chemcoef_set_val(thecoef,cont(icptr)).ne.0)
     $                       return
                        if (chemcoef_add(thecoef,thecoef).ne.0) return
                        if (chemgauss_set(gauss,
     $                       ex(ptprim(iatom,itype)+prim),thecoef,
     $                       c(1,iatom),c(2,iatom),c(3,iatom)).ne.0)
     $                       return
                        if (chemgauss_add(gauss_list,gauss).ne.0)
     $                       return
                        gauss=chemgauss_alloc()
                        icdegree=icdegree+1
                     endif
 50               continue 
c
c get rid of the spare
c
                  call chemgauss_free(gauss)
                  if (chemshell_set_katom(shell,iatom).ne.0) return
                  if (chemshell_set_ktype(shell,itype).ne.0) return
                  if (chemshell_set_kng(shell,icdegree).ne.0) return
                  if (chemshell_set_kstart(shell,mystart).ne. 0) return
                  if (chemshell_set_kloc(shell,bf).ne.0) return
                  if (chemshell_set_kmin(shell,itype).ne.0) return
                  if (chemshell_set_kmax(shell,itype).ne.0) return
                  if (chemshell_set_gauss(shell,gauss_list).ne.0)
     $                 return
                  if (chemshell_add(shell_list,shell).ne.0) return
                  mystart=mystart+icdegree
c kludge -- assume only pure s,p,d,f functions, no spdf junque

                  bf=bf+nocart(minmom(itype))
                  nbasis=nbasis+nocart(minmom(itype))
                  shell=chemshell_alloc()
 40            continue 
            endif
 2       continue 
 1    continue 
      if (bf .ne. nbf+1) then
         call AVSerror('Internal error, nbf != bf')
         return
      endif
c     clean it up
      call chemshell_free(shell)
      if (chemquantum_set_shell(quantum,shell_list).ne.0) return
      if (chemquantum_set_nshell(quantum,nshell).ne.0) return
      if (chemquantum_set_ngauss(quantum,ngauss).ne.0) return

c
c now get mos and mo energies
c
      qpage_list=chemq_page_alloc()
      qpage=qpage_list
      if (qpage_list.eq.0) return
      call read_nrg(qpage,nbf,success)
      if (success.eq.FALSE) return
      if (chemq_page_add(qpage_list,qpage).ne.0) return
      success=FALSE
      call read_mos(qpage_list,nbf,nat,ntypes,start,nocont,success)
      if (success.eq.FALSE) return
      success=FALSE
      if (chemquantum_set_q_page(quantum,qpage_list).ne.0) return
      if (chemquantum_set_npage(quantum,2).ne.0) return
      success=TRUE
      return
      end

      function do_color(indx)
      integer do_color,indx
      include 'chemistry/CHEMlong.inc'
      integer BLUE,GREEN,RED,BLUE_GR,WHITE,BLACK,GREY,YELLOW,MAGNT,
     $     PURPLE,PUNK
      integer colors(103)
      save colors
      logical called
      data called/.false./
      save called

      if (.not. called) then
         err=0
         err=chemgen_util_rgb_to_int(GREEN,0.0,1.0,0.0)
         err=chemgen_util_rgb_to_int(RED,1.0,0.0,0.0)
         err=chemgen_util_rgb_to_int(PUNK,1.0,0.5,0.5)
         err=chemgen_util_rgb_to_int(BLUE,0.0,0.0,1.0)
         err=chemgen_util_rgb_to_int(BLUE_GR,0.0,1.0,1.0)
         err=chemgen_util_rgb_to_int(WHITE,1.0,1.0,1.0)
         err=chemgen_util_rgb_to_int(BLACK,0.0,0.0,0.0)
         err=chemgen_util_rgb_to_int(GREY,0.7,0.7,0.7)
         err=chemgen_util_rgb_to_int(YELLOW,1.0,1.0,0.0)
         err=chemgen_util_rgb_to_int(MAGNT,0.0,1.0,1.0)
         err=chemgen_util_rgb_to_int(PURPLE,1.0,0.0,1.0)
      

C for now, make all H&He white, 1st row grey, second pinkish
C third yellow, fourth magenta, fifth (except lanthanides) purple
c lanthanides red, last row blue-green, actinides red
c
c but special cases: Oxygen red, nitrogen blue, carbon green.


         colors(1)=WHITE
         colors(2)=WHITE
         do 1 i=3,10
            colors(i)=GREY
 1       continue 
         do 2 i=11,18
            colors(i)=PUNK
 2       continue 
         do 3 i=19,36
            colors(i)=YELLOW
 3       continue 
         do 4 i=37,54
            colors(i)=MAGNT
 4       continue 
         do 5 i=55,57
            colors(i)=PURPLE
 5       continue
         
         do 15 i=58,71
            colors(i)=RED
 15      continue 
         do 6 i=72,86
            colors(i)=PURPLE
 6       continue 
         do 7 i=87,89
            colors(i)=BLUE_GR
 7       continue 
         do 8 i=90,103
            colors(i)=RED
 8       continue 
         colors(6)=GREEN
         colors(7)=BLUE
         colors(8)=RED
      endif

      if ( indx.ge.1 .and. indx .le. 103) then
         do_color=colors(indx)
      else
         do_color=0
      endif

      return
      end


      subroutine read_nrg(qpage,nbf,success)
      implicit none
      include 'avs/avs.inc'
      include 'chemistry/CHEMlong.inc'
      integer success,TRUE,FALSE,qpage,nbf,i,ierror,fp
      integer mysetqpagemtx
      real*8 eig(0:nbf-1)
      parameter (TRUE=1)
      parameter (FALSE=0)
      
      success=FALSE
      call iosys('read real "orbital energies" from chk',
     $     nbf,eig,0,' ')
      if (chemq_page_init_page(qpage,nbf,1).ne.0) return
      fp=chemgen_util_stderr()
      if (chemq_page_set_name(qpage,'Alpha Orbital Energies').ne.0)
     $     return
      ierror=mysetqpagemtx(qpage,nbf,1,eig)


      success=TRUE

      return
      end

      subroutine read_mos(qpage_list,nbf,nat,ntypes,start,nocont,
     $     success)
      implicit none
      include 'avs/avs.inc'
      include 'chemistry/CHEMlong.inc'
      integer success,TRUE,FALSE,qpage,nbf,i,j,ierror,itype,ibf,k
      integer qpage_list
      integer mysetqpagemtx
      real*8 mos(0:nbf-1,0:nbf-1),denmat(0:nbf-1,0:nbf-1)
      real*8 d((nbf*(nbf+1))/2)
      integer nat,ntypes,start(nat,ntypes),nocont(nat,ntypes),nnp
      integer nae
      parameter (TRUE=1)
      parameter (FALSE=0)

      success=FALSE
      call iosys('read integer "number of alpha electrons" from chk',1,
     $     nae,0,' ')
      qpage=chemq_page_alloc()
      if (qpage.eq.0) return
      call iosys('read real "scf vector" from chk',
     $     nbf*nbf,mos,0,' ')
      do 10 i=1,nat
         if (start(i,3).ne.0) then
            do 40 j=1,nocont(i,3)
               do 50 k=0,2
                  do 30 ibf=0,nbf-1
                     
                     mos(start(i,3)+k+(j-1)*6,ibf)=
     $                    mos(start(i,3)+k+(j-1)*6,ibf)/sqrt(3.0d0)
 30               continue 
 50            continue 
 40         continue 
         endif
 10   continue 
      if (chemq_page_init_page(qpage,nbf,nbf).ne.0) return
      if (chemq_page_set_name(qpage,'Alpha MO Coefficients').ne.0)
     $     return
      ierror = mysetqpagemtx(qpage,nbf,nbf,mos)

      if (chemq_page_add(qpage_list,qpage).ne.0) return

      do 100 i=0,nbf-1
         do 110 j=0,nbf-1
            denmat(i,j)=0.0
            do 120 k=0,nae-1
               denmat(i,j)=denmat(i,j)+2*mos(i,k)*mos(j,k)
 120        continue 
 110     continue 
 100  continue 
      nnp=nbf*(nbf+1)/2
      call sqtotr(d,denmat,nbf,nnp)
      qpage=chemq_page_alloc()
      if (qpage .eq. 0) return
      if (chemq_page_init_page(qpage,nnp,1).ne.0) return
      if (chemq_page_set_name(qpage,'Total SCF Density') .ne.0)
     $     return
      ierror = mysetqpagemtx(qpage,nnp,1,d)
      
      if (chemq_page_add(qpage_list,qpage).ne.0) return

      success=TRUE

      return
      end


*deck @(#)sqtotr.f	4.1  7/7/93
      subroutine sqtotr(triang,square,num,nnp)
c***begin prologue     sqtotr
c***date written       850601  (yymmdd)
c***revision date      910729  (yymmdd)
c
c   29   july  1991    rlm at lanl
c      modifying so that if the sqare matrix is not symmetric
c      only the first bad element is printed and not the entire
c      matrix.
c***keywords           matrix, pack
c***author             saxe, paul (lanl)
c***source             @(#)sqtotr.f	4.1   7/7/93
c***purpose            packs a real symmetric matrix into a linear array.
c                      call sqtotr(triang,square,num,nnp)
c                        triang   output matrix, nnp words long.
c                        square   input matrix, dimensioned (num,num).
c                        nnp      num*(num+1)/2
c
c***references
c***routines called    abs
c***end prologue       sqtotr
      implicit integer (a-z)
c
      real*8 triang(nnp),square(num,num)
      logical notify
c
      common/io/inp,iout
c
      notify=.true.
      ij=0
      do 2 j=1,num
         do 1 i=1,j
            if (abs(square(i,j)-square(j,i)).gt.1.0d-05) then
               if (abs(square(i,j)-square(j,i))/
     #            max(abs(square(i,j)),abs(square(j,i)))
     #                 .gt.1.0d-05) then
               if(notify) then
                  write (iout,9) i,j,square(i,j),square(j,i)
               endif
               notify=.false.
    9          format (//,' ##### library: sqtotr, the square',
     #                 ' matrix is not symmetric:',2i4,2g18.9,//)
ctemp               call lnkerr('non-symmetric matrix')
               end if
            end if
c
            ij=ij+1
            triang(ij)=square(i,j)
    1    continue
    2 continue
c
      return
      end
