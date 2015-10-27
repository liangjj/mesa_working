*deck drvhmt.f 
c***begin prologue     drvhmt
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             drvhmt
c***purpose            driver for hamiltonian matrix elements
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       drvhmt
      subroutine drvhmt(pt,kmat,coord,diag,n)
c
      implicit integer (a-z)
      character*(*) coord
      character*80 title
      real*8 pt, kmat, ham, dscale
      character*24 typ
      logical diag
      dimension pt(n), kmat(n,n)
      common/io/inp, iout
      data prn/.true./
      pointer (pham,ham(1))
c
      v=1
      need=v+n
      if(diag) then
         mat=need
         eig=mat+n*n
         work=eig+n
         need=work+5*n
      endif
      need=wpadti(need)
      call getmem(need,pham,ngot,'ham',0)
      call vone(coord,pt,ham(v),dscale,typ,n,.true.,.false.)
      if(diag) then
         call copy(kmat,ham(mat),n*n)
         call vec2di(ham(mat),ham(mat),ham(v),'add',n)
         call dsyev('v','l',n,ham(mat),n,ham(eig),ham(work),5*n,info)
         title='eigenvalues'
         call prntrm(title,ham(eig),n,1,n,1,iout)
         call getmem(-ngot,pham,idum,'ham',idum)
      endif
      return
      end


















