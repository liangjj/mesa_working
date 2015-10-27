*deck @(#)soanal.f	5.1  11/6/94
      subroutine soanal(row,ido,c,d,hb,hr,cthc,scr,root,temp,temp2,
     >bflabl)
c
c On option, (0) stores the vector coefficicents in the chk file,
c (1) performs the operation C^t H C from the coefficients in the chk 
c file and writes this to output file, diagonalizes this matrix,
c D^t[C^t H C]D, and writes out the vectors D and eigenvalues
c
c
      implicit real*8 (a-h,o-z)
      character*80 rtitle, ctitle, blabel
      character*128 namchk
      character*16 bflabl(nsef)
      integer and
c
      common /headng/ rtitle, ctitle, blabel
      common/icntrl/ iciwrt, icipun, intwrt, ihamrd, ksym
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
c
*mdc*if cray
*      dimension row(*), hij(4), i(4), j(4)
*      equivalence (next,rnext)
*mdc*else
      dimension c(nroots*nsef),d(nroots,nroots),hb(nsef,nsef),
     >hr(nroots,nroots),cthc(nroots,nroots),scr(*),root(nroots),
     >temp(nroots,nsef),temp2((nroots*(nroots+1))/2)
c
      dimension row(*),hij(4), i(4), j(4), nexta(2)
      equivalence (nexta(1),rnext), (nexta(2),next)
*mdc*endif
c
c  open chk file for read and write
c
      call iosys('read character "checkpoint filename" from rwf',
     >           0,0,0,namchk)
      call iosys('open chk as unknown',0,0,0,namchk)
      call iosys('read character "basis function labels" from rwf',
     >            -1,0,0,bflabl)
c
      if(ido.eq.0)then
c
c write vector to rwf file
c
c take transpose of c before writing out, so it has the form (nsef,nroots)
c
        write(iw,*)'writing coefficients to chk'
        do 5 ii=1,nroots
          do 5 jj=1,ii
            kk=(jj-1)*nroots+ii
            tempc=c(kk)
            ind=(ii-1)*nroots+jj
            c(kk)=c(ind)
            c(ind)=tempc
5       continue
        call iosys('write real "so vectors" to chk',
     >            nsef*nroots,c,0,' ')
        call wmat(c,nsef,nroots,bflabl,bflabl)
      end if
c
      if(ido.gt.0)then
c
c perform the operation C^t H C
c
c read vector from chk 
c
        call iosys('read real "so vectors" from chk',
     >            -1,c,0,' ')
        write(iw,*)'reading so vectors from chk'
c
c write vector to output
c
        call wmat(c,nsef,nroots,bflabl,bflabl)
c
c read hamiltonian
c
      do 10 ii=1,nsef
        do 10 jj=1,nsef
          hb(ii,jj)=0.d0
10    continue
c
      rewind iunt2a
      read (iunt2a) rtitle, blabel, ecore, thresh
      write(iw,1000)rtitle,blabel,ecore,thresh
c
      next = 2
      icount = 0
c
      do 30 ii=1,nsef
c
      call hin(row,next)
c
      last = next - 1
      do 20 iii=1,last
      if(icount.eq.4) then
c
c fill lower part of hamiltonian matrix
c
        do 25 ij=1,4
          hb(i(ij),j(ij))=hij(ij)
25      continue

        write(iw,*)'a',iii,ij,icount,last+1

        write(iw,1500)(i(ij),j(ij),hij(ij),ij=1,4)
c
        icount = 0
      endif
      icount = icount + 1
      i(icount) = ii
      rnext = row(iii)
*mdc*if sun cray
      j(icount) =  and(next,mask)
*mdc*else
*      j(icount) = iand(next,mask)
*mdc*endif
   20 hij(icount) = row(iii)
c
   30 rnext = row(last+1)

      write(iw,*)'b',iii,ij,icount,last+1

      write(iw,1500)(i(ij),j(ij),hij(ij),ij=1,icount)
c
c get last row
c
      do 33 ij=1,icount
        hb(i(ij),j(ij))=hij(ij)
33    continue
c
c fill in whole hamiltonian matrix
c
      do 35 ii=1,nsef
        do 35 jj=1,ii
          hb(jj,ii)=hb(ii,jj)
35    continue
      write(iw,*)'hamiltonian matrix, soanal'
      call wmat(hb,nsef,nsef,bflabl,bflabl)
c
c C^t H C
c
        call ebtc(temp,c,hb,nroots,nsef,nsef)
        call ebc(cthc,temp,c,nroots,nsef,nroots)
c
c if ido.gt.1 adjust diagonal elements
c
        if (ido.gt.1)then
          do 37 ii=1,ido
            read(5,*)adjust
            cthc(ii,ii)=adjust
37        continue
        end if
        
c write out cthc matrix to output
c
        write(iw,*)'cthc matrix'
        call wmat(cthc,nroots,nroots,bflabl,bflabl)
c
c write cthc matrix to chk 
c
        call iosys('write real "cthc matrix" to chk',
     >              nroots*nroots,cthc,0,' ')
c
c diagonalize cthc matrix
c
        nnp=(nroots*(nroots+1))/2
        call sqtotr(temp2,cthc,nroots,nnp)
        call givens(nroots,nroots,nroots,temp2,scr,root,d)
c
c write out eigenvalues
c
        write(iw,*)'eigenvalues of cthc '
        do 40 ii=1,nroots
          write(iw,*)ii,root(ii)
40      continue
c
c write eigenvectors to output
c
        write(iw,*)'eigenvectors of cthc'
        call wmat(d,nroots,nroots,bflabl,bflabl)
c
c write to output in free format
c
        do 45 ii=1,nroots
          write(iw,*)ii,(d(ii,jj),jj=1,3)
45      continue
          
c
c write eigenvectors to chk 
c
        call iosys('write real "dso vectors" to chk',
     >              nroots*nroots,d,0,' ')
      end if
c
c close chechkpoint file
c
      call iosys('close chk',namchk,0,0,' ')

1000  format(///5x,'hamiltonian matrix label - ',a80/
     1       /5x,'integral tape label - ',a80/
     2       /10x,'core energy = ',f20.8////10x,
     3    'theshold for the hamiltonian matrix elements=',1p,d10.1)
1500  format(4(4x,2i5,f15.8))
      return
      end
