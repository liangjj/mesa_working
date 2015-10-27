*deck @(#)sympt.f	5.1  11/6/94
      subroutine sympt(nsym,nnpsym,nnqsym,isym,jsym,ijsym,ijpt,ijklpt,
     #                 symoff,nso,n1int,n2int,nnp,nbf)
c
      implicit integer (a-z)
c
      integer isym(nnpsym),jsym(nnpsym),ijsym(nnpsym),ijklpt(nnqsym)
      integer nso(nsym),ijpt(nnp),symoff(nsym)
c
      logical debug
c
      common /io/ inp ,iout
c
      parameter (debug=.false.)
c
      ioff(i,j)=i*(i-1)/2+j
c
c     ----- store the symmetry types for pairs -----
c
      symoff(1)=0
      do 30 i=2,nsym
         symoff(i)=symoff(i-1)+nso(i-1)
   30 continue
c
      do 40 i=1,nnp
         ijpt(i)=0
   40 continue
c
      n1int=0
      do 2 i=1,nsym
         do 1 j=1,i
            ij=ioff(i,j)
            isym(ij)=i
            jsym(ij)=j
            ijsym(ij)=xor(i-1,j-1)+1
    1    continue
         do 20 ii=1,nso(i)
            ia=ioff(symoff(i)+ii,symoff(i))
            do 19 jj=1,ii
               n1int=n1int+1
               ijpt(ia+jj)=n1int
   19       continue
   20    continue
    2 continue
c
c     ----- now pointers to the first of a symmetry chunk -----
c
      n2int=0
      do 4 ij=1,nnpsym
         do 3 kl=1,ij
            ijkl=ioff(ij,kl)
            if (xor(ijsym(ij),ijsym(kl)).eq.0) then
               ijklpt(ijkl)=n2int
               if (isym(ij).eq.jsym(ij).and.isym(kl).eq.jsym(kl)) then
                  n2int=n2int+ioff(nso(isym(ij)),nso(isym(ij)))*
     #                        ioff(nso(isym(kl)),nso(isym(kl)))
               else if (isym(ij).eq.jsym(ij)) then
                  call lnkerr('1')
                  n2int=n2int+ioff(nso(isym(ij)),nso(jsym(ij)))*
     #                        nso(isym(kl))*nso(jsym(kl))
               else if (isym(kl).eq.jsym(kl)) then
                  call lnkerr('2')
                  n2int=n2int+nso(isym(ij))*nso(jsym(ij))*
     #                        ioff(isym(kl),jsym(kl))
               else
                  n2int=n2int+nso(isym(ij))*nso(jsym(ij))*
     #                        nso(isym(kl))*nso(jsym(kl))
               end if
            else
               ijklpt(ijkl)=-999999999
            end if
    3    continue
    4 continue
c
c     ----- print some pertinent arrays -----
c
      if (debug) then
c
         write (iout,9)
    9    format (/,'  symmetry pointer arrays:')
         do 14 i=1,nsym
            do 13 j=1,i
               do 12 k=1,i
                  if (i.eq.k) then
                     lmax=j
                  else
                     lmax=k
                  end if
                  do 11 l=1,lmax
                     ij=ioff(i,j)
                     kl=ioff(k,l)
                     ijkl=ioff(ij,kl)
                     if (ij.eq.kl.and.isym(ij).eq.jsym(kl)) then
                        n=ioff(nso(isym(ij)),nso(isym(ij)))**2
                     else
                        n=nso(isym(ij))*nso(jsym(ij))*
     #                    nso(isym(kl))*nso(jsym(kl))
                     end if
                     write (iout,10) ijkl,i,j,k,l,ijsym(ij),ijsym(kl),
     #                               ijklpt(ijkl),n
   10                format (1x,i4,2x,4i1,3x,2i3,1x,i10,i5)
   11              continue
   12           continue
   13        continue
   14    continue
c
         write (iout,15) symoff
   15    format (/,' symoff: ',8i4)
         write (iout,16)
   16    format (/,' ijpt:')
         do 90 i=1,nbf
            i1=i*(i-1)/2+1
            i2=i*(i-1)/2+i
            write (iout,80) (ijpt(j),j=i1,i2)
   80       format (1x,20i4,(/,5x,20i4))
   90    continue
         write (iout,17) n1int,n2int
   17    format (/,' number of one-electron integrals in symmetry: ',i8,
     #           /,' number of two-electron integrals in symmetry: ',i8)
c
c
      end if
c
c
      return
      end
