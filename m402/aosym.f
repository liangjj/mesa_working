*deck @(#)aosym.f	5.1  11/6/94
      subroutine aosym(salc,numso,bfns,nirrep,nbf)
c
c***purpose: to determine the atomic orbitals associated with each
c            irreducible representation for kohn calculation
c
c barry schneider lanl
c
      implicit integer (a-z)
      integer numso(nirrep), bfns(nbf)
      real*8 salc(nbf,nbf), tol
      character *1 itoc
c
c
      common /io/ inp, iout
      parameter (tol=1.d-08)
c
c----------------------------------------------------------------------c
c              determine ao's in each irrep                            c
c----------------------------------------------------------------------c
      write (iout,1)
      call iosys ('write integer "kohn symmetries" to rwf',1,
     1             nirrep,0,' ')
      call iosys ('write integer "kohn ao numsym" to rwf',nirrep,numso,
     1             0,' ')
      offset=0
      do 100 irrep=1,nirrep
         call izero(bfns,nbf)
         if (numso(irrep).ne.0) then
             do 200 i=1,numso(irrep)
                do 300 j=1,nbf
                   if (abs(salc(j,i+offset)).ge.tol) then
                       bfns(j)=1
                   endif
  300           continue
  200        continue
             count=0
             do 400 i=1,nbf
                if (bfns(i).ne.0) then 
                    count=count+1
                    bfns(count)=i
                endif
  400        continue
             if (count.ne.0) then
                 call iosys('write integer "symaos-'//itoc(irrep)//
     1                      '" to rwf',count,bfns,0,' ')
                 write (iout,2) irrep, (bfns(ii),ii=1,count)
             endif         
             offset=offset+numso(irrep)
         endif
  100 continue
      write (iout,3)
      return
c
    1 format (//,10x,'atomic orbital symmetry information')
    2 format (//,5x,'irreducible representation',1x,i1,1x,
     1         /,5x,'orbitals',( (/,5x,10(i4,1x))))
    3 format (///)
      end
