*deck @(#)hukmat.f	5.1  11/6/94
      subroutine hukmat(s,smhalf,ediag,corflg,h,c,eigval,triang,u,t1,
     $                  nb,nnp)
c***begin prologue     hukmat
c***date written       850601  yymmdd
c***revision date      yymmdd  yymmdd
c***keywords           huckel, initial guess
c***author             martin, richard (lanl)
c***source             @(#)hukmat.f	5.1   11/6/94
c***purpose            forms and diagonalizes the extended huckel matrix.
c***description
c     call hukmat(s,smhalf,ediag,corflg,h,c,eigval,triang,u,t1,
c                 nb,nnp)
c***references         (none)
c***routines called    sqtotr(math), getmo(m401)
c***end prologue       hukmat
      implicit integer(a-z)
      real*8 s(nb,nb),smhalf(nnp),ediag(nb),h(nnp)
      real*8 c(nb,nb),eigval(nb),triang(nnp),u(nb,nb),t1(nb,nb)
      real*8 pt875,scale,zero
      dimension corflg(nb)
      parameter (pt875=0.875d0,zero=0.0d0)
c
c     form a huckel guess.  the diagonal elements and orbital types
c     formed by minbas are used to construct the full huckel
c     hamiltonian, which is then diagonalized.  arguments:
c
c     u
c     t1
c     h
c     triang ... scratch arrays.
c     c      ... huckel eigenvectors.
c     eigval ... huckel eigenvalues.
c     nb     ... the number of basis functions.
c     nnp    ... nb*(nb+1)/2
c     s      ... the overlap matrix.
c     smhalf ... s(-1/2).
c     ediag  ... diagonal elements from minbas.
c     corflg ... orbital types from minbas.  only the sign, which
c                distinguishes between core and valence, is used.
c
c
      scale=pt875
      do 20 i=1,nb
          u(i,i)=ediag(i)
          limj=i-1
          if(limj.gt.0) then
            do 10 j=1,limj
               u(i,j)=scale*s(i,j)*(ediag(i)+ediag(j))
               if(corflg(i).lt.0.or.corflg(j).lt.0) u(i,j)=zero
               u(j,i)=u(i,j)
   10       continue
         endif
   20 continue
      call sqtotr(h,u,nb,nnp)
c
c     get the orbitals.
      call getmo(smhalf,h,c,eigval,triang,u,t1,nb,nnp)
c
c
      return
      end
