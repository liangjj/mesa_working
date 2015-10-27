*deck htonv.f
c***begin prologue     htonv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           matrix vector
c***author             schneider, barry (nsf)
c***source             tproplib              
c***purpose            matrix vector multiply for time plus space
c***                   dvr hamiltonian. 
c***description        the vector vecin is transformed to the vector
c***                   vecout by operating with the hamiltonian.              
c
c                      n = n3*n2*n1*nt
c                      nc = number of channels
c                      ht(nt,nt) = matrix representation of d
c                                                           --
c                                                           dt
c                      h1(n1,n1) = matrix representation of T(1) + V0(1)
c                      h2(n2,n2) = matrix representation of T(2) + V0(2)
c                      h3(n3,n3) = matrix representation of T(3) + V0(3)
c                      vecin(n,nc,nvc)  = vecin(n3,n2,n1,nt,nvc)
c                      vecout(n,nc,nvc) = vecout(n3,n2,n1,nt,nvc)
c                      v123t(n,nc,nc)   = interaction potential, diagonal in coordinate space
c                                         possible, off-diagonal in channels
c
c                            Structure of Subroutine
c
c                      [ -H  - d/dt ] [ VinReal ]     [ VoutReal ]
c                      [            ] [         ]  =  [          ]
c                      [ d/dt   -H  ] [ VinImag ]     [ VoutImag ]
c                                                                  
c***references         
c***routines called    
c***end prologue       htonv
      subroutine htonv(h1,h2,h3,ht,v123t,vecin,vecout,n,nc,
     1                 n1,n2,n3,nt,nvc,dim)
      implicit integer (a-z)
      real*8 ht, h1, h2, h3, v123t, vecin, vecout
      dimension ht(nt,nt), h1(n1,n1), h2(n2,n2), h3(n3,n3)
      dimension v123t(n,nc,nc), vecin(n,nc,2,nvc), vecout(n,nc,2,nvc)
      common/io/inp, iout
c
      call rzero(vecout,n*nc*2*nvc)
      call vonv(vecout,v123t,vecin,n,nc,2*nvc)
c
c     now we need special cases
c
      if(dim.eq.1) then
         call h1tonv(h1,ht,vecout,vecin,n1,nt,nc,nvc)
      elseif(dim.eq.2) then
         call h2tonv(h1,h2,ht,vecout,vecin,n1,n2,nt,nc,nvc)
      elseif(dim.eq.3) then
         call h3tonv(h1,h2,h3,ht,vecout,vecin,n1,n2,n3,nt,nc,nvc)
      endif
      return
      end       



