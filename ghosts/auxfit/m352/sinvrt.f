*deck @(#)sinvrt.f	1.1  11/20/92
      subroutine sinvrt(s,sm1,u,eigval,t1,t2,num,nnp,triang,
     $                  iprint)
c***begin prologue     sinvrt
c***date written       850601  (yymmdd)
c***revision date      921101  (yymmdd)
c
c   2 March 1993,      russo at lanl
c      modify to throw away linear combinations with small eigenvalues instead
c      of just puking
c   1 november 1992,   rlm at lanl
c      modifying sinv to compute inverse of overlap matrix
c
c***keywords           matrix, invert, square root
c***author             saxe, paul (lanl)
c***source
c***purpose                                              -1
c                      vectorized matrix inversion  v=w     .
c***description
c                      call sinvrt(s,sm1,u,eigval,t1,t2,num,nnp,
c                                triang,iprint)
c
c                        s       input matrix (nnp).
c                        sm1     output matrix (nnp)
c                        u       scratch (num,num).
c                        eigval  eigenvalues of s (num).
c                        t1      scratch (num,num).
c                        t2      scratch (num,num).
c                        num     matrix dimension.
c                        nnp     num*(num+1)/2
c                        triang  scratch (nnp).
c                        iprint  print flag:  2**0    print overlap matrix.
c                                             2**1    print eigenvectors of s.
c                                             2**2    print s**-1.
c
c***references
c***routines called    print(math),vmove(math),rsp(clams),
c                      vecout(math),vsqrt(math),zero(math),ebct(math),
c                      ebc(math),sqtotr(math),vinv(math)
c***end prologue       sinvrt
      implicit integer (a-z)
c
c     ----- input arrays (unmodified) -----
      real*8 s(nnp)
      integer num,nnp,iprint
c
c     ----- output arrays -----
      real*8 sm1(nnp)
c
c     ----- scratch arrays -----
      real*8 u(num,num),t1(num,num)
      real*8 t2(num,num),eigval(num),triang(nnp)
c
c     ----- local variables -----
      real*8 small,one
      logical debug
c
      
      parameter (small=1.0d-06)
      parameter (one=1.0d+00)
      parameter (debug=.true.)
c
      common /io/     inp,iout
c
c     ----- print the overlap matrix if requested -----
c
      if (iprint.eq.1) then
         write (iout,101)
 101     format (//,' the a.o. overlap matrix ',/)
         call print (s,nnp,num,iout)
      end if
c
c     ----- diagonalize s and get eigenvalues and vectors -----
c
      call vmove(triang,s,nnp)
      call degrsp(num,nnp,triang,eigval,1,u,t1,t2)
c
c     ----- print eigenvectors and eigenvalues of s if requested -----
c
      if (iprint.eq.4) then
         write (iout,102)
 102     format (//,' eigenvectors and values of the overlap matrix',/)
         call vecout(u,eigval,num)
       end if
c
c     ----- test for linear dependence and invert -----
      call rzero(t1,num**2)
      do 10 i=1,num
         if(abs(eigval(i)).le.small) then
            if (debug) then
               write(iout,*) 'i,eigval',i,eigval(i)
               write(iout,*) 'sinvrt: overlap matrix nearly singular'
            endif
            t1(i,i)=0.0
         else
            t1(i,i)=1.0/eigval(i)
         endif
   10 continue
c
c
c$$$c
c$$$c     ----- form s**(-1) in the basis which diagonalizes s-----
c$$$      call vinv(eigval,eigval,num)
c$$$      call rzero(t1,num**2)
c$$$      do 11 i=1,num
c$$$         t1(i,i)=eigval(i)
c$$$   11 continue
c
c     ----- transform back into original basis -----
      call ebct(t2,t1,u,num,num,num)
      call ebc(t1,u,t2,num,num,num)
      call sqtotr(sm1,t1,num,nnp)
c
c     ----- print s**(-1) if requested -----
c
      if (iprint.eq.4) then
         write (iout,103)
 103     format (//,' s**(-1) matrix ',/)
         call print (sm1,nnp,num,iout)
         call trtosq(t2,s,num,nnp)
         call ebc(u,t1,t2,num,num,num)
         write(iout,*)"inverse*s="
         do 1666 i=1,num
            write(iout,*) (u(i,j),j=1,num)
 1666    continue 
      end if
c
c
      return
      end
