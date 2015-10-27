*deck @(#)eig940.f	1.2  7/30/91
      subroutine eig940(n,nroots,eigval,eigvec,diag,repcor,prtflg,
     #                  calc,code)
c
      implicit integer (a-z)
c
      character*(*) prtflg,calc,code
      real*8 eigval(n),eigvec(n,n)
      real*8 t,eci
      real*8 repcor
      real*8 diag(n)
      character*4 itoc
c
      common /io/ inp,iout
      call pakstr(code,ilen)
c
      if(calc.ne.'mcscf') then
         write(iout,*)'      incore diagonalization'
      end if
c
      if (prtflg.ne.'minimum') then
         write (iout,90)
 90      format (/,5x,' root   reference  guess energy    ci energy  ',
     #        '  c(0)')
      end if
c
      do 10 iroot=1,nroots
         t=0.0d+00
         do 1 i=1,n
            if (abs(eigvec(i,iroot)).gt.t) then
               t=abs(eigvec(i,iroot))
               ref=i
            end if
    1    continue
         if (prtflg.ne.'minimum') then
            write (iout,2) iroot,ref,diag(ref)+repcor,
     $           eigval(iroot)+repcor,eigvec(ref,iroot)
    2       format (5x,i3,i10,2g18.9,f8.4)
         end if
         if(calc.eq.'mcscf') then
            call iosys('write real "mc root '//itoc(iroot)//'" to rwf',
     #                  n,eigvec(1,iroot),0,' ')
         else
            call iosys('write real "'//code(1:ilen)//' root '
     1                              //itoc(iroot)
     2                              //'" to rwf',n,eigvec(1,iroot),
     3                                   0,' ')
         end if
         eci=eigval(iroot)+repcor
         call iosys('write real "'//code(1:ilen)//' energy '
     1                           //itoc(iroot)
     2                           //'" to rwf',1,eci,0,' ')
 10   continue
c
      if(calc.ne.'mcscf') then
         write(iout,*)'      ci energies have been stored on rwf file'
      end if
      call iosys('write real energy to rwf',1,eigval(nroots)+repcor,
     $            0,' ')
c
      return
      end




