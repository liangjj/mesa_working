*deck  @(#)nsalc.f	5.1 11/6/94
      subroutine nsalc(salco,salcn,pkindx,nsym,mosmo,mosmn,oldnum,num,
     1                 prnt)
c***begin prologue     nsalc
c***date written       910103   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m606, link 606
c***author             schneider barry (lanl)
c***source             @(#)nsalc.f	5.1 11/6/94
c***purpose            to form a linearly indpendent set of symmetry
c***                   orbitals in a reduced basis.
c***
c***routines called
c***end prologue       m606
      implicit integer(a-z)
      logical prnt
      integer pkindx(oldnum),mosmo(oldnum), mosmn(num)
      real*8 salco(oldnum,oldnum), salcn(num,num)
c
      common /io/ inp, iout
c
c----------------------------------------------------------------------c
c           print full salc transformation matrix if desired           c
c----------------------------------------------------------------------c
      if (prnt) then
          write(iout,*) 'old salc transformation matrix'
          call matout(salco,oldnum,oldnum,oldnum,oldnum,iout)
      endif
c----------------------------------------------------------------------c
c                read packing vector                                   c
c----------------------------------------------------------------------c
      call iosys('read integer "packing index vector" from rwf',
     1                oldnum,pkindx,0,' ')
c----------------------------------------------------------------------c
c                check for zero vector and toss out salc               c
c----------------------------------------------------------------------c
      newtot=0
      count=0
      do 10 i=1,nsym
         mocnt=0
         do 20 j=1,mosmo(i)
            count=count+1
            sum=0.d+00
            do 30 k=1,oldnum
               sum=sum+abs(salco(k,count))*pkindx(k)
   30       continue
            if (sum.ne.0.d+00) then
                mocnt=mocnt+1
                newtot=newtot+1
                call rzero(salcn(1,newtot),num)
                do 40 k=1,oldnum
                   if (pkindx(k).ne.0) then
                       salcn(pkindx(k),newtot)=salco(k,count)
                   endif
   40           continue
            endif
   20    continue
         mosmn(i)=mocnt
   10 continue
      if (newtot.ne.num) call lnkerr('error in determining new '//
     1                               'symmetry orbitals')
      if (prnt) then
          write(iout,*) 'new salc transformation matrix'
          call matout(salcn,num,num,num,num,iout)
      endif
c
c
      return
      end
