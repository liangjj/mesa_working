*deck @(#)nsalc.f	1.2  7/30/91
      subroutine nsalc(salco,salcn,pkindx,nsym,mosmn,mosmo,oldnum,num,
     1                 prn)
c***begin prologue     nsalc
c***date written       910103   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m604, link 604
c***author             schneider barry (lanl)
c***source             m604
c***purpose            to form a linearly indpendent set of symmetry
c***                   orbitals in a reduced basis.
c***
c***routines called
c***end prologue       m604
      implicit integer(a-z)
      real*8 salco, salcn
      logical prn
      character *3 ans
      dimension salco(oldnum,oldnum), salcn(num,num), pkindx(oldnum)
      dimension mosmn(nsym), mosmo(nsym), prn(2)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c           bring in full salc transformation matrix                   c
c----------------------------------------------------------------------c
      call iosys('read integer "number of symmetry orbitals" from rwf',
     1             nsym,mosmo,0,' ')
      call iosys('read real "salc transformation matrix" from rwf',
     1            oldnum*oldnum,salco,0,' ')
      if (prn(1)) then
          write(iout,*) '     old salc transformation matrix'
          call matout(salco,oldnum,oldnum,oldnum,oldnum,iout)
      endif
c----------------------------------------------------------------------c
c                read packing vector                                   c
c----------------------------------------------------------------------c
      call iosys ('does "packing index vector" exist on rwf',0,0,0,ans)
      if (ans.ne.'no') then
          call iosys('read integer "packing index vector" from rwf',
     1                oldnum,pkindx,0,' ')
      else
          do 5 i=1,oldnum
             pkindx(i)=i
    5     continue
      endif
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
      if (prn(2)) then
          write(iout,*) '     new salc transformation matrix'
          call matout(salcn,num,num,num,num,iout)
      endif
      call iosys('write integer "number of new symmetry orbitals" to '
     1            //'rwf',nsym,mosmn,0,' ')
      call iosys('write real "new salc transformation matrix" to rwf',
     1            num*num,salcn,0,' ')
      return
      end
