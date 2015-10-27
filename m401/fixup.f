*deck @(#)fixup.f	5.1  11/6/94
c***begin prologue     $vectors
c***date written       850810   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m401, link 401, guess, $vectors, rdinp
c***author             martin, richard    (lanl)
c***source             @(#)fixup.f	5.1   11/6/94
c***purpose            reads scf guess vectors from the input stream.
c***description
c     the initial guess vectors may be read in a list-directed
c     stream from the input file.  only the occupied vectors are expected.
c***references
c***routines called
c***end prologue       $vectors
      subroutine fixup(n,u,u0,t,temp)
c
      implicit integer (a-z)
c
      real*8 u(n,n),u0(n,n),t(n,n),temp(n,n),large
c
      common /io/ inp,iout
c
      call ebtc(t,u,u0,n,n,n)
c
c      write (iout,90)
c   90 format (' u(t) . u0')
c      call matout(t,n,n,n,n,iout)
c
c     ----- reorder u so that u(t) . uo is close to unit matrix -----
c            phase is also important !
c
      call rzero(temp,n**2)
      do 20 j=1,n
         large=abs(t(j,1))
         column=1
         do 19 i=2,n
            if (abs(t(j,i)).gt.large) then
               large=abs(t(j,i))
               column=i
            end if
   19    continue
         if (large.lt.0.8) call lnkerr('u(t) . u0 is not diagonal')
         if (t(j,column).gt.0.0d+00) then
            do 17 i=1,n
               temp(i,j)=u(i,column)
   17       continue
         else
            do 16 i=1,n
               temp(i,j)=-u(i,column)
   16       continue
         end if
   20 continue
c
      call vmove(u,temp,n**2)
c
      call ebtc(t,u,u0,n,n,n)
      do 30 j=1,n
         if(abs(t(j,j)).lt.0.707) call lnkerr('problems with u(t) . u0')
         if (t(j,j).lt.0.0d+00) then
            do 29 i=1,n
               u(i,j)=-u(i,j)
   29       continue
         end if
   30 continue
c
c     ----- check the results -----
c
      call ebtc(t,u,u0,n,n,n)
c      write (iout,91)
c   91 format (' modified u(t) . u0')
c      call matout(t,n,n,n,n,iout)
c
c
      return
      end
