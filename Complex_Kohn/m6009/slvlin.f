*deck @(#)slvlin.f	1.1 9/8/91
c***begin prologue     slvlin
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           slvlin, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            solve complex linear systems
c***description        solve a set(s) of linear equations in which
c***                   the matrix is real or complex and the 
c***                   right hand side is real or complex.
c***                  
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       slvlin
      subroutine slvlin (hambb,hambbc,ipvt,hambfr,hambfc,n,m,ifac,
     1                   typem,typerh,prnt)
      implicit integer(a-z)
      real*8 hambb, hambfr
      complex*16  hambbc, hambfc
      character *(*) ifac, typem, typerh
      character *80 title
      logical prnt
      dimension hambb(n,n), hambbc(n,n), hambfr(n,m), hambfc(n,m)
      dimension ipvt(n)
      common /io/ inp, iout
c----------------------------------------------------------------------c
c     compute lu factorization of bound bound matrix                   c
c----------------------------------------------------------------------c
      if (ifac.eq.'factor') then
          if (typem.eq.'real') then
              call sgefa (hambb,n,n,ipvt,info)
          elseif (typem.eq.'complex') then
              call cgefa(hambbc,n,n,ipvt,info)
          endif
      else
c----------------------------------------------------------------------c
c                solve equations                                       c
c----------------------------------------------------------------------c
         if (typem.eq.'real') then
             if (typerh.eq.'real') then
                 do 20 i=1,m
                    call sgesl (hambb,n,n,ipvt,hambfr(1,i),0)
   20            continue
                 if (prnt) then
                     title='real right hand side'
                     call prntrm(title,hambfr,n,m,n,m,iout)
                 endif
             elseif(typerh.eq.'complex') then
                do 30 i=1,m
                   call rgslc(hambb,n,n,ipvt,hambfc(1,i))
   30           continue
                if (prnt) then
                    title='complex right hand side'
                    call prntcm(title,hambfc,n,m,n,m,iout)
                endif
             endif
         elseif (typem.eq.'complex') then
                 do 60 i=1,m
                    call cgesl (hambbc,n,n,ipvt,hambfc(1,i),0)
   60            continue
                 if (prnt) then
                     title='complex right hand side'
                     call prntcm(title,hambfc,n,m,n,m,iout)
                 endif
         endif
      endif
      return
      end


