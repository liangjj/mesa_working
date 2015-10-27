*deck order.f
c***begin prologue     order
c***date written       930612   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           order
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            order integers or reals in ascending order
c***description        
c***references         none
c
c***routines called
c***end prologue       order
      subroutine order(vals,ivals,n,type)
      implicit integer (a-z)
      real*8 vals, p
      character*(*) type
      dimension vals(n), ivals(n)
c                rearrange a vector of real or integer
c                values in increasing order
      if (type.eq.'real') then
          do 10 ii=2,n
             i=ii-1
             k=i
             p=vals(i)
             do 20 j=ii,n
                if (vals(j).lt.p) then
                    k=j
                    p=vals(j)
                endif
   20        continue
             if (k.ne.i) then
                 vals(k)=vals(i)
                 vals(i)=p
             endif
   10     continue
      elseif (type.eq.'integer') then
          do 30 ii=2,n
             i=ii-1
             k=i
             ip=ivals(i)
             do 40 j=ii,n
                if (ivals(j).lt.ip) then
                    k=j
                    ip=ivals(j)
                endif
   40        continue
             if (k.ne.i) then
                 ivals(k)=ivals(i)
                 ivals(i)=ip
             endif
   30      continue
      else
           call lnkerr('error in call to order of variable type')
      endif
      return
      end

