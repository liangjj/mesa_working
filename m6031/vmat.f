*deck vmat
      subroutine vmat(vp,v,vec,s,scr,rbox,power,v1,nprim,ncon,
     1                pmax,type)
      implicit integer (a-z)
      real*8 vp, v, vec, s, scr, rbox, fpkey, v1, v2, a1, a2
      real*8 f1, f2, f3
      character*8 cpass
      character*(*) type
      character*1600 card
      logical posinp
      dimension vp(nprim,nprim), v(ncon,ncon), power(nprim)
      dimension vec(nprim,*), scr(0:pmax), s(nprim,nprim)
      common/io/ inp, iout
c     potential energy matrix elements
c
      if ( posinp('$vparam',cpass) ) then
           call cardin(card)
      endif
      if (type.eq.'square-well') then
          v1=fpkey(card,'first-depth',-5.d0,' ')
          v2=fpkey(card,'second-depth',0.d0,' ')
          a1=fpkey(card,'first-length',rmat,' ')
          a2=fpkey(card,'second-length',rmat,' ')
          write(iout,1) v1,a1,v2,a2
      elseif(type.eq.'exponential') then
          v1=fpkey(card,'exponential-prefactor',-1.d0,' ')
          a1=fpkey(card,'exponential-strength',1.d0,' ')
          write(iout,2) v1,a1
      else
          call lnkerr('error in potential type')
      endif
      if (type.eq.'square-well') then                  
          do 10 i=1,nprim
             do 20 j=1,i
                lsum=power(i)+power(j)+1
                vp(i,j)=v1*(a1**lsum/lsum)+v2*( a2**lsum/lsum -
     1                                            a1**lsum/lsum )
                vp(i,j)=vp(i,j)*s(i,i)*s(j,j)
                vp(j,i)=vp(i,j)
 20          continue   
 10       continue
      else
           f1=1.d0/a1
           f2=rbox
           f3=exp(-a1*rbox)
           scr(0)=(1.d0-f3)*f1
           do 30 n=1,pmax
              scr(n)=-f2*f3*f1 + n*f1*scr(n-1)
              f2=f2*rbox 
   30      continue
           do 40 i=1,nprim
              do 50 j=1,i
                 vp(i,j)=v1*s(i,i)*s(j,j)*scr(power(i)+power(j))
                 vp(j,i)=vp(i,j)
   50         continue
   40      continue                               
      endif
c     transform
c   
      call ebtc(scr,vec,vp,ncon,nprim,nprim)
      call ebc(v,scr,vec,ncon,nprim,ncon)
      return
    1 format(/,5x,'potential parameters',//,1x,
     1            'depth of well     = ',e15.8,2x,
     2            'length of well    = ',e15.8,/
     3            'height of barrier  = ',e15.8,2x,
     4            'length of barrier = ',e15.8)
    2 format(/,5x,'exponential potential prefactor = ',e15.8
     1       /,5x,'exponential potential strength  = ',e15.8)
      end


