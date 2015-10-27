*deck quad
      subroutine quad(p,w,npts,lmax)
      implicit none
      real*8 p(3,npts),w(npts)
      integer npts,lmax
      real*8 theta,phi,acos,atan2
      real*8 sr,si,diff
      complex*16 ylm,t1
      integer l,m,lp,mp,i
c
c
      do 100 l=0,lmax
         do 100 m=-l,l
            do 100 lp=0,l
               if(l+lp.gt.lmax) goto 100
               do 90 mp=-lp,lp
                  sr=0.0d+00
                  si=0.0d+00
                  do 70 i=1,npts
                     theta=acos(p(3,i))
                   if((p(1,i).eq.0.0d+00).and.(p(2,i).eq.0.0d+00)) then
                        phi=0.0d+00
                     else
                        phi=atan2(p(2,i),p(1,i))
                     endif
                     t1=conjg(ylm(l,m,theta,phi))*ylm(lp,mp,theta,phi)
                     sr=sr+w(i)*real(t1)
                     si=si+w(i)*imag(t1)
   70             continue
                  if((l.eq.lp).and.(m.eq.mp)) then
                     diff=1.0d+00-abs(sr)
                  else
                     diff=abs(sr)
                  endif
		  write(6,*)"difference in real part=",diff
		  write(6,*)"  imag part=",si
                  if((abs(diff).ge.1.0d-7).or.(abs(si).ge.1.0d-7)) then
                     write(0,*) 'error',l,m,lp,mp,sr,si
                  endif
90	continue
  100 continue
c
c
      return
      end
