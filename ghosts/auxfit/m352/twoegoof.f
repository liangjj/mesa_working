*deck %W% %G%
      subroutine twoegoof(fna0,tmpsiz,lval,mycart,yn0,mmax,mindx)
c
c***begin prologue     %M%
c***date written       930205  
c***revision date      %G%
c
c***keywords           
c***author             russo, thomas (lanl)
c***source             %W% %G%
c***purpose            
c***description
c      calculate the arrays "fna0" and "tmpsiz" for the auxiliary two-e
c      code.  fna0(m) is the number of (a,0) integrals to calculate when
c      goes from 0 to m.
c      tmpsiz(m) is size of gena0int's temp array when doing the (m,0).  
c***references
c
c***routines called
c
c***end prologue       %M%

      integer fna0(0:mmax),tmpsiz(0:mmax),mycart(0:mmax),yn0(0:mmax)
      integer lval(0:mmax,mindx,3)
      logical debug
      parameter (debug=.false.)
      common /io/ inp,iout
c
c fna0 is just sum of (l+1)(l+2)/2 for l=0 to m
c
c mycart(m) is the number of unique triples (lx,ly,lz) with lx+ly+lz=m
c yn0(m) is the location in my ordering of the first lx=0,ly!=0 integral
c
      yn0(0)=0
      do 10 m=0,mmax
         fna0(m)=0
         mycart(m)=((m+1)*(m+2))/2
         if (m.gt.0) yn0(m)=mycart(m-1)
         do 20 l=0,m
            fna0(m)=fna0(m)+mycart(l)
 20      continue 
 10   continue 
c
c tmpsiz(m) is the max of (m-l+1)(l+1)(l+2)/2 where l goes from 0 to m
c
      do 30 m=0,mmax
         tmpsiz(m)=1
         do 40 l=0,m
            tmpsiz(m)=max(tmpsiz(m),(m-l+1)*((l+1)*(l+2))/2)
 40      continue 
 30   continue 
c
c lval(m,i,coord) is the l value of the coord'th component of the i'th integral
c ie for the (100)(010)(001) list, lval(1,1,1)=1,lval(1,1,2)=0,lval(1,1,3)=0,
c etc.
c

      do 50 m=0,mmax
         i=1
         do 51 lx=m,0,-1
            do 52 ly=m-lx,0,-1
               lz=m-lx-ly
               lval(m,i,1)=lx
               lval(m,i,2)=ly
               lval(m,i,3)=lz
               i=i+1
 52         continue 
 51      continue 
 50   continue 

      if (debug) then
         write(iout,*)"mycart:"
         write(iout,*)(mycart(i),i=0,mmax)
         write(iout,*)"yn0:"
         write(iout,*)(yn0(i),i=0,mmax)
         write(iout,*)"tmpsiz:"
         write(iout,*)(tmpsiz(i),i=0,mmax)
         write(iout,*)"fna0:"
         write(iout,*)(fna0(i),i=0,mmax)
         write(iout,*)"lval:"
         do 60 m=0,mmax
            write(iout,*)"m=",m,(lval(m,i,1),i=1,mycart(m))
            write(iout,*)"     ",(lval(m,i,2),i=1,mycart(m))
            write(iout,*)"     ",(lval(m,i,3),i=1,mycart(m))
 60      continue 
      endif
      return
      end
