*deck %W% %G%
      subroutine pm999(z,a)
c
c***begin prologue     %M%
c***date written       940104  (yymmdd)
c***revision date      %G%
c
c***keywords           least squares
c***author             russo, thomas (lanl)
c***source             %W% %G%
c***purpose            solves the least squares problem
c***description
c
c
c***references
c                      memory
c***routines called
c***end prologue       %M%
c
      implicit none
      integer a(*)
      real*8 z(*)
      integer maxdat,thisdat,canget,mat1,mat2,vec1,vec2,ipvt,top,i,
     $     terms
      integer wpadti
      parameter(maxdat=128)
      integer polord,pos,start,end
      real*8 x(maxdat),y(maxdat),f(maxdat)
      integer inp,iout
      integer intkey
      real*8 ctofp
      character*4096 ops
      character*80 card
      character*80 token
      character*16 ffnext
      logical positn
      common /io/ inp,iout
c
c read xy pairs from inp
c 
      call iosys('read character options from rwf',-1,0,0,ops)
      polord=intkey(ops,'polord',3,' ')
      if (.not.positn('$data',card,inp)) then
         call lnkerr('xm666: no $data section found')
      endif
      thisdat=1
 1    continue 
      read (inp,1000,end=2) card
 1000 format (a)
      if (card .eq. ' ' .or. index(card,'$').ne.0) goto 2
      pos=0
      token=ffnext(card,pos,start,end)
      if (token .ne. 'floating point') then
         call lnkerr('data must be floating point')
      endif
      x(thisdat)=ctofp(card(start:end))
      token=ffnext(card,pos,start,end)
      if (token .ne. 'floating point') then
         call lnkerr('data must be floating point')
      endif
      y(thisdat)=ctofp(card(start:end))
      token=ffnext(card,pos,start,end)
      if (token .ne. 'floating point') then
         call lnkerr('data must be floating point')
      endif
      f(thisdat)=ctofp(card(start:end))
      thisdat=thisdat+1
      goto 1
 2    continue 
      thisdat=thisdat-1

      write (iout,*)' input data:'
      do 6 i=1,thisdat
         write(iout,*)x(i),f(i)
 6    continue 
c allocate core
      call getscm(0,z,canget,'xm999:',0)

c
c we need to fit to a polord surface, which means we have 
c (polord*(polord+1))/2+polord+1 terms in the fit.
c
      terms=(polord*(polord+1))/2 + polord+1
      mat1=1
      mat2=mat1+thisdat*terms
      vec1=mat2+terms**2
      vec2=vec1+terms
      ipvt=wpadti(vec2+terms)
      top=ipvt+terms

      if (top .gt. canget) then
         call lnkerr('xm666: wow, that is a big problem.')
      endif

      call dolsq(thisdat,polord,x,y,f,z(mat1),z(mat2),z(vec1),z(vec2),
     $     a(ipvt),terms)

      stop
      end
