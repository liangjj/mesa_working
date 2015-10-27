*     ------------------------------------------------------------------
*               m e n u
*     ------------------------------------------------------------------
*
*
*     this routine evaluates a variety of atomic parameters as
*     requested by the user.
*
*
      SUBROUTINE menu
     
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
      INTEGER ouf
      COMMON /inout/ iuf,ouf,iscw
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
     
      character*3 el1,el2,el3,el4,func*1,funct*1
     
4     write(iscw,5)
5     format(//5x,
     :  ' these various functions are available:',//
     :  10x,'1 - expectation values of r**k'/
     :  10x,'2 - slater or magnetic integrals'/
     :  10x,'3 - electron density at the nucleus'/
     :  10x,'4 - spin-orbit parameter'/
     :  10x,'5 - transition integrals'/
     :  10x,'6 - exit to main program'/)
      write(iscw,'(5x,a)')'input number 
     :  corresponding to your selection:'
      read(5,'(i1)') ifunc
      go to (10,20,30,40,50,60) ifunc
     
*  ****  compute expectation values
     
10    write(iscw,'(/5x,a,/a,t22,a)')
     :  'input label for electron followed by k: example',
     :  '   2p  3','format(1x,a3,i3) '
      read(5,'(1x,a3,i3)') el1,k
      call eptr(el,el1,i,*10)
      rkev = quadr(i,i,k)
      write(3,12) el1,k,el1,rkev
      write(iscw,12) el1,k,el1,rkev
12    format(/15x,' value of <',a3,'|r**',i2,'|',a3,'> = ',1pd14.7,
     :         ' a.u.'/)
      go to 4
     
*  ****  determine slater integrals  fk, gk, rk, nk, mk, vk
     
     
20    write(iscw,'(/5x,a/a,t22,a)')
     : 'input parameters for  fk,gk,rk,nk,mk or vk integral: example',
     : ' f 0( 1s, 2s)','format:  (a1,i2,1x,4(a3,1x)) '
      read(5,'(a1,i2,1x,4(a3,1x))') funct,k,el1,el2,el3,el4
      if (funct .ge. 'a' .and. funct .le. 'z')
     :   func = char(ichar(funct) + ichar('a') - ichar('a'))
      call eptr(el,el1,i1,*20)
      call eptr(el,el2,i2,*20)
      if ( el3 .ne. ' ') call eptr(el,el3,i3,*20)
      if ( el4 .ne. ' ') call eptr(el,el4,i4,*20)
      if (func .eq. 'f' ) then
        si = fky(i1,i2,k,rel)
      else if (func .eq. 'g') then
        si = gky(i1,i2,k,rel)
      else if (func .eq. 'r') then
        si = rky(i1,i2,i3,i4,k,rel)
      else if (func .eq. 'n') then
*       si = sn(i1,i2,i2,i1,k)
      else if (func .eq. 'm') then
*       si = sn(i1,i2,i1,i2,k)
      else if (func .eq. 'v') then
*       si = vk(i1,i2,i2,i1,k) - vk(i2,i1,i1,i2,k)
      else
        write(iscw,41)
41      format(15x,'integral unknown: re-enter')
        go to 20
      end if
      if (func .ne. 'r') then
         write(3,25) func,k,el1,el2,si,219474.d0*si
         write(iscw,25) func,k,el1,el2,si,219474.d0*si
25       format(/15x,
     :  'integral  ',a1,i2,'(',a3,',',a3,') = ',1pd14.7,' a.u.'/
     :   40x,0pf14.3,' cm-1'/)
      else
         write(3,26) func,k,el1,el2,el3,el4,si,219474.d0*si
         write(iscw,26) func,k,el1,el2,el3,el4,si,219474.d0*si
26       format(/15x,
     :  'integral  ',a1,i2,'(',2a3,',',2a3,') = ',1pd14.7,' a.u.'/
     :   46x,0pf14.3,' cm-1'/)
      end if
      go to 4
     
*  ****  compute electron density at the nucleus
     
30    write(iscw,'(/5x,a/a,t22,a)')
     :  'input identifying label for electron: example',
     :  '   1s','format(1x,a3) '
      read(5,'(1x,a3)') el1
      call eptr (el,el1,i,*30)
      ll = l(i)
      if (ll .eq. 0) then
         d = az(z,h,ks,l(i)+1,p(l(i)+2,i))**2
      else
         d = 0
      end if
      write(3,32) el1,d
      write(iscw,32) el1,d
32    format(/15x,'density at the nucleus for ',a3,' = ',1pd14.7,
     :     ' a.u.'/)
      go to 4
     
*  ****  compute spin-orbit parameter
     
40    write(iscw,'(/,5x,a/a,t22,a)')
     :  'input identifying label for electron: example',
     :  '   2p','format(1x,a3) '
      read(5,'(1x,a3)') el1
      call eptr (el,el1,i,*40)
      zeta = 0.d0
      if (l(i) .ne. 0)  zeta = bwzeta(i)
      zetacm = 219474*zeta
      write(3,43) el1,zeta,zetacm
      write(iscw,43) el1,zeta,zetacm
43    format(/15x,'spin-orbit parameter for ',a3,' = ',1pd14.7,
     :        ' a.u.'/46x,0pf14.3,' cm-1'/)
      go to 4
     
*  ****  compute transition integrals
     
50    write(iscw,'(/5x,a/a,t22,a)')
     :  'input identifying labels and power of r: example',
     :  ' t 1( 2s, 2p)','format:  (a1,i2,2(1x,a3)) '
      read(5,'(a1,i2,1x,a3,1x,a3)') func,k,el1,el2
      call eptr (el,el1,i1,*50)
      call eptr (el,el2,i2,*50)
      ti = quadr(i1,i2,k)
      write(3,52) func,k,el1,el2,ti
      write(iscw,52) func,k,el1,el2,ti
52    format(/15x,
     :  'integral  ',a1,i2,'(',a3,',',a3,') = ',1pd14.7,' a.u.'/)
      go to 4
60    return
      end
