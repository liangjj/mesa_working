*deck deabm
      subroutine deabm(t,psi,psi0,t0,tf,delt,info,eps1,eps2,
     1                 neq,nstp,h01,h02,h03,eig1,eig2,eig3,
     2                 v,vt,tpert,efield,omega,width,shift,
     3                 hbar,scr,iscr,dim,nmax,n)           
c***begin prologue  deabm
c***purpose   solve a set of coupled first-order differential equations
c             equations using the adams-bashforth-moulton method.
c***keywords  adams-bashforth method, depac, initial value problems,
c             ode, ordinary differential equations, predictor-corrector
c***author    shampine, l. f., (snla), watts, h. a., (snla)
c***description this is a modified version of the code written by the above
c***            authors.  the description can be found in excruciating
c***            detail in the readme file in diffeqlib.  
c
c***references  l. f. shampine and h. a. watts, depac - design of a user
c                 oriented package of ode solvers, report sand79-2374,
c                 sandia laboratories, 1979.
c***routines called  ddes, xermsg
c***revision history  (yymmdd)
c***end prologue  ddeabm
c
      implicit integer(a-z)
      real*8 t, psi, psi0, eps1, eps2, rdum, scr
      real*8 t0, tf, delt, tin, tout, ttemp
      logical start, phase1, nornd, stiff, intout
c
      dimension psi0(neq), info(15), eps1(*), eps2(*)
      dimension t(*), psi(neq,nstp)
      dimension iloc(24), scr(*), iscr(*)
      real*8 h01, h02, h03, eig1, eig2, eig3, v, vt
      real*8 omega, efield, width, shift, hbar
      character*80 title
      character*16 fptoc      
      character*(*) tpert
      dimension nmax(3)
      dimension h01(nmax(1),nmax(1)), h02(nmax(2),nmax(2))
      dimension h03(nmax(3),nmax(3))
      dimension eig1(nmax(1)), eig2(nmax(2)), eig3(nmax(3))
      dimension v(n), vt(n)
      common/io/inp,iout
      call copy(psi0,psi(1,1),neq)
c
c
c
c     check for an apparent infinite loop
c
c***first executable statement  ddeabm
c
      iwds=51
c
c     initialize the first value of x, the starting vector and info
c
      info(1) = 0
      info(2) = 0
      info(3) = 0
      info(4) = 0      
      iscr(iwds) = 0
c      
c     compute the indices for the arrays to be stored in the work array
c

      iloc(1) = 21                     
      iloc(2) = iloc(1) + neq           
      iloc(3) = iloc(2) + 1            
      iloc(4) = iloc(3) + neq          
      iloc(5) = iloc(4) + neq          
      iloc(6) = iloc(5) + neq          
      iloc(7) = iloc(6) + neq          
      iloc(8) = iloc(7) + neq*16       
      iloc(9) = iloc(8) + 12           
      iloc(10) = iloc(9) + 12          
      iloc(11) = iloc(10) + 12         
      iloc(12) = iloc(11) + 12         
      iloc(13) = iloc(12) + 12         
      iloc(14) = iloc(13) + 13         
      iloc(15) = iloc(14) + 13         
      iloc(16) = iloc(15) + 11         
      iloc(17) = iloc(16) + 1          
      iloc(18) = iloc(17) + 1          
      iloc(19) = iloc(18) + 1          
      iloc(20) = iloc(19) + 1           
      iloc(21) = iloc(20) + 1          
c
      tin = t(1)
      do 100 i=1,nstp
         tout = t(i+1)
         scr(iloc(2)) = tin      
         if (iscr(iwds) .ge. 5) then
             if (tin .eq. scr(21 + neq)) then
                 call xermsg ('diffeqlib', 'ddeabm',
     *                        'an apparent infinite loop has '//
     *                        'been detected.$$' //
     *                        'you have made repeated calls '//
     *                        'at t = ' // fptoc(t) //
     *                        ' and the integration has not '//
     *                        'advanced.  check the ' //
     *                        'way you have set parameters for '//
     *                        'the call to the ' //
     *                        'code, particularly info(1).', 13, 2)
                 idid=-33
                 call lnkerr('quit')         
                 return
             endif
         endif
c
         if (info(1).ne.0) then
             start = iscr(21) .ne. (-1)
             phase1 = iscr(22) .ne. (-1)
             nornd = iscr(23) .ne. (-1)
             stiff = iscr(24) .ne. (-1)
             intout = iscr(25) .ne. (-1)
         endif   
c
         ttemp=tin
         call des(neq,tin,psi0,tout,info,eps1,eps2,idid,
     1            scr(iloc(1)),scr(iloc(3)),scr(iloc(4)),scr(iloc(5)),
     2            scr(iloc(6)),scr(iloc(7)),scr(iloc(8)),
     3            scr(iloc(9)),scr(iloc(10)),scr(iloc(11)),
     4            scr(iloc(12)),scr(iloc(13)),scr(iloc(14)),
     5            scr(iloc(15)),scr(11),scr(12),scr(13),
     6            scr(iloc(16)),scr(iloc(17)),scr(iloc(18)),
     7            scr(iloc(19)),scr(1),scr(iloc(20)),
     8            scr(iloc(21)),start,phase1,nornd,stiff,
     9            intout,iscr(26),iscr(27),iscr(28),iscr(29),
     x            iscr(30),iscr(31),iscr(32),iscr(33),iscr(34),
     x            iscr(35),iscr(45),rdum,ipar,
     x            h01,h02,h03,eig1,eig2,eig3,v,vt,tpert,efield,omega,
     x            width,shift,hbar,dim,nmax,n)
         write(iout,1) iscr(30)
c
         iscr(21) = -1
         if (start) then
             iscr(21) = 1
         endif             
         iscr(22) = -1
         if (phase1) then
             iscr(22) = 1
         endif    
         iscr(23) = -1
         if (nornd) then
             iscr(23) = 1
         endif    
         iscr(24) = -1
         if (stiff) then
             iscr(24) = 1
         endif   
         iscr(25) = -1
         if (intout) then
             iscr(25) = 1
         endif    
c
         if (idid .ne. (-2)) then
             iscr(iwds) = iscr(iwds) + 1
         endif    
         if (tin .ne. scr(iloc(2))) then
             iscr(iwds) = 0
         endif
         if(idid.lt.0) then
            write(iout,2)
            call lnkerr('quit')
         else
              call copy(psi0,psi(1,i+1),neq)
              title='psi for t(in) = '//fptoc(ttemp)//' t(out) = '
     1                               //fptoc(tout)
              call prntrm(title,psi0,neq,1,neq,1,iout)
         endif               
 100  continue   
c
      return
 1    format(/,1x,'number of steps taken = ',i3)
 2    format(/,1x,'idid = ',i4,' cannot continue')
      end
