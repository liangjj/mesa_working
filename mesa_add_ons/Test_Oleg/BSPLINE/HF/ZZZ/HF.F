*     ------------------------------------------------------------------
*
*     A Spline Hartree-Fock Program
*
*     by C. Froese Fischer
*        Vanderbilt University
*        Nashville, Tn 37235 USA
*
*     December, 1992
*
*     ------------------------------------------------------------------
*
      PROGRAM HF
*
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
      CHARACTER config*50,el*3,atom*6,term*6
      COMMON /LABEL/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
      LOGICAL fail,omit,rel,all,trace,print,strong
      COMMON /TEST/fail,omit,rel,all,trace
*
      CHARACTER ans*1, aster
      LOGICAL old
      data aster/'*'/
*
          iscw = 0
*
*
*  *****  write out header
 
      write(6,9)
9     format(///////22x,'=============================',
     :             /22x,' H a r t r e e - F o c k . 92',
     :             /22x,'=============================')
*
*  *****  write out dimension information
*
      write(6,99) 'nwf',nwd,'nt',nt
99    format(//15x,'The dimensions for the current version are:'/13x,
     :       2(10x,2(a6,'=',i3,4x)/)/)
*
*  *****  initialize
*
      call init
*
*  ***** set unit numbers and open files
*
1     write(6,'(//a/a//)') ' Start of Case',' ============='
      inquire(file='wfn.inp',exist=old)
      if (old) then
         iuf =  21
         open(unit=iuf,file='wfn.inp',status='old',
     :        form='unformatted')
      else
         iuf = 0
      end if
      ouf = 31
      open(unit=ouf,file='wfn.out',status='unknown',form='unformatted')
      open(unit=3,file='hf.log',status='unknown')
*
      fail = .false.
      do 4 i=1,(nwd)
      dpm(i) = 10.d0
      do 4 j=1,(nwd)
      e(i,j) = 0.d0
4     continue
*
*  *****  determine data about the problem
*
      call data
*
*  *****  set parameters to their default value
*
13    print = .false.
      scftol = 1.d-12
      nscf = 12
      ic = 2 + (nwf+1-ib)/4
      trace = .true.
      if (ib .le. nwf) then
         write(iscw,'(/a)')
     :            ' Default values for remaining parameters? (y/n/h) '
         read(5,'(a)') ans
         if (ans .eq. 'h' .or. ans .eq. 'h') then
            call help(4)
            go to 13
         end if
         if (ans .ne. 'y' .and. ans .ne. 'y') then
*
*  *****  additional parameters
*
   50       write(iscw,'(/a)') ' Default values (no,strong) ? (y/n/h) '
            read(5,'(a)') ans
            if (ans .eq. 'h' .or. ans .eq. 'h') then
               call help(3)
               go to 50
            end if
            if (ans .ne. 'y' .and. ans .ne. 'y') then
               write(iscw,*) ' Enter values in format(i3,1x,l1) '
               read(5,'(i3,1x,l1)') no, strong
               omit = .not. strong
            end if
 16         write(iscw,*)' Default values for print, scftol? (y/n/h)'
            read(5,'(a)') ans
            if ( ans .eq. 'h' .or. ans .eq. 'h' ) then
               call help(5)
               go to 16
            endif
            if ( ans .ne. 'y' .and. ans .ne. 'y'  ) then
               write(iscw,'(a)') ' Input FORMAT(l1, 1x, e6.1) '
               read(5,'(l1,1x,e6.1)') print, scftol
            end if
17          write(iscw,*) ' Default values for nscf, ic ? (y/n/h) '
            read(5,'(a)') ans
            if (ans .eq. 'h' .or. ans .eq. 'h' ) then
               call help(6)
               go to 17
            end if
            if (ans .ne. 'y' .and. ans .ne. 'y' ) then
               write(iscw,*) ' input format(i2, 1x, i1) '
               read(5,'(i2,1x,i1)') nscf, ic
            end if
18          write(iscw,*) ' Default values for trace ? (y/n/h) '
            read(5,'(a)') ans
            if (ans .eq. 'h' .or. ans .eq. 'h' ) then
               call help(7)
               go to 18
            end if
            if (ans .eq. 'n' .or. ans .eq. 'n') trace = .true.
         end if
      end if
*
*
*  *****  perform the mchf iteration
*
      call scf(etotal,scftol,erel)
*
*  *****  output results if print = .true.
*
      call output(print)
15    if (fail) go to 6
      call summry(etotal,erel)
19    write(iscw,'(/a)') ' additional parameters ? (y/n/h) '
      read(5,'(a)') ans
      if (ans .eq. 'h' .or. ans .eq. 'h') then
         call help(8)
         go to 19
      end if
      if (ans .eq. 'y' .or. ans .eq. 'y') call menu
*
*  *****  check for isoelectronic sequence or end of case.
*
20    write(iscw,'(/a)')' do you wish to continue along the sequence ? '
      read(5,'(a1)') ans
      if (ans .eq. 'h' .or. ans .eq. 'h') then
         call help(9)
         go to 20
      end if
*
*  *****  determine end of case
*
6     write(6,'(//a/a//)') ' end of case',' ==========='
      stop
      end
