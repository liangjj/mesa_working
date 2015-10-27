*
*     ------------------------------------------------------------------
*                       e n e x p r
*     ------------------------------------------------------------------
*
*     determine the deviations to the average energy for the following:
*        i) an open p- or d-shell
*       ii) a single electron or hole, any l
*      iii) an s-electron and a single electron, any l
*       iv) an s-electron and an open p- or d-shell
*        v) an open p-shell and a single electron, any l
*
      SUBROUTINE enexpr(term, done)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
      CHARACTER term*6,sl*2,senor*1,psl*2,slm*2,slp*2
      LOGICAL done
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      INTEGER sumtab(5),partab(11),ptrtab(11),ltab(54),
     :   nos(2),plval(11),pacval,sp,ps1,ps2
*
*     ... fint, gint1, and gint2 are coefficients of polynomials
*         in l, tabulated by slater,
*
*
      INTEGER fint(3,54),gint1(3,54),gint2(3,54)
      CHARACTER*1 parch(11)
*
*     ... coefficients of f2 integrals for p(n)l(1) configurations
*
      data fint/2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,
     :   -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3,
     :   4,-2,0,-2,-11,6,-4,-4,15,-2,7,15,4,10,6,0,0,0,
     :   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :   2,-1,0,-4,-4,3,2,5,3,2,-1,0,-4,-4,3,2,5,3,
     :   -4,2,0,2,11,-6,4,4,-15,2,-7,-15,-4,-10,-6,0,0,0,
     :   -2,1,0,4,4,-3,-2,-5,-3,-2,1,0,4,4,-3,-2,-5,-3/
*
*     ... coefficients of g(l-1) integrals
*
      data gint1/-10,5,0,2,11,-6,2,-1,-6,14,-7,0,2,-13,6,2,-1,6,
     :   -8,4,0,-8,4,0,4,10,0,10,-5,0,10,-5,0,4,-8,0,
     :   -8,4,0,-2,13,-6,2,11,-12,4,4,-12,4,-2,-12,0,0,0,-6,3,0,10,-5,0,
     :   -6,3,0,-6,3,0,-4,8,3,0,12,3,6,9,0,-6,3,0,6,21,12,12,12,-27,
     :   12,-6,27,6,-15,-12,-6,3,0,0,6,-3,0,0,-3,6,-3,0,0,-6,3,12,6,3,
     :   -4,2,0,-4,2,0,-4,2,0,-4,2,0,14,11,-9,14,-7,-9,-4,-2,0,-4,2,0,
     :   -2,7,3,2,11,3,8,8,0,0,0,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,-2,1,0,
     :   22,13,0/
*
*     ... coefficients of g(l+1) integrals
*
      data gint2/2,5,-3,2,-7,-15,-10,-25,-15,2,5,9,2,17,21,14,35,21,
     :   4,-2,-6,-8,-20,-12,-8,-20,-12,4,16,12,10,25,15,10,25,15,
     :   4,10,0,4,4,-12,2,-7,-21,-2,-17,-21,-8,-20,-12,0,0,0,
     :   -6,-15,-9,10,25,15,6,3,-3,0,-12,-9,-4,-16,-9,-6,-15,-9,
     :   -6,-15,-9,6,27,9,12,30,-9,12,12,-27,6,-9,-27,-6,-15,-9,
     :   0,0,-3,0,-6,-9,-6,-15,-9,12,18,9,0,6,9,6,15,9,
     :   -4,-10,-6,-4,-10,-6,-4,-10,-6,14,35,12,14,17,-6,-4,-10,-6,
     :   8,8,0,2,-7,-6,-2,-11,-6,-4,-10,-6,-4,-10,-6,0,0,0,
     :   -2,-5,-3,-2,-5,-3,-2,-5,-3,22,31,9,-2,-5,-3,-2,-5,-3/
*
*     ... encoded term value -- s = ltab/10
*                         lterm = l + (ltab mod 10 - 5)
*         example:  ltab = 36 with l = 2  is 3f
*
      data ltab/36,35,34,16,15,14,46,45,44,26,25,24,27,26,25,24,
     :   23,25,55,35,37,36,35,34,33,17,16,15,14,13,36,35,34,16,15,14,
     :   46,45,44,26,25,24,27,26,25,24,23,25,36,35,34,16,15,14/
      data sumtab/1,4,7,10,11/
      data partab/2,3,1,1,4,2,2,3,1,1,2/
      data ptrtab/6,12,17,18,20,30,36,42,47,48,54/
      data plval/1,1,2,0,0,2,1,1,2,0,1/
      data parch/'p','p','d','s','s','d','p','p','d','s','p'/
*
      ip = 1
 1    if (term(ip:ip) .eq. ' ') then
         ip = ip+1
         go to 1
      end if
      sl = term(ip:ip+1)
      senor = ' '
      if (ip.le.4) senor = term(ip+2:ip+2)
*
*   ---  convert lowercase l symbol to uppercase
*
      if (sl(2:2).gt.'a' .and. sl(2:2).lt.'z')
     :    sl(2:2) = char(ichar(sl(2:2)) + ichar('a') - ichar('a'))
*
*  ---  determine if fk or gk data needs to be input
*
      il = 0
      is = 0
      j = 1
      do 2 i = nclosd+1, nwf
         if (sum(i) .ne. 4*l(i)+2 .and. sum(i) .ne. 0.d0) then
            if (j.gt.2) then
               if (sl .ne. 'av' .and. sl .ne. 'av') then
                  done=.false.
               else
                  done = .true.
               end if
               return
            endif
            nos(j) = i
            j = j + 1
            if (l(i) .eq. 0 .and. is .eq. 0) then
               is = is + 1
               iis = i
            else
               il = il + 1
               iil = i
            end if
         end if
 2    continue
      if (sl .ne. 'av' .and. sl .ne. 'av' .and. is+il.ne.0) then
         done = .false.
         c = 0.d0
         if (is+il .le. 2 .and. il .le. 1) then
            if (is .eq. 0 .and. il .eq. 1) then
3              call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	       if (nsl .gt. 1) then
		  write(iscw,*)' ambiguous term: enter seniority'
		  read (5,'(a1)') senor
		  go to 3
	       end if
               call dev(iil,l(iil),sum(iil),ip,done)
            else if (is .eq. 1 .and. il .eq. 1) then
               slm = sl
	       slp = sl
               slm(1:1) = char(ichar(slm(1:1)) - 1)
               slp(1:1) = char(ichar(slp(1:1)) + 1)
	       call looktm(l(iil),slm,senor,sum(iil),ipm,nslm)
	       call looktm(l(iil),slp,senor,sum(iil),ipp,nslp)
	       if (nslm+nslp .eq. 0) then
	          done = .false.
	          return
	       else if (nslm .eq. 1 .and. nslp .eq. 0) then
		  sl = slm
		  ip = ipm
	       else if (nslm .eq. 0 .and. nslp .eq. 1) then
		  sl = slp
		  ip = ipp
	       else if (nslm .eq. 1 .and. nslp .eq. 1) then
4                 write(iscw,'(a,a3,a,a3)')
     :                 ' ambiguous l**n term: enter',slm, ' or ',slp 
	          read(5,'(a2)') sl
		  if (sl .eq. slm) then
		     ip = ipm
		  else if (sl .eq. slp) then
		     ip = ipp
		  else
	             write(iscw,*) ' term not allowed: re-enter'
		     go to 4
		  end if
	       else 
5                 write(iscw,'(a,a)') ' ambiguous l**n parent term:',
     :                   'enter term and seniority'
	          read(5,'(a2,a1)') sl, senor
	          if (senor .eq. ' ') then
			write(iscw,*) 'seniority is needed'
			go to 5
	          end if
	          call looktm(l(iil),sl,senor,sum(iil),ip,nsl)
	          if (nsl .ne. 1) then
		     write(iscw,'(a,a3,a,a3,a)') ' allowed terms are ',
     :                    slm, ' or ', slp,' plus seniority'
		     go to 5
	          end if
	       end if
               call dev(iil,l(iil),sum(iil),ip,done)
               if (done ) then
                  sp = ichar(sl(1:1)) - ichar('0')
                  csp = (sp - 1)/2.
		  if (sl .eq. slm) then
                     c = -csp/(2*l(iil)+1)
                  else
                     c = (csp + 1)/(2*l(iil)+1)
                  end if
                  call add(c,l(iil),iis,iil,.false.)
                  call add(c,l(iil),iil,iis,.false.)
               end if
            else if (is .eq. 1 .and. il .eq. 0) then
               done = .true.
            end if
         else
            if (((l(nos(1)).eq.1).and.(sum(nos(2)).eq.1.d0)).or.
     :          ((l(nos(2)).eq.1).and.(sum(nos(1)).eq.1.d0))) then
               if (l(nos(1)).eq.1.and.sum(nos(2)).eq.1.d0) then
                  isump=sum(nos(1))
                  np = nos(1)
                  nl = nos(2)
               else
                  isump = sum(nos(2))
                  np = nos(2)
                  nl = nos(1)
               endif
               sp=ichar(sl(1:1))-ichar('0')
               lp=lval(sl(2:2))
               ps1=sp+1
               ps2=sp-1
               if (isump.eq.1) then
                  iptr1=1
               else
                  iptr1=sumtab(isump-1)+1
               end if
               iptr2=sumtab(isump)
               nomach=0
               call lookup(partab,iptr1,iptr2,ind,nomach,ps1)
               call lookup(partab,iptr1,iptr2,ind,nomach,ps2)
               psl(1:1)=char(partab(ind)+ ichar('0'))
               psl(2:2)=parch(ind)
               if (nomach.gt.1) then
                  write(iscw,*)' ambiguous parent case'
 10               write(iscw,*)' enter the sl term for p(n) subshell'
                  read (5,'(a)')psl
                  if (psl(2:2).gt.'a'.and.psl(2:2).lt.'z')
     :                   psl(2:2)=char(ichar(psl(2:2))+ichar('a')
     :                                                -ichar('a'))
                  ps1=ichar(psl(1:1))-ichar('0')
                  ps2=lval(psl(2:2))
                  call lookup(plval,iptr1,iptr2,ind,nomach,ps2)
                  if ((nomach.ne.1).and.(partab(ind).ne.ps1))
     :                    go to 10
               end if
               if (isump.eq.1) then
                  iptr1=1
               else
                  iptr1=ptrtab(ind-1)+1
               end if
               iptr2 = ptrtab(ind)
               lv=l(nl)
               pacval=sp*10+lp-lv+5
               nomach=0
               call lookup(ltab,iptr1,iptr2,ind,nomach,pacval)
               if (nomach.ne.1) then
                  done=.false.
                  return
               endif
               val1=((fint(1,ind)*lv+fint(2,ind))*lv+fint(3,ind))
     :                  /(5.d0*(2*lv-1)*(2*lv+3))
               val2=((gint1(1,ind)*lv+gint1(2,ind))*lv+gint1(3,ind))
     :                  /(2.d0*(2*lv+1)*(2*lv-1)**2)
               val3=((gint2(1,ind)*lv+gint2(2,ind))*lv+gint2(3,ind))
     :                  /(2.d0*(2*lv+1)*(2*lv+3)**2)
*
*     ...  add contributions from between p-subshell and l-electron
*
               call add(val1,2,np,nl,.true.)
               call add(val1,2,nl,np,.true.)
               call add(val2,lv-1,np,nl,.false.)
               call add(val2,lv-1,nl,np,.false.)
               call add(val3,lv+1,np,nl,.false.)
               call add(val3,lv+1,nl,np,.false.)
*
*     ... add deviations for p-subshell
*
               call looktm(1,psl,' ',sum(np),ip,nsl)
               call dev(np,1,sum(np),ip,done)
            else
               done = .false.
            end if
         end if
      else
         done = .true.
      end if
      return
      end
