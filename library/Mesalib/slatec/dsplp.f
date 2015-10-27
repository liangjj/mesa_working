*deck dsplp
      subroutine dsplp (dusrmt, mrelas, nvars, costs, prgopt, dattrv,
     +   bl, bu, ind, info, primal, duals, ibasis, work, lw, iwork, liw)
c***begin prologue  dsplp
c***purpose  solve linear programming problems involving at
c            most a few thousand constraints and variables.
c            takes advantage of sparsity in the constraint matrix.
c***library   slatec
c***category  g2a2
c***type      double precision (splp-s, dsplp-d)
c***keywords  linear constraints, linear optimization,
c             linear programming, lp, sparse constraints
c***author  hanson, r. j., (snla)
c           hiebert, k. l., (snla)
c***description
c
c     these are the short usage instructions; for details about
c     other features, options and methods for defining the matrix
c     a, see the extended usage instructions which are contained in
c     the long description section below.
c
c   |------------|
c   |introduction|
c   |------------|
c     the subprogram dsplp( ) solves a linear optimization problem.
c     the problem statement is as follows
c
c                         minimize (transpose of costs)*x
c                         subject to a*x=w.
c
c     the entries of the unknowns x and w may have simple lower or
c     upper bounds (or both), or be free to take on any value.  by
c     setting the bounds for x and w, the user is imposing the con-
c     straints of the problem.  the matrix a has mrelas rows and
c     nvars columns.  the vectors costs, x, and w respectively
c     have nvars, nvars, and mrelas number of entries.
c
c     the input for the problem includes the problem dimensions,
c     mrelas and nvars, the array costs(*), data for the matrix
c     a, and the bound information for the unknowns x and w, bl(*),
c     bu(*), and ind(*).  only the nonzero entries of the matrix a
c     are passed to dsplp( ).
c
c     the output from the problem (when output flag info=1) includes
c     optimal values for x and w in primal(*), optimal values for
c     dual variables of the equations a*x=w and the simple bounds
c     on x in  duals(*), and the indices of the basic columns,
c     ibasis(*).
c
c  |------------------------------|
c  |fortran declarations required:|
c  |------------------------------|
c
c     dimension costs(nvars),prgopt(*),dattrv(*),
c    *bl(nvars+mrelas),bu(nvars+mrelas),ind(nvars+mrelas),
c    *primal(nvars+mrelas),duals(mrelas+nvars),ibasis(nvars+mrelas),
c    *work(lw),iwork(liw)
c
c     external dusrmt
c
c     the dimensions of prgopt(*) and dattrv(*) must be at least 1.
c     the exact lengths will be determined by user-required options and
c     data transferred to the subprogram dusrmt( ).
c
c     the values of lw and liw, the lengths of the arrays work(*)
c     and iwork(*), must satisfy the inequalities
c
c               lw .ge. 4*nvars+ 8*mrelas+lamat+  lbm
c               liw.ge.   nvars+11*mrelas+lamat+2*lbm
c
c     it is an error if they do not both satisfy these inequalities.
c     (the subprogram will inform the user of the required lengths
c     if either lw or liw is wrong.)  the values of lamat and lbm
c     nominally are
c
c               lamat=4*nvars+7
c     and       lbm  =8*mrelas
c
c     lamat determines the length of the sparse matrix storage area.
c     the value of lbm determines the amount of storage available
c     to decompose and update the active basis matrix.
c
c  |------|
c  |input:|
c  |------|
c
c     mrelas,nvars
c     ------------
c     these parameters are respectively the number of constraints (the
c     linear relations a*x=w that the unknowns x and w are to satisfy)
c     and the number of entries in the vector x.  both must be .ge. 1.
c     other values are errors.
c
c     costs(*)
c     --------
c     the nvars entries of this array are the coefficients of the
c     linear objective function.  the value costs(j) is the
c     multiplier for variable j of the unknown vector x.  each
c     entry of this array must be defined.
c
c     dusrmt
c     ------
c     this is the name of a specific subprogram in the dsplp( ) package
c     used to define the matrix a.  in this usage mode of dsplp( )
c     the user places the nonzero entries of a in the
c     array dattrv(*) as given in the description of that parameter.
c     the name dusrmt must appear in a fortran external statement.
c
c     dattrv(*)
c     ---------
c     the array dattrv(*) contains data for the matrix a as follows:
c     each column (numbered j) requires (floating point) data con-
c     sisting of the value (-j) followed by pairs of values.  each pair
c     consists of the row index immediately followed by the value
c     of the matrix at that entry.  a value of j=0 signals that there
c     are no more columns.  the required length of
c     dattrv(*) is 2*no. of nonzeros + nvars + 1.
c
c     bl(*),bu(*),ind(*)
c     ------------------
c     the values of ind(*) are input parameters that define
c     the form of the bounds for the unknowns x and w.  the values for
c     the bounds are found in the arrays bl(*) and bu(*) as follows.
c
c     for values of j between 1 and nvars,
c          if ind(j)=1, then x(j) .ge. bl(j); bu(j) is not used.
c          if ind(j)=2, then x(j) .le. bu(j); bl(j) is not used.
c          if ind(j)=3, then bl(j) .le. x(j) .le. bu(j),(bl(j)=bu(j) ok)
c          if ind(j)=4, then x(j) is free to have any value,
c          and bl(j), bu(j) are not used.
c
c     for values of i between nvars+1 and nvars+mrelas,
c          if ind(i)=1, then w(i-nvars) .ge. bl(i); bu(i) is not used.
c          if ind(i)=2, then w(i-nvars) .le. bu(i); bl(i) is not used.
c          if ind(i)=3, then bl(i) .le. w(i-nvars) .le. bu(i),
c          (bl(i)=bu(i) is ok).
c          if ind(i)=4, then w(i-nvars) is free to have any value,
c          and bl(i), bu(i) are not used.
c
c     a value of ind(*) not equal to 1,2,3 or 4 is an error.  when
c     ind(i)=3, bl(i) must be .le. bu(i).  the condition bl(i).gt.
c     bu(i) indicates infeasibility and is an error.
c
c     prgopt(*)
c     ---------
c     this array is used to redefine various parameters within dsplp( ).
c     frequently, perhaps most of the time, a user will be satisfied
c     and obtain the solutions with no changes to any of these
c     parameters.  to try this, simply set prgopt(1)=1.d0.
c
c     for users with more sophisticated needs, dsplp( ) provides several
c     options that may be used to take advantage of more detailed
c     knowledge of the problem or satisfy other utilitarian needs.
c     the complete description of how to use this option array to
c     utilize additional subprogram features is found under the
c     heading  of dsplp( ) subprogram options in the extended
c     usage instructions.
c
c     briefly, the user should note the following value of the parameter
c     key and the corresponding task or feature desired before turning
c     to that document.
c
c     value     brief statement of purpose for option
c     of key
c     ------    -------------------------------------
c     50        change from a minimization problem to a
c               maximization problem.
c     51        change the amount of printed output.
c               normally, no printed output is obtained.
c     52        redefine the line length and precision used
c               for the printed output.
c     53        redefine the values of lamat and lbm that
c               were discussed above under the heading
c               fortran declarations required.
c     54        redefine the unit number where pages of the sparse
c               data matrix a are stored.  normally, the unit
c               number is 1.
c     55        a computation, partially completed, is
c               being continued.  read the up-to-date
c               partial results from unit number 2.
c     56        redefine the unit number where the partial results
c               are stored.  normally, the unit number is 2.
c     57        save partial results on unit 2 either after
c               maximum iterations or at the optimum.
c     58        redefine the value for the maximum number of
c               iterations.  normally, the maximum number of
c               iterations is 3*(nvars+mrelas).
c     59        provide dsplp( ) with a starting (feasible)
c               nonsingular basis.  normally, dsplp( ) starts
c               with the identity matrix columns corresponding
c               to the vector w.
c     60        the user has provided scale factors for the
c               columns of a.  normally, dsplp( ) computes scale
c               factors that are the reciprocals of the max. norm
c               of each column.
c     61        the user has provided a scale factor
c               for the vector costs.  normally, dsplp( ) computes
c               a scale factor equal to the reciprocal of the
c               max. norm of the vector costs after the column
c               scaling for the data matrix has been applied.
c     62        size parameters, namely the smallest and
c               largest magnitudes of nonzero entries in
c               the matrix a, are provided.  values noted
c               outside this range are to be considered errors.
c     63        redefine the tolerance required in
c               evaluating residuals for feasibility.
c               normally, this value is set to relpr,
c               where relpr = relative precision of the arithmetic.
c     64        change the criterion for bringing new variables
c               into the basis from the steepest edge (best
c               local move) to the minimum reduced cost.
c     65        redefine the value for the number of iterations
c               between recalculating the error in the primal
c               solution.  normally, this value is equal to ten.
c     66        perform "partial pricing" on variable selection.
c               redefine the value for the number of negative
c               reduced costs to compute (at most) when finding
c               a variable to enter the basis.  normally this
c               value is set to nvars.  this implies that no
c               "partial pricing" is used.
c     67        adjust the tuning factor (normally one) to apply
c               to the primal and dual error estimates.
c     68        pass  information to the  subprogram  dfulmt(),
c               provided with the dsplp() package, so that a fortran
c               two-dimensional array can be used as the argument
c               dattrv(*).
c     69        pass an absolute tolerance to use for the feasibility
c               test when the usual relative error test indicates
c               infeasibility.  the nominal value of this tolerance,
c               tolabs, is zero.
c
c
c  |---------------|
c  |working arrays:|
c  |---------------|
c
c     work(*),lw,
c     iwork(*),liw
c     ------------
c     the arrays work(*) and iwork(*) are respectively floating point
c     and type integer working arrays for dsplp( ) and its
c     subprograms.  the lengths of these arrays are respectively
c     lw and liw.  these parameters must satisfy the inequalities
c     noted above under the heading "fortran declarations required:"
c     it is an error if either value is too small.
c
c  |----------------------------|
c  |input/output files required:|
c  |----------------------------|
c
c     fortran unit 1 is used by dsplp( ) to store the sparse matrix a
c     out of high-speed memory.  a crude
c     upper bound for the amount of information written on unit 1
c     is 6*nz, where nz is the number of nonzero entries in a.
c
c  |-------|
c  |output:|
c  |-------|
c
c     info,primal(*),duals(*)
c     -----------------------
c     the integer flag info indicates why dsplp( ) has returned to the
c     user.  if info=1 the solution has been computed.  in this case
c     x(j)=primal(j) and w(i)=primal(i+nvars).  the dual variables
c     for the equations a*x=w are in the array duals(i)=dual for
c     equation number i.  the dual value for the component x(j) that
c     has an upper or lower bound (or both) is returned in
c     duals(j+mrelas).  the only other values for info are .lt. 0.
c     the meaning of these values can be found by reading
c     the diagnostic message in the output file, or by looking for
c     error number = (-info) in the extended usage instructions
c     under the heading:
c
c          list of dsplp( ) error and diagnostic messages.
c
c     bl(*),bu(*),ind(*)
c     ------------------
c     these arrays are output parameters only under the (unusual)
c     circumstances where the stated problem is infeasible, has an
c     unbounded optimum value, or both.  these respective conditions
c     correspond to info=-1,-2 or -3.    see the extended
c     usage instructions for further details.
c
c     ibasis(i),i=1,...,mrelas
c     ------------------------
c     this array contains the indices of the variables that are
c     in the active basis set at the solution (info=1).  a value
c     of ibasis(i) between 1 and nvars corresponds to the variable
c     x(ibasis(i)).  a value of ibasis(i) between nvars+1 and nvars+
c     mrelas corresponds to the variable w(ibasis(i)-nvars).
c
c *long description:
c
c     subroutine dsplp(dusrmt,mrelas,nvars,costs,prgopt,dattrv,
c    *           bl,bu,ind,info,primal,duals,ibasis,work,lw,iwork,liw)
c
c     |------------|
c     |introduction|
c     |------------|
c     the subprogram dsplp( ) solves a linear optimization problem.
c     the problem statement is as follows
c
c                         minimize (transpose of costs)*x
c                         subject to a*x=w.
c
c     the entries of the unknowns x and w may have simple lower or
c     upper bounds (or both), or be free to take on any value.  by
c     setting the bounds for x and w, the user is imposing the con-
c     straints of the problem.
c
c     (the problem may also be stated as a maximization
c     problem.  this is done by means of input in the option array
c     prgopt(*).)  the matrix a has mrelas rows and nvars columns.  the
c     vectors costs, x, and w respectively have nvars, nvars, and
c     mrelas number of entries.
c
c     the input for the problem includes the problem dimensions,
c     mrelas and nvars, the array costs(*), data for the matrix
c     a, and the bound information for the unknowns x and w, bl(*),
c     bu(*), and ind(*).
c
c     the output from the problem (when output flag info=1) includes
c     optimal values for x and w in primal(*), optimal values for
c     dual variables of the equations a*x=w and the simple bounds
c     on x in  duals(*), and the indices of the basic columns in
c     ibasis(*).
c
c  |------------------------------|
c  |fortran declarations required:|
c  |------------------------------|
c
c     dimension costs(nvars),prgopt(*),dattrv(*),
c    *bl(nvars+mrelas),bu(nvars+mrelas),ind(nvars+mrelas),
c    *primal(nvars+mrelas),duals(mrelas+nvars),ibasis(nvars+mrelas),
c    *work(lw),iwork(liw)
c
c     external dusrmt (or 'name', if user provides the subprogram)
c
c     the dimensions of prgopt(*) and dattrv(*) must be at least 1.
c     the exact lengths will be determined by user-required options and
c     data transferred to the subprogram dusrmt( ) ( or 'name').
c
c     the values of lw and liw, the lengths of the arrays work(*)
c     and iwork(*), must satisfy the inequalities
c
c               lw .ge. 4*nvars+ 8*mrelas+lamat+  lbm
c               liw.ge.   nvars+11*mrelas+lamat+2*lbm
c
c     it is an error if they do not both satisfy these inequalities.
c     (the subprogram will inform the user of the required lengths
c     if either lw or liw is wrong.)  the values of lamat and lbm
c     nominally are
c
c               lamat=4*nvars+7
c     and       lbm  =8*mrelas
c
c     these values will be as shown unless the user changes them by
c     means of input in the option array prgopt(*).  the value of lamat
c     determines the length of the sparse matrix "staging" area.
c     for reasons of efficiency the user may want to increase the value
c     of lamat.  the value of lbm determines the amount of storage
c     available to decompose and update the active basis matrix.
c     due to exhausting the working space because of fill-in,
c     it may be necessary for the user to increase the value of lbm.
c     (if this situation occurs an informative diagnostic is printed
c     and a value of info=-28 is obtained as an output parameter.)
c
c  |------|
c  |input:|
c  |------|
c
c     mrelas,nvars
c     ------------
c     these parameters are respectively the number of constraints (the
c     linear relations a*x=w that the unknowns x and w are to satisfy)
c     and the number of entries in the vector x.  both must be .ge. 1.
c     other values are errors.
c
c     costs(*)
c     --------
c     the nvars entries of this array are the coefficients of the
c     linear objective function.  the value costs(j) is the
c     multiplier for variable j of the unknown vector x.  each
c     entry of this array must be defined.  this array can be changed
c     by the user between restarts.  see options with key=55,57 for
c     details of checkpointing and restarting.
c
c     dusrmt
c     ------
c     this is the name of a specific subprogram in the dsplp( ) package
c     that is used to define the matrix entries when this data is passed
c     to dsplp( ) as a linear array.  in this usage mode of dsplp( )
c     the user gives information about the nonzero entries of a
c     in dattrv(*) as given under the description of that parameter.
c     the name dusrmt must appear in a fortran external statement.
c     users who are passing the matrix data with dusrmt( ) can skip
c     directly to the description of the input parameter dattrv(*).
c     also see option 68 for passing the constraint matrix data using
c     a standard fortran two-dimensional array.
c
c     if the user chooses to provide a subprogram 'name'( ) to
c     define the matrix a, then dattrv(*) may be used to pass floating
c     point data from the user's program unit to the subprogram
c     'name'( ). the content of dattrv(*) is not changed in any way.
c
c     the subprogram 'name'( ) can be of the user's choice
c     but it must meet fortran standards and it must appear in a
c     fortran external statement.  the first statement of the subprogram
c     has the form
c
c          subroutine 'name'(i,j,aij, indcat, prgopt, dattrv, iflag)
c
c     the variables i,j, indcat, iflag(10) are type integer,
c          while  aij, prgopt(*),dattrv(*) are type real.
c
c     the user interacts with the contents of iflag(*) to
c     direct the appropriate action.  the algorithmic steps are
c     as follows.
c
c          test iflag(1).
c
c             if(iflag(1).eq.1) then
c
c               initialize the necessary pointers and data
c               for defining the matrix a.  the contents
c               of iflag(k), k=2,...,10, may be used for
c               storage of the pointers.  this array remains intact
c               between calls to 'name'( ) by dsplp( ).
c               return
c
c             end if
c
c             if(iflag(1).eq.2) then
c
c               define one set of values for i,j,aij, and indcat.
c               each nonzero entry of a must be defined this way.
c               these values can be defined in any convenient order.
c               (it is most efficient to define the data by
c               columns in the order 1,...,nvars; within each
c               column define the entries in the order 1,...,mrelas.)
c               if this is the last matrix value to be
c               defined or updated, then set iflag(1)=3.
c               (when i and j are positive and respectively no larger
c               than mrelas and nvars, the value of aij is used to
c               define (or update) row i and column j of a.)
c               return
c
c             end if
c
c               end
c
c     remarks:  the values of i and j are the row and column
c               indices for the nonzero entries of the matrix a.
c               the value of this entry is aij.
c               set indcat=0 if this value defines that entry.
c               set indcat=1 if this entry is to be updated,
c                            new entry=old entry+aij.
c               a value of i not between 1 and mrelas, a value of j
c               not between 1 and nvars, or a value of indcat
c               not equal to 0 or 1 are each errors.
c
c               the contents of iflag(k), k=2,...,10, can be used to
c               remember the status (of the process of defining the
c               matrix entries) between calls to 'name'( ) by dsplp( ).
c               on entry to 'name'( ), only the values 1 or 2 will be
c               in iflag(1).  more than 2*nvars*mrelas definitions of
c               the matrix elements is considered an error because
c               it suggests an infinite loop in the user-written
c               subprogram 'name'( ).  any matrix element not
c               provided by 'name'( ) is defined to be zero.
c
c               the real arrays prgopt(*) and dattrv(*) are passed as
c               arguments directly from dsplp( ) to 'name'( ).
c               the array prgopt(*) contains any user-defined program
c               options.  in this usage mode the array dattrv(*) may
c               now contain any (type real) data that the user needs
c               to define the matrix a.  both arrays prgopt(*) and
c               dattrv(*) remain intact between calls to 'name'( )
c               by dsplp( ).
c     here is a subprogram that communicates the matrix values for a,
c     as represented in dattrv(*), to dsplp( ).  this subprogram,
c     called dusrmt( ), is included as part of the dsplp( ) package.
c     this subprogram 'decodes' the array dattrv(*) and defines the
c     nonzero entries of the matrix  a for dsplp( ) to store.  this
c     listing is presented here as a guide and example
c     for the users who find it necessary to write their own subroutine
c     for this purpose.  the contents of dattrv(*) are given below in
c     the description of that parameter.
c
c     subroutine dusrmt(i,j,aij, indcat,prgopt,dattrv,iflag)
c     dimension prgopt(*),dattrv(*),iflag(10)
c
c     if(iflag(1).eq.1) then
c
c     this is the initialization step.  the values of iflag(k),k=2,3,4,
c     are respectively the column index, the row index (or the next col.
c     index), and the pointer to the matrix entry's value within
c     dattrv(*).  also check (dattrv(1)=0.) signifying no data.
c          if(dattrv(1).eq.0.) then
c          i = 0
c          j = 0
c          iflag(1) = 3
c          else
c          iflag(2)=-dattrv(1)
c          iflag(3)= dattrv(2)
c          iflag(4)= 3
c          end if
c
c          return
c     else
c          j=iflag(2)
c          i=iflag(3)
c          l=iflag(4)
c          if(i.eq.0) then
c
c     signal that all of the nonzero entries have been defined.
c               iflag(1)=3
c               return
c          else if(i.lt.0) then
c
c     signal that a switch is made to a new column.
c               j=-i
c               i=dattrv(l)
c               l=l+1
c          end if
c
c          aij=dattrv(l)
c
c     update the indices and pointers for the next entry.
c          iflag(2)=j
c          iflag(3)=dattrv(l+1)
c          iflag(4)=l+2
c
c     indcat=0 denotes that entries of the matrix are assigned the
c     values from dattrv(*).  no accumulation is performed.
c          indcat=0
c          return
c     end if
c     end
c
c     dattrv(*)
c     ---------
c     if the user chooses to use the provided subprogram dusrmt( ) then
c     the array dattrv(*) contains data for the matrix a as follows:
c     each column (numbered j) requires (floating point) data con-
c     sisting of the value (-j) followed by pairs of values.  each pair
c     consists of the row index immediately followed by the value
c     of the matrix at that entry.  a value of j=0 signals that there
c     are no more columns.  (see "example of dsplp( ) usage," below.)
c     the dimension of dattrv(*) must be 2*no. of nonzeros
c     + nvars + 1 in this usage.  no checking of the array
c     length is done by the subprogram package.
c
c     if the save/restore feature is in use (see options with
c     key=55,57 for details of checkpointing and restarting)
c     dusrmt( ) can be used to redefine entries of the matrix.
c     the matrix entries are redefined or overwritten.  no accum-
c     ulation is performed.
c     any other nonzero entry of a, defined in a previous call to
c     dsplp( ), remain intact.
c
c     bl(*),bu(*),ind(*)
c     ------------------
c     the values of ind(*) are input parameters that define
c     the form of the bounds for the unknowns x and w.  the values for
c     the bounds are found in the arrays bl(*) and bu(*) as follows.
c
c     for values of j between 1 and nvars,
c          if ind(j)=1, then x(j) .ge. bl(j); bu(j) is not used.
c          if ind(j)=2, then x(j) .le. bu(j); bl(j) is not used.
c          if ind(j)=3, then bl(j) .le. x(j) .le. bu(j),(bl(j)=bu(j) ok)
c          if ind(j)=4, then x(j) is free to have any value,
c          and bl(j), bu(j) are not used.
c
c     for values of i between nvars+1 and nvars+mrelas,
c          if ind(i)=1, then w(i-nvars) .ge. bl(i); bu(i) is not used.
c          if ind(i)=2, then w(i-nvars) .le. bu(i); bl(i) is not used.
c          if ind(i)=3, then bl(i) .le. w(i-nvars) .le. bu(i),
c          (bl(i)=bu(i) is ok).
c          if ind(i)=4, then w(i-nvars) is free to have any value,
c          and bl(i), bu(i) are not used.
c
c     a value of ind(*) not equal to 1,2,3 or 4 is an error.  when
c     ind(i)=3, bl(i) must be .le. bu(i).  the condition bl(i).gt.
c     bu(i) indicates infeasibility and is an error.  these
c     arrays can be changed by the user between restarts.  see
c     options with key=55,57 for details of checkpointing and
c     restarting.
c
c     prgopt(*)
c     ---------
c     this array is used to redefine various parameters within dsplp( ).
c     frequently, perhaps most of the time, a user will be satisfied
c     and obtain the solutions with no changes to any of these
c     parameters.  to try this, simply set prgopt(1)=1.d0.
c
c     for users with more sophisticated needs, dsplp( ) provides several
c     options that may be used to take advantage of more detailed
c     knowledge of the problem or satisfy other utilitarian needs.
c     the complete description of how to use this option array to
c     utilize additional subprogram features is found under the
c     heading "usage of dsplp( ) subprogram options."
c
c     briefly, the user should note the following value of the parameter
c     key and the corresponding task or feature desired before turning
c     to that section.
c
c     value     brief statement of purpose for option
c     of key
c     ------    -------------------------------------
c     50        change from a minimization problem to a
c               maximization problem.
c     51        change the amount of printed output.
c               normally, no printed output is obtained.
c     52        redefine the line length and precision used
c               for the printed output.
c     53        redefine the values of lamat and lbm that
c               were discussed above under the heading
c               fortran declarations required.
c     54        redefine the unit number where pages of the sparse
c               data matrix a are stored.  normally, the unit
c               number is 1.
c     55        a computation, partially completed, is
c               being continued.  read the up-to-date
c               partial results from unit number 2.
c     56        redefine the unit number where the partial results
c               are stored.  normally, the unit number is 2.
c     57        save partial results on unit 2 either after
c               maximum iterations or at the optimum.
c     58        redefine the value for the maximum number of
c               iterations.  normally, the maximum number of
c               iterations is 3*(nvars+mrelas).
c     59        provide dsplp( ) with a starting (feasible)
c               nonsingular basis.  normally, dsplp( ) starts
c               with the identity matrix columns corresponding
c               to the vector w.
c     60        the user has provided scale factors for the
c               columns of a.  normally, dsplp( ) computes scale
c               factors that are the reciprocals of the max. norm
c               of each column.
c     61        the user has provided a scale factor
c               for the vector costs.  normally, dsplp( ) computes
c               a scale factor equal to the reciprocal of the
c               max. norm of the vector costs after the column
c               scaling for the data matrix has been applied.
c     62        size parameters, namely the smallest and
c               largest magnitudes of nonzero entries in
c               the matrix a, are provided.  values noted
c               outside this range are to be considered errors.
c     63        redefine the tolerance required in
c               evaluating residuals for feasibility.
c               normally, this value is set to the value relpr,
c               where relpr = relative precision of the arithmetic.
c     64        change the criterion for bringing new variables
c               into the basis from the steepest edge (best
c               local move) to the minimum reduced cost.
c     65        redefine the value for the number of iterations
c               between recalculating the error in the primal
c               solution.  normally, this value is equal to ten.
c     66        perform "partial pricing" on variable selection.
c               redefine the value for the number of negative
c               reduced costs to compute (at most) when finding
c               a variable to enter the basis.  normally this
c               value is set to nvars.  this implies that no
c               "partial pricing" is used.
c     67        adjust the tuning factor (normally one) to apply
c               to the primal and dual error estimates.
c     68        pass  information to the  subprogram  dfulmt(),
c               provided with the dsplp() package, so that a fortran
c               two-dimensional array can be used as the argument
c               dattrv(*).
c     69        pass an absolute tolerance to use for the feasibility
c               test when the usual relative error test indicates
c               infeasibility.  the nominal value of this tolerance,
c               tolabs, is zero.
c
c
c  |---------------|
c  |working arrays:|
c  |---------------|
c
c     work(*),lw,
c     iwork(*),liw
c     ------------
c     the arrays work(*) and iwork(*) are respectively floating point
c     and type integer working arrays for dsplp( ) and its
c     subprograms.  the lengths of these arrays are respectively
c     lw and liw.  these parameters must satisfy the inequalities
c     noted above under the heading "fortran declarations required."
c     it is an error if either value is too small.
c
c  |----------------------------|
c  |input/output files required:|
c  |----------------------------|
c
c     fortran unit 1 is used by dsplp( ) to store the sparse matrix a
c     out of high-speed memory.  this direct access file is opened
c     within the package under the following two conditions.
c     1. when the save/restore feature is used.  2. when the
c     constraint matrix is so large that storage out of high-speed
c     memory is required.  the user may need to close unit 1
c     (with deletion from the job step) in the main program unit
c     when several calls are made to dsplp( ).  a crude
c     upper bound for the amount of information written on unit 1
c     is 6*nz, where nz is the number of nonzero entries in a.
c     the unit number may be redefined to any other positive value
c     by means of input in the option array prgopt(*).
c
c     fortran unit 2 is used by dsplp( ) only when the save/restore
c     feature is desired.  normally this feature is not used.  it is
c     activated by means of input in the option array prgopt(*).
c     on some computer systems the user may need to open unit
c     2 before executing a call to dsplp( ).  this file is type
c     sequential and is unformatted.
c
c     fortran unit=i1mach(2) (check local setting) is used by dsplp( )
c     when the printed output feature (key=51) is used.  normally
c     this feature is not used.  it is activated by input in the
c     options array prgopt(*).  for many computer systems i1mach(2)=6.
c
c  |-------|
c  |output:|
c  |-------|
c
c     info,primal(*),duals(*)
c     -----------------------
c     the integer flag info indicates why dsplp( ) has returned to the
c     user.  if info=1 the solution has been computed.  in this case
c     x(j)=primal(j) and w(i)=primal(i+nvars).  the dual variables
c     for the equations a*x=w are in the array duals(i)=dual for
c     equation number i.  the dual value for the component x(j) that
c     has an upper or lower bound (or both) is returned in
c     duals(j+mrelas).  the only other values for info are .lt. 0.
c     the meaning of these values can be found by reading
c     the diagnostic message in the output file, or by looking for
c     error number = (-info) under the heading "list of dsplp( ) error
c     and diagnostic messages."
c     the diagnostic messages are printed using the error processing
c     subprogram xermsg( ) with error category level=1.
c     see the document "brief instr. for using the sandia math.
c     subroutine library," sand79-2382, nov., 1980, for further inform-
c     ation about resetting the usual response to a diagnostic message.
c
c     bl(*),bu(*),ind(*)
c     ------------------
c     these arrays are output parameters only under the (unusual)
c     circumstances where the stated problem is infeasible, has an
c     unbounded optimum value, or both.  these respective conditions
c     correspond to info=-1,-2 or -3.  for info=-1 or -3 certain comp-
c     onents of the vectors x or w will not satisfy the input bounds.
c     if component j of x or component i of w does not satisfy its input
c     bound because of infeasibility, then ind(j)=-4 or ind(i+nvars)=-4,
c     respectively.  for info=-2 or -3 certain
c     components of the vector x could not be used as basic variables
c     because the objective function would have become unbounded.
c     in particular if component j of x corresponds to such a variable,
c     then ind(j)=-3.  further, if the input value of ind(j)
c                      =1, then bu(j)=bl(j);
c                      =2, then bl(j)=bu(j);
c                      =4, then bl(j)=0.,bu(j)=0.
c
c     (the j-th variable in x has been restricted to an appropriate
c     feasible value.)
c     the negative output value for ind(*) allows the user to identify
c     those constraints that are not satisfied or those variables that
c     would cause unbounded values of the objective function.  note
c     that the absolute value of ind(*), together with bl(*) and bu(*),
c     are valid input to dsplp( ).  in the case of infeasibility the
c     sum of magnitudes of the infeasible values is minimized.  thus
c     one could reenter dsplp( ) with these components of x or w now
c     fixed at their present values.  this involves setting
c     the appropriate components of ind(*) = 3, and bl(*) = bu(*).
c
c     ibasis(i),i=1,...,mrelas
c     ------------------------
c     this array contains the indices of the variables that are
c     in the active basis set at the solution (info=1).  a value
c     of ibasis(i) between 1 and nvars corresponds to the variable
c     x(ibasis(i)).  a value of ibasis(i) between nvars+1 and nvars+
c     mrelas corresponds to the variable w(ibasis(i)-nvars).
c
c     computing with the matrix a after calling dsplp( )
c     --------------------------------------------------
c     following the return from dsplp( ), nonzero entries of the mrelas
c     by nvars matrix a are available for usage by the user.  the method
c     for obtaining the next nonzero in column j with a row index
c     strictly greater than i in value, is completed by executing
c
c         call dpnnzr(i,aij,iplace,work,iwork,j)
c
c     the value of i is also an output parameter.  if i.le.0 on output,
c     then there are no more nonzeroes in column j.  if i.gt.0, the
c     output value for component number i of column j is in aij.  the
c     parameters work(*) and iwork(*) are the same arguments as in the
c     call to dsplp( ).  the parameter iplace is a single integer
c     working variable.
c
c     the data structure used for storage of the matrix a within dsplp()
c     corresponds to sequential storage by columns as defined in
c     sand78-0785.  note that the names of the subprograms lnnzrs(),
c     lchngs(),linitm(),lloc(),lrwpge(), and lrwvir() have been
c     changed to dpnnzr(),dpchng(),pinitm(),iploc(),dprwpg(), and
c     dprwvr() respectively.  the error processing subprogram lerror()
c     is no longer used; xermsg() is used instead.
c
c  |--------------------------------|
c  |subprograms required by dsplp( )|
c  |--------------------------------|
c     called by dsplp() are dplpmn(),dplpup(),dpinit(),dpopt(),
c                          dplpdm(),dplpce(),dpincw(),dplpfl(),
c                          dplpfe(),dplpmu().
c
c     error processing subprograms xermsg(),i1mach(),d1mach()
c
c     sparse matrix subprograms dpnnzr(),dpchng(),dprwpg(),dprwvr(),
c                               pinitm(),iploc()
c
c     mass storage file subprograms sopenm(),sclosm(),dreadp(),dwritp()
c
c     basic linear algebra subprograms dcopy(),dasum(),ddot()
c
c     sparse matrix basis handling subprograms la05ad(),la05bd(),
c                                             la05cd(),la05ed(),mc20ad()
c
c     vector output subprograms dvout(),ivout()
c
c     machine-sensitive subprograms i1mach( ),d1mach( ),
c                                   sopenm(),sclosm(),dreadp(),dwritp().
c     common block used
c     -----------------
c     /la05dd/ small,lp,lenl,lenu,ncp,lrow,lcol
c     see the document aere-r8269 for further details.
c    |-------------------------|
c    |example of dsplp( ) usage|
c    |-------------------------|
c     program lpex
c     the optimization problem is to find x1, x2, x3 that
c     minimize x1 + x2 + x3, x1.ge.0, x2.ge.0, x3 unconstrained.
c
c     the unknowns x1,x2,x3 are to satisfy constraints
c
c        x1 -3*x2 +4*x3 = 5
c        x1 -2*x2     .le.3
c            2*x2 - x3.ge.4
c
c     we first define the dependent variables
c          w1=x1 -3*x2 +4*x3
c          w2=x1- 2*x2
c          w3=    2*x2 -x3
c
c     we now show how to use dsplp( ) to solve this linear optimization
c     problem.  each required step will be shown in this example.
c     dimension costs(03),prgopt(01),dattrv(18),bl(06),bu(06),ind(06),
c    *primal(06),duals(06),ibasis(06),work(079),iwork(103)
c
c     external dusrmt
c     mrelas=3
c     nvars=3
c
c     define the array costs(*) for the objective function.
c     costs(01)=1.
c     costs(02)=1.
c     costs(03)=1.
c
c     place the nonzero information about the matrix in dattrv(*).
c     define col. 1:
c     dattrv(01)=-1
c     dattrv(02)=1
c     dattrv(03)=1.
c     dattrv(04)=2
c     dattrv(05)=1.
c
c     define col. 2:
c     dattrv(06)=-2
c     dattrv(07)=1
c     dattrv(08)=-3.
c     dattrv(09)=2
c     dattrv(10)=-2.
c     dattrv(11)=3
c     dattrv(12)=2.
c
c     define col. 3:
c     dattrv(13)=-3
c     dattrv(14)=1
c     dattrv(15)=4.
c     dattrv(16)=3
c     dattrv(17)=-1.
c
c     dattrv(18)=0
c
c     constrain x1,x2 to be nonnegative. let x3 have no bounds.
c     bl(1)=0.
c     ind(1)=1
c     bl(2)=0.
c     ind(2)=1
c     ind(3)=4
c
c     constrain w1=5,w2.le.3, and w3.ge.4.
c     bl(4)=5.
c     bu(4)=5.
c     ind(4)=3
c     bu(5)=3.
c     ind(5)=2
c     bl(6)=4.
c     ind(6)=1
c
c     indicate that no modifications to options are in use.
c     prgopt(01)=1
c
c     define the working array lengths.
c     lw=079
c     liw=103
c     call dsplp(dusrmt,mrelas,nvars,costs,prgopt,dattrv,
c    *bl,bu,ind,info,primal,duals,ibasis,work,lw,iwork,liw)
c
c     calculate val, the minimal value of the objective function.
c     val=ddot(nvars,costs,1,primal,1)
c
c     stop
c     end
c    |------------------------|
c    |end of example of usage |
c    |------------------------|
c
c    |-------------------------------------|
c    |usage of dsplp( ) subprogram options.|
c    |-------------------------------------|
c
c     users frequently have a large variety of requirements for linear
c     optimization software.  allowing for these varied requirements
c     is at cross purposes with the desire to keep the usage of dsplp( )
c     as simple as possible. one solution to this dilemma is as follows.
c     (1) provide a version of dsplp( ) that solves a wide class of
c     problems and is easy to use. (2) identify parameters within
c     dsplp() that certain users may want to change.  (3) provide a
c     means of changing any selected number of these parameters that
c     does not require changing all of them.
c
c     changing selected parameters is done by requiring
c     that the user provide an option array, prgopt(*), to dsplp( ).
c     the contents of prgopt(*) inform dsplp( ) of just those options
c     that are going to be modified within the total set of possible
c     parameters that can be modified.  the array prgopt(*) is a linked
c     list consisting of groups of data of the following form
c
c          link
c          key
c          switch
c          data set
c
c     that describe the desired options.  the parameters link, key and
c     switch are each one word and are always required.  the data set
c     can be comprised of several words or can be empty.  the number of
c     words in the data set for each option depends on the value of
c     the parameter key.
c
c     the value of link points to the first entry of the next group
c     of data within prgopt(*).  the exception is when there are no more
c     options to change.  in that case, link=1 and the values for key,
c     switch and data set are not referenced.  the general layout of
c     prgopt(*) is as follows:
c          ...prgopt(1)=link1 (link to first entry of next group)
c          .  prgopt(2)=key1 (key to the option change)
c          .  prgopt(3)=switch1 (on/off switch for the option)
c          .  prgopt(4)=data value
c          .       .
c          .       .
c          .       .
c          ...prgopt(link1)=link2 (link to first entry of next group)
c          .  prgopt(link1+1)=key2 (key to option change)
c          .  prgopt(link1+2)=switch2 (on/off switch for the option)
c          .  prgopt(link1+3)=data value
c          ...     .
c          .       .
c          .       .
c          ...prgopt(link)=1 (no more options to change)
c
c     a value of link that is .le.0 or .gt. 10000 is an error.
c     in this case dsplp( ) returns with an error message, info=-14.
c     this helps prevent using invalid but positive values of link that
c     will probably extend beyond the program limits of prgopt(*).
c     unrecognized values of key are ignored.  if the value of switch is
c     zero then the option is turned off.  for any other value of switch
c     the option is turned on.  this is used to allow easy changing of
c     options without rewriting prgopt(*).  the order of the options is
c     arbitrary and any number of options can be changed with the
c     following restriction.  to prevent cycling in processing of the
c     option array prgopt(*), a count of the number of options changed
c     is maintained.  whenever this count exceeds 1000 an error message
c     (info=-15) is printed and the subprogram returns.
c
c     in the following description of the options, the value of
c     latp indicates the amount of additional storage that a particular
c     option requires.  the sum of all of these values (plus one) is
c     the minimum dimension for the array prgopt(*).
c
c     if a user is satisfied with the nominal form of dsplp( ),
c     set prgopt(1)=1 (or prgopt(1)=1.d0).
c
c     options:
c
c -----key = 50.  change from a minimization problem to a maximization
c     problem.
c     if switch=0  option is off; solve minimization problem.
c              =1  option is on; solve maximization problem.
c     data set =empty
c     latp=3
c
c -----key = 51.  change the amount of printed output.  the nominal form
c     of dsplp( ) has no printed output.
c     the first level of output (switch=1) includes
c
c     (1) minimum dimensions for the arrays costs(*),bl(*),bu(*),ind(*),
c         primal(*),duals(*),ibasis(*), and prgopt(*).
c     (2) problem dimensions mrelas,nvars.
c     (3) the types of and values for the bounds on x and w,
c         and the values of the components of the vector costs.
c     (4) whether optimization problem is minimization or
c         maximization.
c     (5) whether steepest edge or smallest reduced cost criteria used
c         for exchanging variables in the revised simplex method.
c
c     whenever a solution has been found, (info=1),
c
c     (6) the value of the objective function,
c     (7) the values of the vectors x and w,
c     (8) the dual variables for the constraints a*x=w and the
c         bounded components of x,
c     (9) the indices of the basic variables,
c    (10) the number of revised simplex method iterations,
c    (11) the number of full decompositions of the basis matrix.
c
c     the second level of output (switch=2) includes all for switch=1
c     plus
c
c    (12) the iteration number,
c    (13) the column number to enter the basis,
c    (14) the column number to leave the basis,
c    (15) the length of the step taken.
c
c     the third level of output (switch=3) includes all for switch=2
c     plus
c    (16) critical quantities required in the revised simplex method.
c          this output is rather voluminous.  it is intended to be used
c          as a diagnostic tool in case of a failure in dsplp( ).
c
c     if switch=0  option is off; no printed output.
c              =1  summary output.
c              =2  lots of output.
c              =3  even more output.
c     data set =empty
c     latp=3
c
c -----key = 52.  redefine the parameter, idigit, which determines the
c     format and precision used for the printed output.  in the printed
c     output, at least abs(idigit) decimal digits per number is
c     printed.  if idigit.lt.0, 72 printing columns are used.  if
c     idigit.gt.0, 133 printing columns are used.
c     if switch=0  option is off; idigit=-4.
c              =1  option is on.
c     data set =idigit
c     latp=4
c
c -----key = 53.  redefine lamat and lbm, the lengths of the portions of
c     work(*) and iwork(*) that are allocated to the sparse matrix
c     storage and the sparse linear equation solver, respectively.
c     lamat must be .ge. nvars+7 and lbm must be positive.
c     if switch=0  option is off; lamat=4*nvars+7
c                                 lbm  =8*mrelas.
c              =1  option is on.
c     data set =lamat
c               lbm
c     latp=5
c
c -----key = 54. redefine ipagef, the file number where the pages of the
c     sparse data matrix are stored.  ipagef must be positive and
c     different from isave (see option 56).
c     if switch=0  option is off; ipagef=1.
c              =1  option is on.
c     data set =ipagef
c     latp=4
c
c -----key = 55.  partial results have been computed and stored on unit
c     number isave (see option 56), during a previous run of
c     dsplp( ).  this is a continuation from these partial results.
c     the arrays costs(*),bl(*),bu(*),ind(*) do not have to have
c     the same values as they did when the checkpointing occurred.
c     this feature makes it possible for the user to do certain
c     types of parameter studies such as changing costs and varying
c     the constraints of the problem.  this file is rewound both be-
c     fore and after reading the partial results.
c     if switch=0  option is off; start a new problem.
c              =1  option is on; continue from partial results
c                  that are stored in file isave.
c     data set = empty
c     latp=3
c
c -----key = 56.  redefine isave, the file number where the partial
c     results are stored (see option 57).  isave must be positive and
c     different from ipagef (see option 54).
c     if switch=0  option is off; isave=2.
c              =1  option is on.
c     data set =isave
c     latp=4
c
c -----key = 57.  save the partial results after maximum number of
c     iterations, maxitr, or at the optimum.  when this option is on,
c     data essential to continuing the calculation is saved on a file
c     using a fortran binary write operation.  the data saved includes
c     all the information about the sparse data matrix a.  also saved
c     is information about the current basis.  nominally the partial
c     results are saved on fortran unit 2.  this unit number can be
c     redefined (see option 56).  if the save option is on,
c     this file must be opened (or declared) by the user prior to the
c     call to dsplp( ).  a crude upper bound for the number of words
c     written to this file is 6*nz.  here nz= number of nonzeros in a.
c     if switch=0  option is off; do not save partial results.
c              =1  option is on; save partial results.
c     data set = empty
c     latp=3
c
c -----key = 58.  redefine the maximum number of iterations, maxitr, to
c     be taken before returning to the user.
c     if switch=0  option is off; maxitr=3*(nvars+mrelas).
c              =1  option is on.
c     data set =maxitr
c     latp=4
c
c -----key = 59.  provide dsplp( ) with exactly mrelas indices which
c     comprise a feasible, nonsingular basis.  the basis must define a
c     feasible point: values for x and w such that a*x=w and all the
c     stated bounds on x and w are satisfied.  the basis must also be
c     nonsingular.  the failure of either condition will cause an error
c     message (info=-23 or =-24, respectively).  normally, dsplp( ) uses
c     identity matrix columns which correspond to the components of w.
c     this option would normally not be used when restarting from
c     a previously saved run (key=57).
c     in numbering the unknowns,
c     the components of x are numbered (1-nvars) and the components
c     of w are numbered (nvars+1)-(nvars+mrelas).  a value for an
c     index .le. 0 or .gt. (nvars+mrelas) is an error (info=-16).
c     if switch=0  option is off; dsplp( ) chooses the initial basis.
c              =1  option is on; user provides the initial basis.
c     data set =mrelas indices of basis; order is arbitrary.
c     latp=mrelas+3
c
c -----key = 60.  provide the scale factors for the columns of the data
c     matrix a.  normally, dsplp( ) computes the scale factors as the
c     reciprocals of the max. norm of each column.
c     if switch=0  option is off; dsplp( ) computes the scale factors.
c              =1  option is on; user provides the scale factors.
c     data set =scaling for column j, j=1,nvars; order is sequential.
c     latp=nvars+3
c
c -----key = 61.  provide a scale factor, costsc, for the vector of
c     costs.  normally, dsplp( ) computes this scale factor to be the
c     reciprocal of the max. norm of the vector costs after the column
c     scaling has been applied.
c     if switch=0  option is off; dsplp( ) computes costsc.
c              =1  option is on; user provides costsc.
c     data set =costsc
c     latp=4
c
c -----key = 62.  provide size parameters, asmall and abig, the smallest
c     and largest magnitudes of nonzero entries in the data matrix a,
c     respectively.  when this option is on, dsplp( ) will check the
c     nonzero entries of a to see if they are in the range of asmall and
c     abig.  if an entry of a is not within this range, dsplp( ) returns
c     an error message, info=-22. both asmall and abig must be positive
c     with asmall .le. abig.  otherwise,  an error message is returned,
c     info=-17.
c     if switch=0  option is off; no checking of the data matrix is done
c              =1  option is on; checking is done.
c     data set =asmall
c               abig
c     latp=5
c
c -----key = 63.  redefine the relative tolerance, tolls, used in
c     checking if the residuals are feasible.  normally,
c     tolls=relpr, where relpr is the machine precision.
c     if switch=0  option is off; tolls=relpr.
c              =1  option is on.
c     data set =tolls
c     latp=4
c
c -----key = 64. use the minimum reduced cost pricing strategy to choose
c     columns to enter the basis.  normally, dsplp( ) uses the steepest
c     edge pricing strategy which is the best local move.  the steepest
c     edge pricing strategy generally uses fewer iterations than the
c     minimum reduced cost pricing, but each iteration costs more in the
c     number of calculations done.  the steepest edge pricing is
c     considered to be more efficient.  however, this is very problem
c     dependent.  that is why dsplp( ) provides the option of either
c     pricing strategy.
c     if switch=0  option is off; steepest option edge pricing is used.
c              =1  option is on; minimum reduced cost pricing is used.
c     data set =empty
c     latp=3
c
c -----key = 65.  redefine mxitbr, the number of iterations between
c     recalculating the error in the primal solution.  normally, mxitbr
c     is set to 10.  the error in the primal solution is used to monitor
c     the error in solving the linear system.  this is an expensive
c     calculation and every tenth iteration is generally often enough.
c     if switch=0  option is off; mxitbr=10.
c              =1  option is on.
c     data set =mxitbr
c     latp=4
c
c -----key = 66.  redefine npp, the number of negative reduced costs
c     (at most) to be found at each iteration of choosing
c     a variable to enter the basis.  normally npp is set
c     to nvars which implies that all of the reduced costs
c     are computed at each such step.  this "partial
c     pricing" may very well increase the total number
c     of iterations required.  however it decreases the
c     number of calculations at each iteration.
c     therefore the effect on overall efficiency is quite
c     problem-dependent.
c
c     if switch=0 option is off; npp=nvars
c              =1 option is on.
c     data set =npp
c     latp=4
c
c -----key =  67.  redefine the tuning factor (phi) used to scale the
c     error estimates for the primal and dual linear algebraic systems
c     of equations.  normally, phi = 1.d0, but in some environments it
c     may be necessary to reset phi to the range 0.001-0.01.  this is
c     particularly important for machines with short word lengths.
c
c     if switch = 0 option is off; phi=1.d0.
c               = 1 option is on.
c     data set  = phi
c     latp=4
c
c -----key = 68.  used together with the subprogram dfulmt(), provided
c     with the dsplp() package, for passing a standard fortran two-
c     dimensional array containing the constraint matrix.  thus the sub-
c     program dfulmt must be declared in a fortran external statement.
c     the two-dimensional array is passed as the argument dattrv.
c     the information about the array and problem dimensions are passed
c     in the option array prgopt(*).  it is an error if dfulmt() is
c     used and this information is not passed in prgopt(*).
c
c     if switch = 0 option is off; this is an error is dfulmt() is
c                                  used.
c               = 1 option is on.
c     data set  = ia = row dimension of two-dimensional array.
c                 mrelas = number of constraint equations.
c                 nvars  = number of dependent variables.
c     latp = 6
c -----key = 69.  normally a relative tolerance (tolls, see option 63)
c     is used to decide if the problem is feasible.  if this test fails
c     an absolute test will be applied using the value tolabs.
c     nominally tolabs = zero.
c     if switch = 0 option is off; tolabs = zero.
c               = 1 option is on.
c     data set  = tolabs
c     latp = 4
c
c    |-----------------------------|
c    |example of option array usage|
c    |-----------------------------|
c     to illustrate the usage of the option array, let us suppose that
c     the user has the following nonstandard requirements:
c
c          a) wants to change from minimization to maximization problem.
c          b) wants to limit the number of simplex steps to 100.
c          c) wants to save the partial results after 100 steps on
c             fortran unit 2.
c
c     after these 100 steps are completed the user wants to continue the
c     problem (until completed) using the partial results saved on
c     fortran unit 2.  here are the entries of the array prgopt(*)
c     that accomplish these tasks.  (the definitions of the other
c     required input parameters are not shown.)
c
c     change to a maximization problem; key=50.
c     prgopt(01)=4
c     prgopt(02)=50
c     prgopt(03)=1
c
c     limit the number of simplex steps to 100; key=58.
c     prgopt(04)=8
c     prgopt(05)=58
c     prgopt(06)=1
c     prgopt(07)=100
c
c     save the partial results, after 100 steps, on fortran
c     unit 2; key=57.
c     prgopt(08)=11
c     prgopt(09)=57
c     prgopt(10)=1
c
c     no more options to change.
c     prgopt(11)=1
c     the user makes the call statement for dsplp( ) at this point.
c     now to restart, using the partial results after 100 steps, define
c     new values for the array prgopt(*):
c
c     again inform dsplp( ) that this is a maximization problem.
c     prgopt(01)=4
c     prgopt(02)=50
c     prgopt(03)=1
c
c     restart, using saved partial results; key=55.
c     prgopt(04)=7
c     prgopt(05)=55
c     prgopt(06)=1
c
c     no more options to change.  the subprogram dsplp( ) is no longer
c     limited to 100 simplex steps but will run until completion or
c     max.=3*(mrelas+nvars) iterations.
c     prgopt(07)=1
c     the user now makes a call to subprogram dsplp( ) to compute the
c     solution.
c    |--------------------------------------------|
c    |end of usage of dsplp( ) subprogram options.|
c    |--------------------------------------------|
c
c     |-----------------------------------------------|
c     |list of dsplp( ) error and diagnostic messages.|
c     |-----------------------------------------------|
c      this section may be required to understand the meanings of the
c     error flag =-info  that may be returned from dsplp( ).
c
c -----1. there is no set of values for x and w that satisfy a*x=w and
c     the stated bounds.  the problem can be made feasible by ident-
c     ifying components of w that are now infeasible and then rede-
c     signating them as free variables.  subprogram dsplp( ) only
c     identifies an infeasible problem; it takes no other action to
c     change this condition.  message:
c     dsplp( ). the problem appears to be infeasible.
c     error number =         1
c
c     2. one of the variables in either the vector x or w was con-
c     strained at a bound.  otherwise the objective function value,
c     (transpose of costs)*x, would not have a finite optimum.
c     message:
c     dsplp( ). the problem appears to have no finite soln.
c     error number =         2
c
c     3.  both of the conditions of 1. and 2. above have occurred.
c     message:
c     dsplp( ). the problem appears to be infeasible and to
c     have no finite soln.
c     error number =         3
c
c -----4.  the real and integer working arrays, work(*) and iwork(*),
c     are not long enough. the values (i1) and (i2) in the message
c     below will give you the minimum length required.  also redefine
c     lw and liw, the lengths of these arrays.  message:
c     dsplp( ). work or iwork is not long enough. lw must be (i1)
c     and liw must be (i2).
c               in above message, i1=         0
c               in above message, i2=         0
c     error number =        4
c
c -----5. and 6.  these error messages often mean that one or more
c     arguments were left out of the call statement to dsplp( ) or
c     that the values of mrelas and nvars have been over-written
c     by garbage.  messages:
c     dsplp( ). value of mrelas must be .gt.0. now=(i1).
c               in above message, i1=         0
c     error number =         5
c
c     dsplp( ). value of nvars must be .gt.0. now=(i1).
c               in above message, i1=         0
c     error number =         6
c
c -----7.,8., and 9.  these error messages can occur as the data matrix
c     is being defined by either dusrmt( ) or the user-supplied sub-
c     program, 'name'( ). they would indicate a mistake in the contents
c     of dattrv(*), the user-written subprogram or that data has been
c     over-written.
c     messages:
c     dsplp( ). more than 2*nvars*mrelas iters. defining or updating
c     matrix data.
c     error number =        7
c
c     dsplp( ). row index (i1) or column index (i2) is out of range.
c               in above message, i1=         1
c               in above message, i2=        12
c     error number =        8
c
c     dsplp( ). indication flag (i1) for matrix data must be
c     either 0 or 1.
c               in above message, i1=        12
c     error number =        9
c
c -----10. and 11.  the type of bound (even no bound) and the bounds
c     must be specified for each independent variable. if an independent
c     variable has both an upper and lower bound, the bounds must be
c     consistent.  the lower bound must be .le. the upper bound.
c     messages:
c     dsplp( ). independent variable (i1) is not defined.
c               in above message, i1=         1
c     error number =        10
c
c     dsplp( ).  lower bound (r1) and upper bound (r2) for indep.
c     variable (i1) are not consistent.
c               in above message, i1=         1
c               in above message, r1=    0.
c               in above message, r2=    -.1000000000e+01
c     error number =        11
c
c -----12. and 13.  the type of bound (even no bound) and the bounds
c     must be specified for each dependent variable.  if a dependent
c     variable has both an upper and lower bound, the bounds must be
c     consistent. the lower bound must be .le. the upper bound.
c     messages:
c     dsplp( ). dependent variable (i1) is not defined.
c               in above message, i1=         1
c     error number =        12
c
c     dsplp( ).  lower bound (r1) and upper bound (r2) for dep.
c      variable (i1) are not consistent.
c               in above message, i1=         1
c               in above message, r1=    0.
c               in above message, r2=    -.1000000000e+01
c     error number =        13
c
c -----14. - 21.  these error messages can occur when processing the
c     option array, prgopt(*), supplied by the user.  they would
c     indicate a mistake in defining prgopt(*) or that data has been
c     over-written.  see heading usage of dsplp( )
c     subprogram options, for details on how to define prgopt(*).
c     messages:
c     dsplp( ). the user option array has undefined data.
c     error number =        14
c
c     dsplp( ). option array processing is cycling.
c     error number =        15
c
c     dsplp( ). an index of user-supplied basis is out of range.
c     error number =        16
c
c     dsplp( ). size parameters for matrix must be smallest and largest
c     magnitudes of nonzero entries.
c     error number =        17
c
c     dsplp( ). the number of revised simplex steps between check-points
c     must be positive.
c     error number =        18
c
c     dsplp( ). file numbers for saved data and matrix pages must be
c     positive and not equal.
c     error number =        19
c
c     dsplp( ). user-defined value of lamat (i1)
c     must be .ge. nvars+7.
c               in above message, i1=         1
c     error number =         20
c
c     dsplp( ). user-defined value of lbm must be .ge. 0.
c     error number =         21
c
c -----22.  the user-option, number 62, to check the size of the matrix
c     data has been used.  an element of the matrix does not lie within
c     the range of asmall and abig, parameters provided by the user.
c     (see the heading: usage of dsplp( ) subprogram options,
c     for details about this feature.)  message:
c     dsplp( ). a matrix element's size is out of the specified range.
c     error number =        22
c
c -----23.  the user has provided an initial basis that is singular.
c     in this case, the user can remedy this problem by letting
c     subprogram dsplp( ) choose its own initial basis.  message:
c     dsplp( ). a singular initial basis was encountered.
c     error number =         23
c
c -----24.  the user has provided an initial basis which is infeasible.
c     the x and w values it defines do not satisfy a*x=w and the stated
c     bounds.  in this case, the user can let subprogram dsplp( )
c     choose its own initial basis.  message:
c     dsplp( ). an infeasible initial basis was encountered.
c     error number =        24
c
c -----25.subprogram dsplp( ) has completed the maximum specified number
c     of iterations.  (the nominal maximum number is 3*(mrelas+nvars).)
c     the results, necessary to continue on from
c     this point, can be saved on fortran unit 2 by activating option
c     key=57.  if the user anticipates continuing the calculation, then
c     the contents of fortran unit 2 must be retained intact.  this
c     is not done by subprogram dsplp( ), so the user needs to save unit
c     2 by using the appropriate system commands.  message:
c     dsplp( ). max. iters. (i1) taken. up-to-date results
c     saved on file (i2). if(i2)=0, no save.
c               in above message, i1=       500
c               in above message, i2=         2
c     error number =        25
c
c -----26.  this error should never happen.  message:
c     dsplp( ). moved to a singular point. this should not happen.
c     error number =        26
c
c -----27.  the subprogram la05a( ), which decomposes the basis matrix,
c     has returned with an error flag (r1).  (see the document,
c     "fortran subprograms for handling sparse linear programming
c     bases", aere-r8269, j.k. reid, jan., 1976, h.m. stationery office,
c     for an explanation of this error.)  message:
c     dsplp( ). la05a( ) returned error flag (r1) below.
c               in above message, r1=    -.5000000000e+01
c     error number =        27
c
c -----28.  the sparse linear solver package, la05*( ), requires more
c     space.  the value of lbm must be increased.  see the companion
c     document, usage of dsplp( ) subprogram options, for details on how
c     to increase the value of lbm.  message:
c     dsplp( ). short on storage for la05*( ) package. use prgopt(*)
c     to give more.
c     error number =        28
c
c -----29.  the row dimension of the two-dimensional fortran array,
c     the number of constraint equations (mrelas), and the number
c     of variables (nvars), were not passed to the subprogram
c     dfulmt().  see key = 68 for details.  message:
c     dfulmt() of dsplp() package. row dim., mrelas, nvars are
c     missing from prgopt(*).
c     error number =        29
c
c     |-------------------------------------------------------|
c     |end of list of dsplp( ) error and diagnostic messages. |
c     |-------------------------------------------------------|
c***references  r. j. hanson and k. l. hiebert, a sparse linear
c                 programming subprogram, report sand81-0297, sandia
c                 national laboratories, 1981.
c***routines called  dplpmn, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  corrected references to xerrwv.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dsplp
      double precision bl(*),bu(*),costs(*),dattrv(*),duals(*),
     * prgopt(*),primal(*),work(*),zero
c
      integer ibasis(*),ind(*),iwork(*)
      character*8 xern1, xern2
c
      external dusrmt
c
c***first executable statement  dsplp
      zero=0.d0
      iopt=1
c
c     verify that mrelas, nvars .gt. 0.
c
      if (mrelas.le.0) then
         write (xern1, '(i8)') mrelas
         call xermsg ('slatec', 'dsplp', 'value of mrelas must be ' //
     *      '.gt. 0.  now = ' // xern1, 5, 1)
         info = -5
         return
      endif
c
      if (nvars.le.0) then
         write (xern1, '(i8)') nvars
         call xermsg ('slatec', 'dsplp', 'value of nvars must be ' //
     *      '.gt. 0.  now = ' // xern1, 6, 1)
         info = -6
         return
      endif
c
      lmx=4*nvars+7
      lbm=8*mrelas
      last = 1
      iadbig=10000
      ictmax=1000
      ictopt= 0
c
c     look in option array for changes to work array lengths.
20008 next=prgopt(last)
      if (.not.(next.le.0 .or. next.gt.iadbig)) go to 20010
c
c     the checks for small or large values of next are to prevent
c     working with undefined data.
      nerr=14
      call xermsg ('slatec', 'dsplp',
     +   'the user option array has undefined data.', nerr, iopt)
      info=-nerr
      return
20010 if (.not.(next.eq.1)) go to 10001
      go to 20009
10001 if (.not.(ictopt.gt.ictmax)) go to 10002
      nerr=15
      call xermsg ('slatec', 'dsplp',
     +   'option array processing is cycling.', nerr, iopt)
      info=-nerr
      return
10002 continue
      key = prgopt(last+1)
c
c     if key = 53, user may specify lengths of portions
c    of work(*) and iwork(*) that are allocated to the
c     sparse matrix storage and sparse linear equation
c     solving.
      if (.not.(key.eq.53)) go to 20013
      if (.not.(prgopt(last+2).ne.zero)) go to 20016
      lmx=prgopt(last+3)
      lbm=prgopt(last+4)
20016 continue
20013 ictopt = ictopt+1
      last = next
      go to 20008
c
c     check length validity of sparse matrix staging area.
c
20009 if (lmx.lt.nvars+7) then
         write (xern1, '(i8)') lmx
         call xermsg ('slatec', 'dsplp', 'user-defined value of ' //
     *      'lamat = ' // xern1 // ' must be .ge. nvars+7.', 20, 1)
         info = -20
         return
      endif
c
c     trivial check on length of la05*() matrix area.
c
      if (.not.(lbm.lt.0)) go to 20022
      nerr=21
      call xermsg ('slatec', 'dsplp',
     +   'user-defined value of lbm must be .ge. 0.', nerr, iopt)
      info=-nerr
      return
20022 continue
c
c     define pointers for starts of subarrays used in work(*)
c     and iwork(*) in other subprograms of the package.
      lamat=1
      lcsc=lamat+lmx
      lcolnr=lcsc+nvars
      lerd=lcolnr+nvars
      lerp=lerd+mrelas
      lbasma=lerp+mrelas
      lwr=lbasma+lbm
      lrz=lwr+mrelas
      lrg=lrz+nvars+mrelas
      lrprim=lrg+nvars+mrelas
      lrhs=lrprim+mrelas
      lww=lrhs+mrelas
      lwork=lww+mrelas-1
      limat=1
      libb=limat+lmx
      librc=libb+nvars+mrelas
      lipr=librc+2*lbm
      liwr=lipr+2*mrelas
      liwork=liwr+8*mrelas-1
c
c     check array length validity of work(*), iwork(*).
c
      if (lw.lt.lwork .or. liw.lt.liwork) then
         write (xern1, '(i8)') lwork
         write (xern2, '(i8)') liwork
         call xermsg ('slatec', 'dsplp', 'work or iwork is not long ' //
     *      'enough. lw must be = ' // xern1 // ' and liw must be = ' //
     *      xern2, 4, 1)
         info = -4
         return
      endif
c
      call dplpmn(dusrmt,mrelas,nvars,costs,prgopt,dattrv,
     * bl,bu,ind,info,primal,duals,work(lamat),
     * work(lcsc),work(lcolnr),work(lerd),work(lerp),work(lbasma),
     * work(lwr),work(lrz),work(lrg),work(lrprim),work(lrhs),
     * work(lww),lmx,lbm,ibasis,iwork(libb),iwork(limat),
     * iwork(librc),iwork(lipr),iwork(liwr))
c
c     call dplpmn(dusrmt,mrelas,nvars,costs,prgopt,dattrv,
c    1 bl,bu,ind,info,primal,duals,amat,
c    2 csc,colnrm,erd,erp,basmat,
c    3 wr,rz,rg,rprim,rhs,
c    4 ww,lmx,lbm,ibasis,ibb,imat,
c    5 ibrc,ipr,iwr)
c
      return
      end
