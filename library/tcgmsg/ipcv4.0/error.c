/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/error.c,v 1.1 1994/02/23 16:54:16 d3g681 Exp $ */

#include <stdio.h>
#include <setjmp.h>
#include <signal.h>

#include "sndrcvP.h"
#include "sndrcv.h"
#include "signals.h"
#include "sockets.h"

#ifdef SHMEM
#include "sema.h"
#include "shmem.h"
#endif

extern jmp_buf SR_jmp_buf;   /* Jumped to on soft error */ 

#ifdef CRAY
#include <errno.h>
#else
extern int errno;
#endif

extern void exit();
extern int SR_caught_sigint;

void Error(string, integer)
     char *string;
     long integer;
{
  (void) signal(SIGCHLD, SIG_DFL); /* Death of children to be expected */
  (void) signal(SIGINT, SIG_IGN);

  (void) fflush(stdout);
  if (SR_caught_sigint) {
    (void) fprintf(stderr,"%2ld: interrupt\n",NODEID_());
  }
  else {
    (void) fprintf(stderr,"%2ld: %s %ld (%#lx).\n", NODEID_(), string,
		   integer,integer);
    if (errno != 0)
      perror("system error message");
    if (DEBUG_)
      PrintProcInfo();
  }
  (void) fflush(stderr);

  /* Shut down the sockets and remove shared memory and semaphores to
     propagate an error condition to anyone that is trying to communicate
     with me */

  ZapChildren();  /* send interrupt to children which should trap it
		     and call Error in the handler */

#ifdef SHMEM
#if defined(NOSPIN)
  (void) SemSetDestroyAll();
#endif
  (void) DeleteSharedRegion(SR_proc_info[NODEID_()].shmem_id);
#endif
  ShutdownAll();    /* Close sockets for machines with static kernel */

  if (SR_exit_on_error)
    exit(1);
  else {
    SR_error = 1;
    (void) longjmp(SR_jmp_buf, 1); /* For NXTVAL server */
  }
}

void PARERR_(code)
   long *code;
/*
  Interface from fortran to c error routine
*/
{
  Error("User detected error in FORTRAN", *code);
}
