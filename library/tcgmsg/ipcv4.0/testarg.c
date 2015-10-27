/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/testarg.c,v 1.1 1994/02/23 16:55:24 d3g681 Exp $ */

/*
  This checks the functioning of the include file farg.h
*/

#include "farg.h"

void parg()
{
  int i;

  for (i=0; i<ARGC_; i++)
    (void) printf("argv(%d)=%s\n", i, ARGV_[i]);
}
