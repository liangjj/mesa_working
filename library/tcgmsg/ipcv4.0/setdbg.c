/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/setdbg.c,v 1.1 1994/02/23 16:55:01 d3g681 Exp $ */

#include "sndrcv.h"
#include "sndrcvP.h"

void SETDBG_(value)
    long *value;
/*
  set global debug flag for this process to value
*/
{
  SR_debug = *value;
}

