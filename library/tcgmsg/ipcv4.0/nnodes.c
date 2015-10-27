/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/nnodes.c,v 1.1 1994/02/23 16:54:44 d3g681 Exp $ */

#include "sndrcv.h"
#include "sndrcvP.h"

long NNODES_()
/*
  return total no. of processes
*/
{
  return SR_n_proc;
}

