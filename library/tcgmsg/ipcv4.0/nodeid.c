/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/nodeid.c,v 1.1 1994/02/23 16:54:45 d3g681 Exp $ */

#include "sndrcv.h"
#include "sndrcvP.h"

long NODEID_()
/*
  return logical node no. of current process
*/
{
  return SR_proc_id;
}
