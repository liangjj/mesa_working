/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/msgtypesc.h,v 1.1 1994/02/23 16:54:41 d3g681 Exp $ */

#ifndef MSGTYPES_H_
#define MSGTYPES_H_

/* 
   This defines bit masks that can be OR'ed with user types (1-32767)
   to indicate the nature of the data to the message passing system
*/ 
#ifdef IPSC
#define MSGDBL 0
#define MSGINT 0
#define MSGCHR 0
#else
#define MSGDBL  65536
#define MSGINT 131072
#define MSGCHR 262144
#endif

#endif
