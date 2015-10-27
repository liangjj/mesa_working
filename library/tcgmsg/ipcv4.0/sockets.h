/* $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/sockets.h,v 1.1 1994/02/23 16:55:14 d3g681 Exp $ */

extern void ShutdownAll();
extern int ReadFromSocket();
extern int WriteToSocket();
extern void CreateSocketAndBind();
extern int ListenAndAccept();
extern int CreateSocketAndConnect();
extern long PollSocket();
