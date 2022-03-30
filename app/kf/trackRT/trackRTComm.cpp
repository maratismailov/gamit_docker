/*
   trackRT real-time communication module.
   This code is derived from the test_bnc_qt.cpp developed by
   the German Federal Agency for Cartography and Geodesy (BKG)

   This module receives data through the SocketObs object.  Once
   and epoch of data is available, the data is passed to the trackRT
   processing.

   Vers 1.1: Data on socket compatable with BNC ver 2.5 ASCII streams

   Copyright (C) MIT 2009.

*/
// Part of BNC, a utility for retrieving decoding and
// converting GNSS data streams from NTRIP broadcasters.
//
// Copyright (C) 2007
// German Federal Agency for Cartography and Geodesy (BKG)
// http://www.bkg.bund.de
// Czech Technical University Prague, Department of Geodesy
// http://www.fsv.cvut.cz
//
// Email: euref-ip@bkg.bund.de
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, version 2.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

/* -------------------------------------------------------------------------
 * BKG NTRIP Client
 * -------------------------------------------------------------------------
 *
 * Class:      test_bnc_qt
 *
 * Purpose:    Example program to read BNC output from IP port.
 *
 * Author:     L. Mervart
 *
 * Created:    24-Jan-2007
 *
 * Changes:    
 *
 * -----------------------------------------------------------------------*/

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>

#include "GPSDecoder.h"

#include <QTcpSocket>
extern "C" {
   void saveobsa_( int*, int*, int*, int [], char* );
   void procobsrt_( int* , int* ); 
   void reportobs_( int* );
/* C-routine to read command line: Results are then passed to trt_cpcoml*/
   void trt_comline( int argc, char *argv[], int *port, char machine[],
                     char trt_cmdfile[], char trt_prtroot[], char refsite[], 
                     int *nums, char prcsite[] ) ; 
/* Fortran routine decode command values and save them in the Fortran common */
   void trt_cpcoml_( int* argc,  int* port, char machine[],
                     char trt_cmdfile[],  char trt_prtroot[], char refsite[], 
                     int *nums, char prcsite[] ) ;
/* Fortran Routine to read command file */
   void read_batrt_( int*, int [] ) ; 
} 

using namespace std;

QTcpSocket* sockconnect( int ) ; /* Prototype */

int main(int argc, char* argv[]) {


  char machine[80] = "127.0.0.1" ;  // Name of machine with stream -m option
  int  port = 3765 ;          // Port number on machine -p option (default for bnc set) 
  char trt_cmdfile[256] = "\0" ; /* Name of the trackRT command file
                        passed with the -f option */
  char trt_prtroot[256] = "\0" ; /* Name of the root print name
                        passed with the -n option */
  int  nums = 0 ;  /* Number of sites to be processed */
  int  ep = 0 ;    /* Epoch number counter */ 
  int  OK ;        /* Counter used if socket connection lost (counts to 60s) */
  int  rerun = 0 ; /* Set zero for first call to read_batRT */
  int  pass_debug[10] ;  /* Debug array passed back from read_batrt */ 
  char refsite[5] = "\0"   ; /* Name of reference site -r option */
  char prcsite[4*20+1] = "\0" ; /* Names of upto 20 sites (max_site in trackRT.h -s option */

/* Get the command line for trackRT.  Most values are saved
  in the fortran common but the machine and port are returned
  here */
  trt_comline( argc, argv, &port, machine,
                trt_cmdfile, trt_prtroot, refsite, &nums, prcsite ) ;

/*
  printf("\nAfter TrackRT %s Port %d\n", machine, port);
  printf("Track command file: %s \n",trt_cmdfile);
  printf("Reference Site    : %s \n",refsite);
  printf("Sites %2d Names    : %s \n", nums, prcsite);
*/

/* Now copy results to Fortran common block */
  trt_cpcoml_( &argc, &port, machine,
              trt_cmdfile, trt_prtroot, refsite, &nums, prcsite ) ;

/* Read the command file (all information passed through common */
  read_batrt_( &rerun, pass_debug );

  QTcpSocket socketObs;

// Connect to the socket
//  socketObs.connectToHost("127.0.0.1", port);
  socketObs.connectToHost(machine, port);

  if (!socketObs.waitForConnected(10000)) {
    cerr << "socketObs: not connected on port " << port << endl;
//  Wait and try again to see of port returns
    OK = 0;
    while ( OK < 60 ) {
       sleep(1);
       OK++;
       socketObs.connectToHost(machine, port); // Try again
       if (socketObs.waitForConnected(10000)) {
          cerr << "trackRT: socketObs: Port start up " << OK << " sec" << endl;
          OK = 100;  // This will cause loop to exit and continue
       }
    }
    if( OK < 100 ) {
//      Socket connection failed
        cerr << "socketObs: Fail after 60-sec port " << port << endl;
        exit(1);
    }
  }


  cout.setf(ios::fixed);


  int numrt=0 ; // Number of real time values
  int Socketcnt = 0 ; // Number of connection tries 

  QByteArray buffer;
  const int maxbytes = 363000 ;  // Must be same in trackRTObs.h
  char  obsstring[maxbytes]; 
  int   nbytes ;   // Number of bytes read from socket
  int   endepoch ;  // Returns == 1 when epoch complete and can be 
                    // processed; otherwise zero.

  while (true) {
    if (socketObs.state() != QAbstractSocket::ConnectedState) {
      cerr << "trackRT: socketObs: disconnected: Wait 60 secs" << endl;
//    Wait and try again to see of port returns
      OK = 0;
      while ( OK < 60 ) {
         sleep(1);
         OK++;
         Socketcnt = 0; 
         socketObs.connectToHost(machine, port); // Try again
         socketObs.waitForConnected(10000); // Wait for connection
         if( socketObs.state() == QAbstractSocket::ConnectedState) {
             cerr << "trackRT: socketObs: Reconnected " << OK << " sec" << endl;
             OK = 100;  // This will cause loop to exit and continue
         }
      }
      if( OK < 100 ) {
//      Socket connection failed
        cerr << "socketObs: Reconnect failed after 60-sec port " << port << endl;
        exit(1);
      }
    }

    nbytes = socketObs.bytesAvailable();
    if( nbytes > maxbytes ) {
       cerr << "trackRT buffer overflow " << maxbytes << " maxbytes exceeded" << endl ;
       exit(2);
    }

    if ( nbytes ) {

       buffer = socketObs.readAll();

//     printf("Bytes %d Size %d \n",nbytes, buffer.size() );
//     Copy over to a string array for processing in Fortran.
       strcpy( obsstring, buffer.constData());

       saveobsa_( &nbytes, &endepoch, &numrt, pass_debug, obsstring) ;

      if( endepoch == 1 ) {
//        The epoch has just changed for go process current data
//        block.  When done, decode number the rest of data in this
//        block
//        printf("Processing Data Block with %d numrt\n",numrt);
          if( pass_debug[9] > 0 ) { 
               reportobs_( &numrt );
          }

          procobsrt_(&numrt, &ep) ;

//        Now save the rest of the block
          numrt = 0 ;
//        saveobsa_( &nbytes, &endepoch, &numrt, pass_debug, obsstring) ;
      }
    }
    else {
      socketObs.waitForReadyRead(1);
    }
  }

  return 0;
}

/* Routine to read command line.  Problems when this is tried directly in
   Fortran so we decode here and then pass to Fortran routine to save the 
   values in common */
void trt_comline( int argc, char *argv[], int *port, char machine[],
                   char trt_cmdfile[], char trt_prtroot[], char refsite[], 
                   int *nums, char prcsite[] ) {

/* Function to read the trackRT command line and return the values to the
   main program.  The values are then copied into the fotran common block */

int  i, done;    // Loop, indicator that commandline is fully read
int  indx ;      // Counter for position of site names
char coml[256] ; // Allow upto 256 characters for file names

/* Loop over the command line arguments decoding the arguments */
  for( i=1 ; i < argc ; i++ ) {
//    copy argument into coml string 
      strcpy(coml,argv[i]);
//    Now start testing
      if( strncmp(coml,"-p",2) == 0 ) {
//        Port option; increment i and get port 
          i++ ; *port = atoi(argv[i]); }
      else if( strncmp(coml,"-m",2) == 0 ) {
//        Get the machine name
          i++; strcpy(machine, argv[i]);}
      else if ( strncmp(coml,"-f",2) == 0 ) {
//        Get name of trackRT commamnd file
          i++; strcpy(trt_cmdfile, argv[i]);}
      else if ( strncmp(coml,"-n",2) == 0 ) {
//        Get name of root to print file 
          i++; strcpy(trt_prtroot, argv[i]);}
      else if ( strncmp(coml,"-r",2) == 0 ) {
//        Get 4-char name reference site
          i++; strncpy(refsite, argv[i],5);} 
      else if ( strncmp(coml,"-d",2) == 0 ) {
/*        Get 4-character codes of sites to be processed
          Here we read the 4-char entries until the end of the
          command line or the next option is found */
          *nums = 0; done = 0 ; 
          while ( done == 0 ) {
/*           Read list of sites until end of arguments or
             another - is found */
             i ++ ; 
             if ( strncmp(argv[i],"-",1) == 0 ) {
                i--;
                done = 1; }
             else {
                indx = *nums * 4; 
                strncpy(prcsite+indx,argv[i],4);
                *nums = *nums + 1; 
                if( i+1 == argc ) {done = 1;}
             }
             indx = *nums * 4; 
             prcsite[indx] = NULL ;  // Null terminate the list
          }
       }
    }
}







 
