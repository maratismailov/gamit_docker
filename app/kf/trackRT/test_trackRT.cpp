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
 * Class:      test_trackRT
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
#include <iostream>
#include <iomanip>

//#include "GPSADecode.h"
//#include "bnctime.h"
extern "C" {
//#include "RTCM3/rtcm3torinex.h"
}
extern "C" {
   void saveobsa_( int*, int*, int*, int [], char* );
   void reportobs_( int* );
}
#include <QTcpSocket>

using namespace std;

int main(int argc, char**  argv ) {

// See if Port passed
  if (argc < 2) {
    cerr << "Usage: test_trackRT port [host]\n";
    exit(1);
  }

  int port = atoi(argv[1]);

// See if machine passed; Set default and then decode command line argument
  char machine[80] = "127.0.0.1" ;  // Name of machine with stream -m option
  if ( argc == 3 ) {
      strcpy(machine, argv[2]);
  }
  cout << endl << "test_trackRT port connection debug program" << endl;
  cout << "Connecting to port " << port << " on host " << machine << endl;

  QTcpSocket socketObs;

  socketObs.connectToHost(machine, port);
//  socketObs.connectToHost("caputo.mit.edu", port);
  if (!socketObs.waitForConnected(10000)) {
    cerr << "socketObs: not connected on port " << port << endl;
    exit(1);
  }

  cout.setf(ios::fixed);


  QByteArray buffer;
  int maxbytes = 363000 ;  // Must be same in trackRTObs.h
  char obsstring[maxbytes]; 
  int nbytes, endepoch, nobs ;
  int pass_debug[10] ; /* Debug array; see to output SaveObsA debug */

  nobs = 0;
  pass_debug[9] = 1 ;

  while (true) {
    if (socketObs.state() != QAbstractSocket::ConnectedState) {
      cerr << "socketObs: disconnected" << endl;
      exit(1);
    }
    nbytes = socketObs.bytesAvailable();
    if( nbytes > maxbytes ) {
       cerr << "Buffer overflow " << maxbytes << " maxbytes exceeded" << endl ;
       exit(2);
    }

    if ( nbytes ) {
      buffer = socketObs.readAll();
      printf("Bytes %d Nobs %d Size %d \n",nbytes, nobs, buffer.size() );
//    Copy over to a string array for processing in Fortran.
      strcpy( obsstring, buffer.constData());
      printf("Rec: |%s",obsstring);

//      for ( k=0 ; k < nbytes ; k++ ) {printf("%c",obsstring[k]);}
//      printf("\n");

//      printf("End %d\n",obsstring[nbytes-1]);

//    Now try a fortran call to see what we find
      saveobsa_( &nbytes, &endepoch, &nobs, pass_debug, obsstring) ;

      if( endepoch == 1 ) {
//        The epoch has just changed for go process current data
//        block.  When done, decode number the rest of data in this
//        block
          printf("Processing Data Block with %d nobs\n",nobs);
          reportobs_( &nobs );
          nobs = 0 ;
          printf("------------------------------------------------------");
          printf("------------------------------------------------------\n");
//        printf("Remaining record size nbytes %d \n",nbytes) ;
//        saveobsa_( &nbytes, &endepoch, &nobs, pass_debug, obsstring) ;

      }
    }
    else {
      socketObs.waitForReadyRead(1);
    }
  }

  return 0;
}
