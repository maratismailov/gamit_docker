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

#ifndef GPSDECODER_H
#define GPSDECODER_H

#include <iostream>
#include <vector>
#include <string>
#include <QPointer>
#include <QList>
#include <QStringList>

//#include "bncconst.h"

class t_obsInternal {
 public:

  t_obsInternal() : 
    Head(0) {
    StatInf[0] = '\x0';
  }
  char   StatInf[363];  // Station ID

  int    Head;          // L2 signal-to noise ratio (mapped to integer)
};

class t_obs : public QObject{
 public:
  enum t_obs_status {initial, posted, received};

  t_obs() {
    _status = initial;

    _o.StatInf[0]     = '\0';
    _o.Head          = 0;
  }

  ~t_obs() {}

  t_obsInternal _o;
  t_obs_status  _status;
};

/*
typedef QPointer<t_obs> p_obs;

class GPSDecoder {
 public:
  virtual t_irc Decode(char* buffer, int bufLen, std::vector<std::string>& errmsg) = 0;

  virtual ~GPSDecoder() {
    QListIterator<p_obs> it(_obsList);
    while (it.hasNext()) {
      p_obs obs = it.next();
      if (!obs.isNull() && obs->_status == t_obs::initial) {
        delete obs;
      }
    }
  }

  virtual int corrGPSEpochTime() const {return -1;}

  struct t_antInfo {
    enum t_type { ARP, APC };

    t_antInfo() {
      xx = yy = zz = height = 0.0;
      type = ARP;
      height_f = false;
      message  = 0;
    };

    double xx;
    double yy;
    double zz;
    t_type type;
    double height;
    bool   height_f;
    int    message;
  };

  QList<p_obs>     _obsList;
  QList<int>       _typeList;  // RTCM   message types
  QStringList      _antType;   // RTCM   antenna descriptor
  QList<t_antInfo> _antList;   // RTCM   antenna XYZ
};
*/

#endif
