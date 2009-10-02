#include "MIData.h"

#include <cmath>
#include <nongui/nonguilib.h>
#include <util/utillib.h>
#include <limits.h>
#include <string.h>

const std::string MIDatum::INVALID_STRING("INVALID_STRING");

MIDatum::MIDatum() {
  b = 0;
  i = INT_MIN;
  u = UINT_MAX;
  s = SHRT_MIN;
  f = FLT_MIN;
  d = DBL_MIN;
  radio = UINT_MAX;
  radio_count = UINT_MAX;
  isColor = false;
  isColorIndex = false;
  str = INVALID_STRING;
  strList.clear();
}

bool MIDatum::WriteDatum(std::string& dat_str) const {
  char buf[4096];

  if (i != INT_MIN) {
    sprintf(buf, "i %d; ", i);
  } else if (u != UINT_MAX) {
    sprintf(buf, "u %d; ", u);
  } else if (s != SHRT_MIN) {
    sprintf(buf, "s %hd; ", s);
  } else if (f != FLT_MIN) {
    sprintf(buf, "f %f; ", f);
  } else if (d != DBL_MIN) {
    sprintf(buf, "d %g; ", d);
  } else if (radio != UINT_MAX) {
    sprintf(buf, "r %d %d; ", radio, radio_count);
  } else if (strList.size() != 0) {
    //FIXME: handle strlist
  } else if (str != std::string("INVALID_STRING")) {
    sprintf(buf, "t%05d %s; ", str.size(), str.c_str());
  } else if (isColor) {
    sprintf(buf, "c %d %d %d; ", color[0], color[1], color[2]);
  } else if (isColorIndex) {
    sprintf(buf, "C %d; ", color[0]);
  } else { // assume boolean
    sprintf(buf, "b %d; ", b);
  }
  dat_str += std::string(buf);
  return true;
}


unsigned int MIDatum::ReadDatum(const std::string& dat) {
  char buf[1024];
  int intColor[3];
  memset(buf, 0, 1024);

  // scan char-by-char until we get semicolon
  unsigned int count;
  for (count = 0; count < dat.size(); ++count) {
    if (dat[count] == ';') {
      ++count; // skip trailing space
      break;
    }
    buf[count] = dat[count];
  }

  int btmp;
  switch (buf[0]) {
    case 'i':
      if (sscanf(buf, "i %d", &i) == 1) {
        return count;
      }
      break;
    case 'u':
      if (sscanf(buf, "u %d", &u) == 1) {
        return count;
      }
      break;
    case 's':
      if (sscanf(buf, "s %hd", &s) == 1) {
        return count;
      }
      break;
    case 'f':
      if (sscanf(buf, "f %f", &f) == 1) {
        return count;
      }
      break;
    case 'd':
      if (sscanf(buf, "d %lf", &d) == 1) {
        return count;
      }
      break;
    case 'r':
      if (sscanf(buf, "r %d %d", &radio, &radio_count) == 2) {
        return count;
      }
      break;
    case 'c':
      if (sscanf(buf, "c %d %d %d", &intColor[0], &intColor[1], &intColor[2]) == 3) {
        color[0] = intColor[0];
        color[1] = intColor[1];
        color[2] = intColor[2];
        isColor=true;
        return count;
      }
      break;
    case 'C':
      if (sscanf(buf, "C %d", &intColor[0]) == 1) {
        color[0] = intColor[0];
        isColorIndex=true;
        return count;
      }
      break;
    case 'b':
      if (sscanf(buf, "b %d", &btmp) == 1) {
        b = btmp != 0; return count;
      }
      break;
    case 't':
      {
        int len;
        if (!sscanf(buf, "t%05d", &len) == 1) {
          return 0;
        }
        char tmp2[1024], fmt[64];
        memset(tmp2, 0, 1024);
        sprintf(fmt, "%%%dc", len);
        sscanf(&buf[7], fmt, tmp2);
        str = std::string(tmp2);
        return count;
        break;
      }
      //FIXME: handle strlist
  }
  return 0;
}



std::string MIDataToString(MIData& data) {
  std::string tmp,str("(");

  for (MIDataConstIter i = data.begin(); i != data.end(); ++i) {
    tmp = i->first + ".";
    if (!i->second.WriteDatum(tmp)) {
      return false;
    }
    str += tmp;
  }
  str+=")";
  return str;
}

bool StringToMIData(const std::string& instr, MIData& data, std::string*) {
  if (!instr.size()) {
    return false;
  }

  unsigned int count = 1;
  std::string str(&instr[1]); // chops off first '('
  data.clear();

  while (count < instr.size()) {
    std::string::size_type pos = str.find_first_of('.');
    if (pos == std::string::npos) {
      return true; // got to end of input, all that remains is trailing ')' or garbage
    }
    std::string name(str, 0, pos);
    count += (pos+1); // skip past name and '.'
    MIDatum d;
    unsigned int tmp_cnt = d.ReadDatum(&str[pos+1]);
    if (!tmp_cnt) {
      return false;
    }
    data[name] = d;
    count += (tmp_cnt+1);
    str = std::string(&instr[count]);
  }
  return true;
}

