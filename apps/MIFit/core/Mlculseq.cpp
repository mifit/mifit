#define NOMINMAX
#include <cmath>
#include <climits>
#include <cstring>

#include "nonguilib.h"
#include "mathlib.h"
#include "chemlib.h"

#include "Molecule.h"
#include "seqtypes.h"
#include "RESIDUE.h"
#include "RESIDUE_.h"

using namespace chemlib;


int Molecule::SeqMax() {
  //int seqmax = 0;
  MIIter<RESIDUE> ri = GetResidues();
  ri.Last();
  RESIDUE* res = ri;
  if (!Residue::isValid(res)) {
    return 0;
  }
  return std::max((size_t)res->seqpos(), alt_seq.length());
}

char Molecule::GetSeq(int index) {
  if (index >= (int)alt_seq.length()) {
    return ' ';
  } else {
    return alt_seq.c_str()[index];
  }
}

void Molecule::SetSeq(char /* aa */, int index) {
  if (index < (int)alt_seq.length()) {
    alt_seq = MIToUpper(alt_seq);
    //if(islower(aa)) aa = toupper(aa);
    //alt_seq.SetAt(index, aa);
  }
}

std::string Molecule::SeqString() {
  return alt_seq;
}

void Molecule::SetSequence(std::string s) {
  int i /*, iseq=0*/;
  char c;
  alt_seq = "";
  for (i = 0; i < (int)s.length(); i++) {
    c = s.c_str()[i];
    if (isalpha(c) || c == '.') {
      if (islower(c)) {
        c = toupper(c);
      }
      alt_seq += c;
    }
  }
}

void Molecule::ReadSequence(std::string path, int /* type */, int /* skiplines */) {
  FILE* fp = fopen(path.c_str(), "r");
  if (!fp) {
    Logger::message("Cannot open sequence file");
  }
  std::string file;
  char buf[1024];
  while (fgets(buf, sizeof buf, fp) != NULL) {
    if (!(buf[0] == '#' || buf[0] == '<' || buf[0] == ';' || buf[0] == '>')) {
      file += buf;
    }
  }
  SetSequence(file);
  fclose(fp);
}

int Molecule::WriteSequence(std::string path, int type) {
  FILE* fp = fopen(path.c_str(), "w");
  char c;
  if (!fp) {
    Logger::message("Cannot open sequence file");
    return 0;
  }
  fprintf(fp, ">%s\n\n", compound.c_str());
  //std::string s = SeqString();
  int l = alt_seq.length();
  int i;
  if (type == SEQ_FORMAT_SIMPLE) {
    for (i = 0; i < l; i++) {
      c = alt_seq.c_str()[i];
      fprintf(fp, "%c", alt_seq.c_str()[i]);
      if (i > 0 && i%10 == 0 && i%100 != 0) {
        fprintf(fp, " ");
      }
      if (i > 0 && i%100 == 0) {
        fprintf(fp, "\n");
      }
    }
  } else if (type == SEQ_FORMAT_MODEL) {
    i = 0;
    for (MIIter<RESIDUE> res = GetResidues(); Residue::isValid(res); ++res) {
      c = (char)res->name1();
      if (!(c == '?' || c == 'B' || c == 'J' || c == 'O' || c == 'U' || c == 'X' || c == 'Z') && isupper(c)) {
        fprintf(fp, "%c", c);
        if (i > 0 && i%10 == 0 && i%100 != 0) {
          fprintf(fp, " ");
        }
        if (i > 0 && i%100 == 0) {
          fprintf(fp, "\n");
        }
        i++;
      }
    }
  } else if (type == SEQ_FORMAT_MSF) {
    int nout = 0;
    for (i = 0; i < l; i++) {
      c = toupper(alt_seq.c_str()[i]);
      if (!(c == '?' || c == 'B' || c == 'J' || c == 'O' || c == 'U' || c == 'X' || c == 'Z') && isalpha(c)) {
        fprintf(fp, "%c", alt_seq.c_str()[i]);
        //if(nout>0 && nout%10==0 && nout%100!=0)fprintf(fp," ");
        if (nout > 0 && nout%80 == 0) {
          fprintf(fp, "\n");
        }
        nout++;
      }
    }
  }
  fprintf(fp, "\n");
  fclose(fp);
  return 1;
}

void Molecule::InsertGap(RESIDUE* gap_point) {
  if (gap_point == NULL) {
    return;
  }
  MIIter<RESIDUE> res = GetResidues();
  while (Residue::isValid(res)) {
    if (gap_point == res) {
      while (res) {
        res->setSeqpos(res->seqpos() + 1);
        ++res;
      }
      return;
    }
    ++res;
  }
}

void Molecule::DeleteGap(RESIDUE* gap_point) {
  if (gap_point == NULL) {
    return;
  }
  MIIter<RESIDUE> res = GetResidues();
  RESIDUE* prev = NULL;
  while (Residue::isValid(res)) {
    if (gap_point == res) {
      if (prev) {
        if (prev->seqpos() < gap_point->seqpos()-1) {
          while (res) {
            res->setSeqpos(res->seqpos()-1);
            ++res;
          }
        }
      } else if (res == residues && res->seqpos() > 0) {
        while (res) {
          res->setSeqpos(res->seqpos()-1);
          ++res;
        }
      }
      return;
    }
    prev = res;
    ++res;
  }
}

void Molecule::InsertLowerGap(RESIDUE* gap_point) {
  if (!gap_point) {
    return;
  }
  int p = gap_point->seqpos();
  int l = alt_seq.length();
  if (p >= l) {
    return;
  }
  alt_seq.insert(p,".");
}

void Molecule::DeleteLowerGap(RESIDUE* gap_point) {
  if (!gap_point) {
    return;
  }
  int p = gap_point->seqpos();
  int l = alt_seq.length();
  if (p >= l) {
    return;
  }
  if (alt_seq.c_str()[p] != '.') {
    return;
  }
  alt_seq.erase(p,1);
}

int Molecule::SequenceIdentities() {
  int n = 0;
  int l = alt_seq.length();
  if (l > 0) {
    for (MIIter<RESIDUE> res = GetResidues(); Residue::isValid(res); ++res) {
      if (res->seqpos() < (unsigned int)l) {
        if (alt_seq[res->seqpos()] == res->name1() && isalpha(res->name1())) {
          n++;
        }
      }
    }
  }
  return n;
}

