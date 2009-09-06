#include "MIBitVector.h"
#include <string.h>

MIBitVector::MIBitVector(const std::string& hex_data) {
  if (hex_data.size()) {
    FromHexString(hex_data);
  }
}

void MIBitVector::Set(size_t i, bool state) {
  if (i >= GetSize()) {
    Resize(i+1);
  }
  _data[i] = state;
}

bool MIBitVector::IsSet(size_t i) const {
  return (i < GetSize() && _data[i] != 0);
}

void MIBitVector::Resize(size_t size) {
  _data.resize(size);
}

size_t MIBitVector::GetSize() const {
  return _data.size();
}

const char* hexDigits = "0123456789ABCDEF";

std::string MIBitVector::ToHexString() const {
  std::string hex_data;
  size_t last = GetSize(), size = last;

  //NOTE: resultant string will be zero-padded at the end
  if (last%8 != 0) {
    last += 8-(last%8);
  }

  for (size_t i = 0; i < last; i += 8) {
    unsigned char ch = 0;

    for (unsigned int j = 0; j < 8; ++j) {
      if (i+j < size) {
        ch |= (_data[i+j] << j);
      }
    }
    hex_data.push_back(hexDigits[(ch   ) & 0x0f]);
    hex_data.push_back(hexDigits[(ch>>4) & 0x0f]);
  }
  return hex_data;
}

bool MIBitVector::FromHexString(const std::string& hex_data) {
  _data.clear();
  if (!hex_data.size()) {
    return true;
  }

  _data.reserve(hex_data.size()*4);

  for (size_t i = 0; i < hex_data.size(); ++i) {
    switch (hex_data[i]) {
      default:
        return false;
        break;
      case '0':
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(0);
        break;
      case '1':
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(0);
        break;
      case '2':
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(0);
        break;
      case '3':
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(0);
        break;
      case '4':
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(0);
        break;
      case '5':
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(0);
        break;
      case '6':
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(0);
        break;
      case '7':
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(0);
      case '8':
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(1);
        break;
      case '9':
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(1);
        break;
      case 'A':
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(1);
        break;
      case 'B':
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(1);
        break;
      case 'C':
        _data.push_back(0);
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(1);
        break;
      case 'D':
        _data.push_back(1);
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(1);
        break;
      case 'E':
        _data.push_back(0);
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(1);
        break;
      case 'F':
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(1);
        _data.push_back(1);
        break;
    }
  }
  return true;
}

std::string MIBitVector::ToBinaryString() const {
  std::string bin_str;
  for (size_t i = 0; i < _data.size(); ++i) {
    bin_str.push_back(hexDigits[_data[i]]);
  }
  return bin_str;
}

size_t MIBitVector::LastBitSet() const {
  if (GetSize() == 0) {
    return 0;
  }

  for (size_t last = GetSize()-1; last != 0; --last) {
    if (IsSet(last)) {
      return last;
    }
  }
  return (size_t)IsSet(0);
}

bool MIBitVector::IsEqual(const MIBitVector& other) const {
  size_t this_last = LastBitSet();
  size_t other_last = other.LastBitSet();

  if (this_last != other_last) {
    return false;
  }

  if (this_last == 0) {
    return true;
  }

  return (memcmp(&_data[0], &other._data[0], this_last) == 0);
}

//#define DO_MIBITVECTOR_TESTS
#ifdef DO_MIBITVECTOR_TESTS
void DoTest(const std::string& instr) {
  std::string outstr;
  MIBitVector bv(instr);
  outstr = bv.ToHexString();
  printf("Does %s = %s ? %d\n", instr.c_str(), outstr.c_str(), instr == outstr);
  outstr = bv.ToBinaryString();
  printf("Binary %s\n\n", outstr.c_str());
}

int main(int argc, char** argv) {
  DoTest("1");
  DoTest("2");
  DoTest("4");
  DoTest("8");
  DoTest("01");
  DoTest("02");
  DoTest("04");
  DoTest("08");

  DoTest("F0F0F0F0");
  DoTest("DEADBEEF");
  DoTest("D34DB33F");
  DoTest("00000001");
  DoTest("00000003");
  DoTest("000000F0");

  MIBitVector a("010000");
  MIBitVector b("0100");
  MIBitVector c("01001");
  MIBitVector d, e;
  printf("010000 = 0100? %d\n", a.IsEqual(b));
  printf("010000 = 01001? %d\n", a.IsEqual(c));
  printf("empty = empty? %d\n", d.IsEqual(e));


  MIBitVector f("00000000000000000000000000001");
  size_t i = f.GetSize();
  f.Set(i, true);
  f.Set(i+1, true);
  f.Resize(i);
  printf("converted %s\n", f.ToHexString().c_str());



  return 1;
}

#endif
