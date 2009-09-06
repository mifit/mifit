#ifndef MI_BITVECTOR_H
#define MI_BITVECTOR_H

#include <string>
#include <vector>

class MIBitVector {
public:
  MIBitVector(const std::string& hex_data = "");

  bool FromHexString(const std::string& hex_data);
  std::string ToHexString() const;
  std::string ToBinaryString() const;

  void Set(size_t i, bool state);
  bool IsSet(size_t i) const;

  size_t GetSize() const;
  void Resize(size_t size);

  size_t LastBitSet() const;
  bool IsEqual(const MIBitVector& other) const;

private:
  std::vector<unsigned char> _data;
};

#endif
