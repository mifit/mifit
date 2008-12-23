#ifndef MIMolIO_H
#define MIMolIO_H

#include <string>

#include "MIMolInfo.h"

namespace chemlib {
class MIAtom;

class MIMolIOBase {
public:
  MIMolIOBase();
  virtual ~MIMolIOBase();

  virtual bool Write(MIMolInfo& mol, const std::string& filename, int writerIndex = -1) const;
  virtual bool Read(MIMolInfo& mol, const std::string& filename, int readerIndex = -1) const;

  void registerReader(Reader* reader);
  void registerWriter(Writer* writer);

  // Returns the index of a Reader/Writer which will read/write the format specified by the extension.
  int getReaderIndex(const std::string& filename) const;
  int getWriterIndex(const std::string& filename) const;

protected:
  unsigned int GetReaderCount() const {
    return readers.size();
  }

  unsigned int GetWriterCount() const {
    return writers.size();
  }

  std::string GetReaderDescription(unsigned int idx) const;
  std::string GetWriterDescription(unsigned int idx) const;

  typedef std::vector<Reader*> ReaderList;
  typedef std::vector<Writer*> WriterList;

  ReaderList readers;
  WriterList writers;
};

class MIColorSetter {
public:
  virtual ~MIColorSetter() {
  }

  virtual bool operator()(MIAtom*) const = 0;
  virtual bool operator()(MIAtom*, char) const = 0;
};

void MIRegisterColorSetter(MIColorSetter*);
MIColorSetter* MIGetColorSetter();

}

#endif // MIMolIO_H
