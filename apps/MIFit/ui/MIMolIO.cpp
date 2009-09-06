#include "MIMolIO.h"
#include "macafxwin.h"
#include "uitest.h"
#include "ui/MIDialog.h"

using namespace chemlib;

MIMolIO::MIMolIO() : MIMolIOBase() {
}

bool MIMolIO::Write(MIMolInfo& mi, const std::string& fname, int writerIndex) const {
  std::string filename = fname;
  if (filename.size() == 0) {
    std::string filter;
    WriterList::const_iterator iter = writers.begin();
    for (iter = writers.begin(); iter != writers.end(); ++iter) {
      Writer* writer = *iter;
      filter += writer->getDescription();
      filter += "|";
      filter += writer->getExtension();
      filter += "|";
    }
    filter += "All files (*.*)|*.*";
    filename = MIFileSelector("Choose a file to save to",
                              "", "", "", filter.c_str(), MI_SAVE_MODE);
  }
  if (filename.size() == 0) {
    return false;
  }

  if (writerIndex == -1) {
    writerIndex = getWriterIndex(filename);
  }
  if (writerIndex < 0) {
    std::vector<std::string> writerDescriptions;
    for (unsigned int i = 0; i < GetWriterCount(); ++i) {
      writerDescriptions.push_back(writers[i]->getDescription().c_str());
    }
    writerIndex = MIGetSingleChoiceIndex("Select Output Format", "Output format",
                    writerDescriptions);
  }
  if (writerIndex < 0) {
    return false;
  }

  return MIMolIOBase::Write(mi, filename, writerIndex);
}

bool MIMolIO::Read(MIMolInfo& mi, const std::string& fname, int readerIndex) const {
  std::string filename = fname;
  if (filename.size() == 0) {
    std::string filter;
    ReaderList::const_iterator iter = readers.begin();
    for (iter = readers.begin(); iter != readers.end(); ++iter) {
      Reader* reader = *iter;
      filter += reader->getDescription();
      filter += "|";
      filter += reader->getExtension();
      filter += "|";
    }
    filter += "All files (*.*)|*.*";
    filename = MIFileSelector("Choose a file to read",
                              "", "", "", filter.c_str(), MI_OPEN_MODE);
  }
  if (filename.size() == 0) {
    return false;
  }

  if (readerIndex == -1) {
    readerIndex = getReaderIndex(filename);
  }
  if (readerIndex < 0) {
    std::vector<std::string> readerDescriptions;
    for (unsigned int i = 0; i < GetReaderCount(); ++i) {
      readerDescriptions.push_back(readers[i]->getDescription().c_str());
    }
    readerIndex = MIGetSingleChoiceIndex("Select Input Format", "Input format",
                    readerDescriptions);
  }
  if (readerIndex < 0) {
    return false;
  }

  return MIMolIOBase::Read(mi, filename, readerIndex);
}

