#ifndef XMLFILES_H
#define XMLFILES_H

#include <vector>
#include <string>

#include <chemlib/chemlib.h>

#include "ATOMLABEL.h"
#include "SURFDOT.h"
#include "Cfiles.h"
#include "Colors.h"

class Annotation;
class CMapHeaderBase;


//@{
// Returns true if the file is in XML format.
//@}
bool IsXML(FILE *fp);

//@{
// An XML file archive object.
//@}
class XMLArchive : public CArchive
{
    std::vector<std::string> tags;
    int inTag;
    std::string f;
    char indent[2048];
    void calc_indent();
public:
    void BeginTag(const char *tag);
    void BeginTag(const char *tag, const char *attr);
    void EndTag();
    XMLArchive(const char *pathname, int m);
    virtual ~XMLArchive();

    void WriteField(const char *tag, int value);
    void WriteField(const char *tag, float value);
    void WriteField(const char *tag, const char *value);
    void WriteField(const char *tag);
    void WriteField(const chemlib::RESIDUE*);
    void WriteField(const chemlib::MIAtom*);
    void WriteHBond(chemlib::Bond&);
    void WriteBond(chemlib::Bond&);
    void WriteConnect(chemlib::Bond&);
    void WriteLabel(const ATOMLABEL&);
    void WriteColor(PaletteColor&);
    void WriteDot(SURFDOT&);
    void WriteSymmCenter(float center[3]);
    void WriteMapHeader(CMapHeaderBase&);
    void WriteAnnotation(Annotation*);
};

#endif // ifndef XMLFILES_H
