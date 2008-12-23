#ifndef MIFIT_MOLECULEXMLHANDLER_H_
#define MIFIT_MOLECULEXMLHANDLER_H_

#include <QXmlDefaultHandler>
#include <vector>

#include "Xguicryst.h"
#include "Annotation.h"
#include "CMapHeader.h"


class Molecule;

/**
 * XML molecular SAX parser for .mlw format.
 */
class MoleculeXmlHandler : public QXmlDefaultHandler {

  char tmpElement[20];
  chemlib::MIAtom* atom;
  chemlib::RESIDUE* residueList;
  chemlib::RESIDUE* residue;
  chemlib::RESIDUE* prev_res;
  Molecule* molecule;
  std::vector<chemlib::MIAtom*> atoms;
  chemlib::Bond bond;
  SURFDOT dot;
  ATOMLABEL* label;
  bool in_MapHeader;
  bool in_Annotation;
  bool in_Defaults;
  bool in_Molecule;
  bool in_EMap;
  Annotation* ann;
  std::string data;
  CMapHeader mapheader;

  QString errorString_;

public:
  MoleculeXmlHandler();
  ~MoleculeXmlHandler();


  virtual bool characters(const QString& ch);
  virtual bool startElement(const QString& namespaceURI, const QString& localName, const QString& qName, const QXmlAttributes& atts);
  virtual bool endElement(const QString& namespaceURI, const QString& localName, const QString& qName);

  virtual QString errorString() const {
    return errorString_;
  }

  virtual bool error(const QXmlParseException& exception);
  virtual bool fatalError(const QXmlParseException& exception);
  virtual bool warning(const QXmlParseException& exception);

  typedef std::vector<Molecule*> MoleculeList;
  MoleculeList molecules;

  typedef std::vector<CMapHeaderBase> MapHeaderList;

  /**
   * These map headers are read from XML session files, but since
   * the maps are not loaded when the XML is parsed, they must be
   * stored for later use when the maps are loaded.
   */
  MapHeaderList mapHeaders;
};

#endif
