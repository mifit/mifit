#ifndef __LSQFitDialog_H__
#define __LSQFitDialog_H__

#include <vector>

#include "chemlib.h"
#include "corelib.h"

class Molecule;
class Displaylist;

#include "ui_LSQFitDialog.h"

//@{
// Struct for cantaining a match between two chains.
//  Used by LSQFit dialog.
//@}
class MATCH {
public:
  chemlib::RESIDUE* target;
  chemlib::RESIDUE* source;
  int length;
  std::string atomtype;
};

//@{
// Dialog box for least-squares overlapping of two chains.
// Allows for multiple segementsw to be matched.
// Segments can either be in the same molecule (For example
// A and B chains, or in two different realated molecules.
//  Someday we should make this able to automatically find the overlaps
// say with A GA optimization of a distance matrix...
//@}
class LSQFitDialog : public QDialog, public Ui::LSQFitDialog {
    Q_OBJECT

public:
  LSQFitDialog(QWidget *parent);

  void InitializeFromData(const MIData &dat);
  void GetData(MIData &dat);

  void SetAtomType(const char*);
  void ListTargetChoices();
  void ListSourceChoices();
  Displaylist* displaylist;
  Molecule* m_source;
  Molecule* m_target;
  bool success;
  double r[3][3];
  double v[3];
  double rms;
  void ListTarget();
  void ListSource();
  void ListMatches();
  bool MatchOK(MATCH*);
  chemlib::RESIDUE* targetres, * sourceres;

  std::vector<MATCH> Matches;

private Q_SLOTS:
  void on_targetComboBox_currentIndexChanged(const QString &);
  void on_sourceComboBox_currentIndexChanged(const QString &);
  void on_sourceListWidget_currentTextChanged(const QString &str);
  void on_targetListWidget_currentTextChanged(const QString &str);
  void on_calcButton_clicked();
  void on_removeButton_clicked();
  void on_addButton_clicked();
  void on_exportButton_clicked();
  void on_importButton_clicked();

private:
  void updateButtons();
  void updateMatrix();
  void setSuccess(bool value);
  Molecule *findMolecule(const std::string &);
};

#endif
