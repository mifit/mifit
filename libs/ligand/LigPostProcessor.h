#ifndef mifit_ligand_LigPostProcessor_h
#define mifit_ligand_LigPostProcessor_h

#include "ligandlib.h"

class LigPostProcessor {

  LigDictEntry& m_mon;

  bool genCoords;
  bool guessBondOrders;
  bool genConstraints;

  void GenerateConstraints(bool allConstraints);

public:
  LigPostProcessor(LigDictEntry& mon, const char* format);

  void Process();
  void GenerateCoordinates();
  void GuessBondOrders();
  void GenerateConstraints();

  void setDoGenerateCoordinates(bool on);
  void setDoGuessBondOrders(bool on);
  void setDoGenerateConstraints(bool on);
};

#endif

