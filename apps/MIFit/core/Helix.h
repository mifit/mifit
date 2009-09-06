#ifndef mifit_model_Helix_h
#define mifit_model_Helix_h

namespace chemlib{
class RESIDUE;
}

class Helix {
public:
  Helix(double radius);
  virtual ~Helix();

  //For list operations
  Helix* m_pNext;

  bool MakeHelix(chemlib::MIMoleculeBase* mol,
      chemlib::RESIDUE* pHelixStart, chemlib::RESIDUE* pHelixStop);
  bool CreateAxis(int axis_len, double caA[4][3], double caB[4][3]);
  bool ProjectPoint(double pt[3], double enda[3], double endb[3], double ptp[3]);
  bool CreateCylinder();
  void SetColor(unsigned char red, unsigned char green, unsigned char blue);

private:

  friend class GLRenderer;

#define Helix_NSEG 8
  double m_dRadius;

  //Ends of the radius vector
  double m_dAxisA[3];
  double m_dAxisB[3];

  //Cylinder end points and normals
  double m_dEndA[Helix_NSEG][3];
  double m_dEndB[Helix_NSEG][3];
  double m_dNormA[Helix_NSEG][3];
  double m_dNormB[Helix_NSEG][3];

  //Cap Normals
  double m_dCapNormA[Helix_NSEG][3];
  double m_dCapNormB[Helix_NSEG][3];

  unsigned char rgb[4];

};

#endif
