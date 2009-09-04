#ifndef mifit_map_CMapHeaderBase_h
#define mifit_map_CMapHeaderBase_h

#include <string>
#include <vector>
#include <boost/signal.hpp>

#include <chemlib/chemlib.h>

#include "MAP_POINT.h"

namespace MISymmop {
  const unsigned int MAXSTRING = 1024;
  const unsigned int MAXSYMMOPS = 192;
  const unsigned int MAXNCRSYMMOPS = 48;
}

//@{
// class encapsulating map header information such as unit cell, spacegroup.
//@}
struct umtzfile_;

class CMapHeaderBase {
public:
  std::string crystal_name;
  float a, b, c, alpha, beta, gamma;
  int spgpno;
  std::string spgpname;
  std::string SymInfoString;
  std::string title;
  std::vector<std::string>SymopsString;
  std::string PG_symbol;
  int npg_ops;
  float NCRSymmops[MISymmop::MAXNCRSYMMOPS][12];
  int nNCRSymmops;
  int nx, ny, nz;
  double scale;   /* scale on density */
  float Bsolvent, Ksolvent;  /* bulk solvent model */
  float sc11, sc22, sc33, sc12, sc13, sc23;   /* anisotropic scaling factors */
  int use_bulksolvent, use_aniso;
  float resmax, resmin;
  int nsym;
  float symops[3][4][MISymmop::MAXSYMMOPS];
  float ctof[3][3];   /* cartesian to fractional */
  float ftoc[3][3];   /* fractional to cartesian */
  int maptype;   /* fft.h contains maptypes info */
  int fc_is_fom;
  int hmin, hmax, kmin, kmax, lmin, lmax;
  float M_Coefficient;
  float N_Coefficient;

  CMapHeaderBase();
  CMapHeaderBase(const CMapHeaderBase& rhs);
  CMapHeaderBase& operator=(const CMapHeaderBase& rhs);
  CMapHeaderBase(const std::string &from_label);
  virtual ~CMapHeaderBase() {
  }

  void updateSymmetryAndCell(const CMapHeaderBase& mapHeader); // sets only spacegroup / cell

  const char* GetHMSymbol();   // international Hermann-Mauguin space group symbol
  int FindHMSymbol(const char*);
  bool AsuBounds(float& xmin, float& xmax, float& ymin, float& ymax, float& zmin, float& zmax);
  bool IsOk();
  int scan_symmops(const char* buf, float symop[3][4][MISymmop::MAXSYMMOPS], char opstring[MISymmop::MAXSYMMOPS][MISymmop::MAXSTRING]);
  int addcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addinversion(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addIcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addFcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addBcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addCcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addAcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int addRcenter(int n, float symop[3][4][MISymmop::MAXSYMMOPS]);
  int SpgpFromSymmops();
  std::string Label() const;
  void getfx(chemlib::MIAtom* CA1, float* fx0, float* fy0, float* fz0);
  void getfx(MAP_POINT* CA1, float* fx0, float* fy0, float* fz0);

  virtual bool LoadCrystal(const char*) {
    return false;
  }

  virtual bool SaveCrystal(const std::string&) {
    return false;
  }

  bool GetSymInfo();
  bool WriteSymmOpsMTZ(umtzfile_*);
  void ScanSymops();

  /* non-crystallographic symmops */
  void symm_mh(float x, float y, float z, float* xp, float* yp, float* zp, int jth);
  void unsymm_mh(float x, float y, float z, float* xp, float* yp, float* zp, int jth);
  int FindSpacegroup(const char* symbol);
  void SetSymmOps();
  void NCRTransform(float* fx, float* fy, float* fz, int jth);
  void CtoF(float* x, float* y, float* z) const {
    float xp = ctof[0][0]*(*x)+ctof[0][1]*(*y)+ctof[0][2]*(*z);
    float yp = ctof[1][0]*(*x)+ctof[1][1]*(*y)+ctof[1][2]*(*z);
    float zp = ctof[2][0]*(*x)+ctof[2][1]*(*y)+ctof[2][2]*(*z);
    *x = xp;
    *y = yp;
    *z = zp;
  }

  void FtoC(float* x, float* y, float* z) const {
    float xp = ftoc[0][0]*(*x)+ftoc[0][1]*(*y)+ftoc[0][2]*(*z);
    float yp = ftoc[1][0]*(*x)+ftoc[1][1]*(*y)+ftoc[1][2]*(*z);
    float zp = ftoc[2][0]*(*x)+ftoc[2][1]*(*y)+ftoc[2][2]*(*z);
    *x = xp;
    *y = yp;
    *z = zp;
  }

  void EchoCrystal();

  boost::signal1<void, CMapHeaderBase*> mapHeaderChanged;

};

bool MapHeaderOK(CMapHeaderBase*);
void MIGetSpacegroups(std::vector<std::string> &spacegropus);
#endif
