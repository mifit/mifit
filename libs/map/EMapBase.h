#ifndef mifit_map_EMapBase_h
#define mifit_map_EMapBase_h

#include <QObject>
#include <string>

#include <math/mathlib.h> // for PLINE

#include "CMapHeaderBase.h"
#include "MapSettingsBase.h"
#include "CREFL.h"
#include "maptypes.h"

namespace chemlib {
class RESIDUE;
class MIAtom;
}

//@{
// Contouring types for the contour map function.
//@}
typedef enum {MAP_CUBE = 0, MAP_SPHERE, MAP_BBLOB, MAP_BOX} CONTOUR_METHOD_TYPE;


struct mmtz_column_;

//@{
// Electron density map class.
// Holds the map data, cell and spacegroup information, as well as methods
// to contour, load and save maps.
// Probably suffering a bit from feature bloat...
//@}
class EMapBase : public QObject {
    Q_OBJECT
  
  /**
   * Copy constructor declared only to prevent access.
   */
  EMapBase(const EMapBase&);
  
  /**
   * Assignment operator declared only to prevent access.
   */
  EMapBase& operator=(const EMapBase&);
  
protected:
  std::string _fostr,_fcstr,_fomstr,_phistr,_sigfstr,_freeRstr;

  //@{
  // hide the map if false
  //@}
  bool visible;
  
  bool modified;
  int contur_sec(int planedirection, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int level, int color, float center[3]);
  long SavemmCIF(FILE* fp);
  long SaveXtalViewPhase(FILE* fp);
  long SaveWarpPhase(FILE* fp);
  long SaveCNSPhase(FILE* fp);
  int sfFFT(chemlib::RESIDUE* res, float& scale);

  std::vector<float> map_points;
  std::vector<float> section;
  std::vector<MAP_POINT> PList;
  float RList[5][10];
  float maxdens(chemlib::MIAtom* CAprev, chemlib::MIAtom* CA1, int l, float fx0, float fy0, float fz0);
  bool FreeRSet;
  bool FOMsValid;
  bool FosValid;
  bool FcsValid;
  bool PhicsValid;

  MapSettingsBase* settings;
  std::string version_string;

  bool GetIndicesFromDefaults(unsigned int num_cols,
                              mmtz_column_* col,
                              int&, int&, int&, int&, int&, int&);

  /**
   * Whether the map is predicted as a difference map.
   * Calculated in CalcRMS.
   */
  bool predictedAsDifferenceMap;
  
public:
  MapSettingsBase *GetSettings() { return settings; }

  void SetVersion(const std::string& version) {
    version_string = version;
  }

  float MapLevel(int index) {
    return settings->MapLevel[index];
  }

  void MapLevel(int index, float level) {
    settings->MapLevel[index] = level;
  }

  int MapLevelColor(int index) {
    return settings->MapLevelColor[index];
  }

  void MapLevelColor(int index, int color) {
    settings->MapLevelColor[index] = color;
  }

  bool MapLevelOn(int index) {
    return settings->MapLevelOn[index];
  }

  void MapLevelOn(int index, bool on) {
    settings->MapLevelOn[index] = on;
  }

  int m_radiusmax() {
    return settings->m_radiusmax;
  }

  void m_radiusmax(int radiusmax) {
    settings->m_radiusmax = radiusmax;
  }

  /**
   * Whether the map is predicted as a difference map.
   * Calculated in CalcRMS which is run when the map is loaded.
   */
  bool isPredictedAsDifferenceMap() {
    return predictedAsDifferenceMap;
  }
  
  //@{
  // the mapheader holds space group, resolution and unit cell information.
  //@}
  CMapHeaderBase* mapheader;
  //@{
  // an STL vector of reflection data.
  //@}
  std::vector<CREFL> refls;
  //@{
  //
  //@}
  float refls_stholmin, refls_stholmax;
  //@{
  // scale is used to scale Fobs to Fcalc.
  //@}
  float scale;
  //@{
  //
  //@}
  float mapmin, mapmax;
  //@{
  // an STL vector of the contoured edges (chickenwire).
  //@}
  std::vector<PLINE> edges;
  //@{
  // if true contour on an orthogonal grid instead of along cell directions
  //@}
  bool orthgrid;
  //@{
  // the spacing of the grid for orthgrid
  //@}
  float gridspacing;
  //@{
  // the current atoms being fit to the map.
  //@}
  std::vector<chemlib::MIAtom*> * CurrentAtoms;
  //@{
  // Use non-crystallographic symmetry.
  // Non NCR is currently implemented, this is for the future.
  //@}
  bool UseNCR;
  //@{
  // a file pointer for logging.
  // used only in debugging mode.
  //@}
  FILE* flog;
  //@{
  // fine, medium or coarse grid.
  // set in the contour dialog box
  // used for FFT'ing a map from structure factors
  //@}
  int mapGrid;
  //@{
  // the name of the map for labeling: can be changed
  //@}
  std::string mapName;
  /**
   * the name of the column used for DirectFFT type maps
   */
  std::string fColumnName;
  //@{
  // the file name for the map; should never change once file is loaded
  //@}
  std::string pathName; 
  //@{
  // true if the map has changed
  //@}
  bool changed;
  //@{
  //
  //@}
  int mapnumber;
  //@{
  // returns the index of a map point given the 3-space coordinates
  //@}
  int mapdex(int ix, int iy, int iz) {
    return (mapheader->nx*(mapheader->ny*((iz+64*mapheader->nz)%mapheader->nz)+
                           (iy+64*mapheader->ny)%mapheader->ny)+(ix+64*mapheader->nx)%mapheader->nx);
  }

  bool IsModified() {
    return modified;
  }

  void SetModified(bool m = true) {
    modified = m;
  }

  long Reindex(int index_mat[3][3]);
  float CorrScore(std::vector<chemlib::MIAtom*> atoms);
  float RFactor;
  const std::vector<CREFL>& GetRefls() {
    return refls;
  }

  long AddFreeRFlag(int percent, bool use_shells);

  bool HasFreeRFlag() { return FreeRSet;  }
  bool HasFOMs() { return FOMsValid; }
  bool HasFos() { return FosValid; }
  bool HasFcs() { return FcsValid; }
  bool HasPhic() { return PhicsValid; }

  void RecalcResolution();
  float Rho(float x, float y, float z);

  static bool IsCif(const char* pathname);
  static bool IsCCP4MTZFile(const char* pathname);
  static bool IsRefFile(const char* pathname);
  static bool IsFinFile(const char* pathname);
  static bool IsScaFile(const char* pathname);
  static bool IsPhsFile(const char* pathname);

  static bool IsCCP4MapFile(const char* pathname);
  static bool IsCNSMapFile(const char* pathname);



  //@{
  // interpret a script command.
  // if it is an Emap command and it has been processed this returns true.
  //@}
  bool ScriptCommand(const char* buf);
  //@{
  // the name of the crystal.
  //@}
  std::string Crystal;
  //@{
  // return a string describing this EMap.
  //@}
  std::string Info();
  //@{
  //
  //@}
  bool ScaleSolventPercent(float percentSolvent);
  //@{
  // Smooth the map by radius r.
  //@}
  bool SmoothMap(float r);
  //@{
  // Calulate a solvent mask given a percent solvent and a radius.
  //@}
  bool CalcSolventMask(float percentSolvent, float radius);
  //@{
  // Copy the current atoms into the EMap container CurrentAtoms.
  //@}
  void SetCurrentAtoms(std::vector<chemlib::MIAtom*> * current);
  //@{
  // Finds the best 6 positions to put the next Calpha carbon given
  // the current and the previous Calpha atoms.
  //@}
  std::vector<MAP_POINT>& PlaceCA(chemlib::MIAtom* CAprev, chemlib::MIAtom* CA1);
  //@{
  // Score the fit of the atoms with the map density.
  // The score is not equivalent to an r-value despite the name of the function.
  // The higher the score, the better the fit.
  //@}
  float RDensity(std::vector<chemlib::MIAtom*>& atoms);
  //@{
  // Score the fit of the atoms with the map density with a correlation function.
  // The higher the score, the better the fit.
  //@}
  float RCorrelation(std::vector<chemlib::MIAtom*>& atoms);
  //@{
  // the center of the map in x, y, z Cartesian coords.
  //@}
  float map_center[3];
  //@{
  // check if the center has moved for recontouring purposes.
  //@}
  void CheckCenter(float x, float y, float z);
  //@{
  // self-explanatory
  //@}
  void GetCrystalFromUser();
  //@{
  // Read a map from a file.
  // @param the name of the file.
  //@}
  int Read(const char* file);
  //@{
  // Read a crystal from a file.
  // @param the name of the file.
  //@}
  //@{
  // return the value of the boolean changed.
  //@}
  bool GetChanged() {
    return changed;
  }

  //@{
  // set teh value of the boolean changed.
  //@}
  void SetChanged(bool c) {
    changed = c;
  }

  //@{
  // return the map number.
  //@}
  int GetMapNumber() {
    return mapnumber;
  }

  //@{
  // set the map number.
  //@}
  void SetMapNumber(int n) {
    mapnumber = n;
  }

  //@{
  // return the map line width.
  //@}
  float GetMapLinewidth() {
    return settings->maplinewidth;
  }

  //@{
  // set the map line width.
  //@}
  void SetMapLinewidth(float n) {
    settings->maplinewidth = n;
  }

  //@{
  // show the map.
  //@}
  void Show() {
    visible = true;
    mapVisibilityChanged(this);
  }

  //@{
  // hide the map.
  //@}
  void Hide() {
    visible = false;
    mapVisibilityChanged(this);
  }

  //@{
  // is the map visible?
  //@}
  bool Visible() {
    return visible;
  }

  //@{
  // calculate structure factors from a model.
  // Crystal must be set properly for this work.
  // May also want to set the resolution limits.
  //@}
  bool SFCalc(chemlib::RESIDUE* res);
  //@{
  // Load a CCP4 mtz map (reflections file).
  // the function guesses the columns that will be used and asks the
  // user to verify.
  // @param pathname the name of the file.
  //@}
  long LoadMTZMap(const char* pathname, bool require_fo=true);
  //@{
  // Save relection data in mtz format.
  // @param pathname the name of the file to be created.
  //@}
  long SaveCCP4Phase(const char* pathname);
  //@{
  // set the crystal.
  // Will load the information from the crystal database entry of the same name.
  //@}
  void SetCrystal(const char* crystal);
  //@{
  // enumeration of the type of data the map was created from.
  // Important for getting a number of functions right.
  // Also used for the save functions.
  //@}
  enum
  {
    mmCIF_phase,
    XtalView_phase,
    Warp_phase,
    CCP4_phase,
    CNS_phase,
    XtalView_map,
    CCP4_map,
    XtalView_fin,
    HKL_scaI,
    HKL_scaF,
    DTrek_ref
  };

  //@{
  // Save the structure factors in the EMap to a file.
  // @param pathname name of the file
  // @param type the format of the file to be formed.
  //@}
  bool SavePhases(const char* pathname, int type);
  //@{
  // Load an XtalView/FSFOUR format density map file.
  //@}
  long LoadFSFOURMapFile(const char* pathname);
  //@{
  // Load the structure factors from a file.
  //@}
  bool LoadMapPhaseFile(const char* s, bool require_fo=true);

  //@{
  // Interpret any column label strings (FO=foo FC=foo, etc) on the given line
  //@}
  bool InterpretColumnLabelString(const char *s);

  //@{
  // loads a map from a density file specified in s.
  // Attempts to determine the map format before loading.
  //@}
  bool LoadMapFile(const char* s);

  //@{
  // load a CCP4 density map.
  // Will expand the symmetry.
  //@}
  long LoadCCP4Map(const char* pathname, float* rms, float* min, float* max);
  //@{
  // load a CNS density map.
  // Will expand the symmetry.
  //@}
  long LoadCNSMap(const char* pathname);
  //@{
  // A string containing the map id, the mapName, type and resolution paramters
  //@}
  std::string MapID();
  //@{
  //
  //@}
  EMapBase();
  //@{
  //
  //@}
  virtual ~EMapBase();
  //  bool ReSize(int, int, int);
  //@{
  // return a pointer to the map header
  //@}
  CMapHeaderBase* GetMapHeader() {
    return mapheader;
  }

  //@{
  // load the map from a phase file
  //@}
  virtual long LoadMap(const char* pathname, int type = XtalView_phase);
  //@{
  // load the map from a CIF format file.
  //@}
  virtual long LoadCIFMap(const char* pathname, int datablock = 0);

  // return maptypes which could possibly be made from the given file, and maptypes possible if Fc is calculated
  static bool GetPossibleMapTypes(const std::string& pathname, 
                                  std::vector<unsigned int>& types, 
                                  std::vector<unsigned int>& with_fc_types);

  // return maptypes which can be made given the currently-defined fields, and assmung Fc is calculated in with_fc_types
  bool GetMapTypes(std::vector<unsigned int>& types, 
                   std::vector<unsigned int>& with_fc_types);

  // with the current state, can we produce a map of the requested type?
  bool CanDoMapType(unsigned int type);

  //@{
  // FFT the map.
  // Gets information from the user with a dialog box then calls FFTCalc.
  //@}
  virtual bool FFTMap(int maptype = -1, int gridlevel = -1,
                      float resMin = -1.0f, float resMax = -1.0f);

  //@{
  // Calculate the structure factors.
  // Call if the refection list changes or the resolution changes.
  //@}
  bool FFTCalc();
  //@{
  // Scale the map so that one rms/sigma is 50.0.
  //@}
  void ScaleMap(float rms);
  float CalcRMS();
  //@{
  // used in solvent flattening to smooth the map.
  //@}
  bool SigmaMap();
  //@{
  // contour within a box.
  //@}
  void ContourBox(float x1, float x2, float y1, float y2, float z1, float z2);
  //@{
  // non-crystallographically symmetry average the map.
  // There is no user interface yet.
  //@}
  int NCRAverage(int ix, int iy, int iz);
  //@{
  // calculate a line segment of the chicken-wire map.
  //@}
  int segment(double* x1, double* y1, double* x2, double* y2, int i, int j,
              float* sec, int nrow, int upper, int level);
  //@{
  // Contour the map given a center and current atoms (for blob).
  //@}
  void Contour(float center[3], std::vector<chemlib::MIAtom*> * current);
  //@{
  // interpolate the value of rho given a point in fractional coords.
  //@}
  float avgrho(float fx, float fy, float fz);
  //@{
  // return true if there are phases/structure factors.
  //@}
  bool HasPhases() {
    return (refls.size() > 0);
  }

  //@{
  // return true if the map has density points.
  //@}
  bool HasDensity() {
    return map_points.size() > 8;
  }

  //@{
  // sprinkle water throughout the map into empty peaks that make sense to be water.
  //@}
  int HydrateMap(int minlevel, int maxadd, int add_water, chemlib::MIMoleculeBase*, float dmin, float dmax,
                 float xmin = 0.0, float xmax = 1.0F, float ymin = 0.0, float ymax = 1.0F, float zmin = 0.0, float zmax = 1.0F);

  // only meaningful in interactive version
  virtual bool PromptForCrystal() {
    return false;
  }

  // Get the column labels from the given file
  static bool GetMTZColumnInfo(const std::string& fname, 
                               std::vector<std::string> &labels,
                               std::vector<char> &types);

  // returns true if the given file is a mtz file and has all of the column labels in the vector
  bool HasColumnLabels(const std::string& fname, const std::vector<std::string> &labels);

  // use the given values for defaults for column labels
  void UseColumnLabels(const std::string& fostr, const std::string& fcstr,
                       const std::string& fomstr, const std::string& phistr,
                       const std::string& sigfstr, const std::string& freeRstr);

  virtual bool PromptForColumnLabels(unsigned int /* num_cols */,
                                     mmtz_column_* /* col */,
                                     int&, int&, int&, int&, int&, int&) {
    return false;
  }

Q_SIGNALS:
  void mapContourLevelsChanged(EMapBase*);
  void mapVisibilityChanged(EMapBase*);

};

// The following is declared only for testing purposes
bool parseCoefficients(const char* str, int& mapType);

#endif
