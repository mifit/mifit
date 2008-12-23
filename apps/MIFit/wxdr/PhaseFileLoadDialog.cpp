#include "PhaseFileLoadDialog.h"
#include "maplib.h"
#include "MIDialog.h"

static const char* F_TEXT = "F = ";
static const char* FO_TEXT = "Fo = ";
#define MODEL_LABEL "Model "

static bool isModelString(const char* str) {
  return startsWith(str, MODEL_LABEL);
}

static int modelStringToNumber(const char* str) {
  int result = -1;
  if (isModelString(str)) {
    result = atoi(&str[strlen(MODEL_LABEL)]) - 1;
  }
  return result;
}

static std::string numberToModelString(unsigned int n) {
  return ::format(MODEL_LABEL "%u", n+1);
}

static void SetChoiceOptions(QComboBox *choice,
                  const std::vector<std::string>& file_fcs,
                  const std::vector<std::string>& models) {
  choice->clear();
  for (size_t i=0; i < file_fcs.size(); ++i) {
    choice->addItem(file_fcs[i].c_str());
  }
  for (size_t i=0; i < models.size(); ++i) {
    choice->addItem(numberToModelString(i).c_str());
  }
  choice->setCurrentIndex(0);
}



PhaseFileLoadDialog::PhaseFileLoadDialog(QWidget *parent) 
  : QDialog(parent) {
  _res_tmp_map=0;
  resmin_userset=false;
  resmax_userset=false;
  setupUi(this);
}

PhaseFileLoadDialog::~PhaseFileLoadDialog() { 
  delete _res_tmp_map;
}

void PhaseFileLoadDialog::ShowExtras(unsigned int mapNum, bool show) {
  if (mapNum == 1) {
    fc1Label->setVisible(show);
    map1FcChoice->setVisible(show);
    fom1Label->setVisible(_showFom && show);
    map1FomChoice->setVisible(_showFom && show);
  } else {
    fc2Label->setVisible(show);
    map2FcChoice->setVisible(show);
    fom2Label->setVisible(_showFom && show);
    map2FomChoice->setVisible(_showFom && show);
  }
}

void PhaseFileLoadDialog::updateCheckBoxes(unsigned int num) {
  if (num == 1) {
    map1Check->setChecked(true);
  } else {
    map2Check->setChecked(true);
  }
}

void PhaseFileLoadDialog::OnMapTypeChange(unsigned int num)
{
  mapType_userset[num-1] = true;
  fo_userset[num-1] = false;
  fc_userset[num-1] = false;
  phi_userset[num-1] = false;
  fom_userset[num-1] = false;
  updateCheckBoxes(num);
  updateDefaults(num);
}

void PhaseFileLoadDialog::OnMapFoChange(unsigned int num)
{
  fo_userset[num-1] = true;
  updateCheckBoxes(num);
  updateDefaults(num);
}

static bool phiLastSet = false;

void PhaseFileLoadDialog::OnMapFcChange(unsigned int num)
{
  fc_userset[num-1] = true;
  phiLastSet = false;
  updateCheckBoxes(num);
  updateDefaults(num);
  if (num==1)
    UpdateResolution();
}

void PhaseFileLoadDialog::OnMapPhiChange(unsigned int num)
{
  phi_userset[num-1] = true;
  phiLastSet = true;
  updateCheckBoxes(num);
  updateDefaults(num);
  if (num==1)
    UpdateResolution();
}

void PhaseFileLoadDialog::OnMapFomChange(unsigned int num)
{
  fom_userset[num-1] = true;
  updateCheckBoxes(num);
  updateDefaults(num);
}

static std::string findPreferredFcItem(const std::string& foStr, const std::string& phiStr,
    const std::vector<std::string>& fcItems, const std::vector<std::string>& models) {
  
  std::string result;
  std::vector<std::string> preferredFcs;
  preferredFcs.push_back(foStr + "C");
  if (isModelString(phiStr.c_str())) {
    preferredFcs.push_back(phiStr);
  }
  preferredFcs.push_back(std::string("FC"));
  std::vector<std::string>::iterator iter;
  for (iter = preferredFcs.begin(); iter != preferredFcs.end(); ++iter) {
    std::string& searchStr = *iter;
    if (std::find(fcItems.begin(), fcItems.end(), searchStr) != fcItems.end()) {
      result = searchStr;
      break;
    }
  }
  // if preferred setting not found and models available, set to first model,
  // otherwise set to first item
  if (result.empty()) {
    if (models.size() > 0) {
      result = numberToModelString(0).c_str();
    } else if (fcItems.size() > 0) {
      result = fcItems[0];
      // if foStr and fcStr are the same, set to second item
      if (foStr == result && fcItems.size() > 1) {
        result = fcItems[1];        
      }
    }
  }
  return result;
}

static std::string findPreferredPhiItem(const std::string& foStr, const std::string& fcStr,
    const std::vector<std::string>& phiItems, const std::vector<std::string>& models) {
  
  std::string result;
  std::vector<std::string> preferredPhis;
  if (!foStr.empty()) {
    preferredPhis.push_back("PH" + foStr);
    if (foStr[0]=='F') {
      preferredPhis.push_back("PH" + std::string(&foStr[1]));
    }
  }
  if (isModelString(fcStr.c_str())) {
    preferredPhis.push_back(fcStr);
  }
  preferredPhis.push_back(std::string("PHI"));
  preferredPhis.push_back(std::string("PHIC"));
  std::vector<std::string>::iterator iter;
  for (iter = preferredPhis.begin(); iter != preferredPhis.end(); ++iter) {
    std::string& searchStr = *iter;
    if (std::find(phiItems.begin(), phiItems.end(), searchStr) != phiItems.end()) {
      result = searchStr;
      break;
    }
  }
  // if preferred setting not found and models available, set to first model,
  // otherwise set to first item
  if (result.empty()) {
    if (models.size() > 0) {
      result = numberToModelString(0).c_str();
    } else if (phiItems.size() > 0) {
      result = phiItems[0];
    }
  }
  return result;
}


static void setStringSelection(QComboBox *cb, const std::string &str){
  for ( int i=0;i<cb->count(); ++i) {
    if (cb->itemText(i) == QString(str.c_str())) {
      cb->setCurrentIndex(i);
      break;
    }
  }
}

void PhaseFileLoadDialog::updateDefaults(unsigned int mapNum) {
  
  // Strings for field values; an empty value will indicate not yet set
  std::string mapTypeStr;
  std::string foStr;
  std::string fcStr;
  std::string phiStr;
  std::string fomStr;

  QComboBox* mapTypeChoice=map1Type;
  QComboBox* foChoice=map1FoChoice;
  QComboBox* fcChoice=map1FcChoice;
  QComboBox* phiChoice=map1PhiChoice;
  QComboBox* fomChoice=map1FomChoice;
  QCheckBox* mapCheck=map1Check;
  QLabel* modelText = map1Model;

  if (mapNum==2) {
    mapTypeChoice=map2Type;
    foChoice=map2FoChoice;
    fcChoice=map2FcChoice;
    phiChoice=map2PhiChoice;
    fomChoice=map2FomChoice;
    mapCheck=map2Check;
    modelText=map2Model;
  }

  // If set by user, get value
  if (mapType_userset[mapNum-1]) {
    mapTypeStr = mapTypeChoice->currentText().toStdString();
  }
  if (fo_userset[mapNum-1]) {
    foStr = foChoice->currentText().toStdString();
  }
  if (fc_userset[mapNum-1]) {
    fcStr = fcChoice->currentText().toStdString();
  }
  if (phi_userset[mapNum-1]) {
    phiStr = phiChoice->currentText().toStdString();
  }
  if (fom_userset[mapNum-1]) {
    fomStr = fomChoice->currentText().toStdString();
  }
  
  // Collect choice items into local vectors
  std::vector<std::string> foItems;
  for (int i=0; i < foChoice->count(); ++i) {
    foItems.push_back(foChoice->itemText(i).toStdString());
  }
  std::vector<std::string> fcItems;
  for (int i=0; i < fcChoice->count(); ++i) {
    fcItems.push_back(fcChoice->itemText(i).toStdString());
  }
  std::vector<std::string> phiItems;
  for (int i=0; i < phiChoice->count(); ++i) {
    phiItems.push_back(phiChoice->itemText(i).toStdString());
  }

  // If only one choice, simply set
  if (foStr.empty() && foItems.size() == 1) {
    foStr = foItems[0];
  }
  if (fcStr.empty() && fcItems.size() == 1) {
    fcStr = fcItems[0];
  }
  if (phiStr.empty() && phiItems.size() == 1) {
    phiStr = phiItems[0];
  }
  
  if (foStr.empty() && phiStr.empty() && fcStr.empty()
      && (mapTypeStr.empty() || mapTypeStr == "Direct FFT")) {
    // Search for precomputed structure factors. If found, 
    // set foStr, phiStr, and mapTypeStr
    std::vector<std::pair<std::string, std::string> > recognizedPrecomps;
    if (mapNum == 1) {
      recognizedPrecomps.push_back(std::make_pair("FWT", "PHWT"));
      recognizedPrecomps.push_back(std::make_pair("2FOFCWT", "PH2FOFCWT"));
      // For Buseter "SigmaA" weighted 2Fo-Fc Electron Density Map
      recognizedPrecomps.push_back(std::make_pair("2FOFCWT", "2PHFOFCWT"));
      // Possible typo in Buster manual
      recognizedPrecomps.push_back(std::make_pair("2FOFWT", "2PHFOFCWT"));
    } else {
      recognizedPrecomps.push_back(std::make_pair("DELFWT", "PHDELFWT"));
      // For Buster "SigmaA" weighted  Fo-Fc Difference Density Map
      recognizedPrecomps.push_back(std::make_pair("FOFCWT", "PHFOFCWT"));
      // "SigmaA" weighted  Fo-Ffrag-Fsolv Difference Density Map
      recognizedPrecomps.push_back(std::make_pair("FOFRGSLV", "PHFOFRGSLV"));
    }
    std::vector<std::pair<std::string, std::string> >::iterator iter;
    for (iter = recognizedPrecomps.begin(); iter != recognizedPrecomps.end(); ++iter) {
      std::string foSearch = iter->first;
      std::string phiSearch = iter->second;
      if (std::find(foItems.begin(), foItems.end(), foSearch) != foItems.end() &&
          std::find(phiItems.begin(), phiItems.end(), phiSearch) != phiItems.end()) {
        foStr = foSearch;
        phiStr = phiSearch;
        if (mapTypeStr.empty()) {
          mapTypeStr = "Direct FFT";
        }
        mapCheck->setChecked(true);
        break;
      }
    }
  }
  if (foStr.empty()) {
    std::vector<std::string> preferredFos;
    if (phiStr.size() && phiStr[0]=='P' && phiStr[1]=='H') {
      preferredFos.push_back(&phiStr[2]);
      preferredFos.push_back("F" + std::string(&phiStr[2]));
    }
    preferredFos.push_back(std::string("FO"));
    preferredFos.push_back(std::string("FNAT"));
    preferredFos.push_back(std::string("FP"));
    preferredFos.push_back(std::string("F"));
    std::vector<std::string>::iterator iter;
    for (iter = preferredFos.begin(); iter != preferredFos.end(); ++iter) {
      std::string& searchStr = *iter;
      if (std::find(foItems.begin(), foItems.end(), searchStr) != foItems.end()) {
        foStr = searchStr;
        break;
      }
    }
    // if preferred setting not found, set to first item
    if (foStr.empty() && foItems.size() > 0) {
      foStr = foItems[0];
    }
  }
  if (phiStr.empty()) {
    phiStr = findPreferredPhiItem(foStr, fcStr, phiItems, _models);
  }
  if (fcStr.empty()) {
    fcStr = findPreferredFcItem(foStr, phiStr, fcItems, _models);
  }

  if (mapTypeStr.empty()) {
    if (mapNum == 1) {
      mapTypeStr = "2Fo-Fc";
    } else {
      mapTypeStr = "Fo-Fc";
    }    
  }
  
  if ((isModelString(phiStr.c_str()) || isModelString(fcStr.c_str())) && phiStr != fcStr) {
    if (phiLastSet) {
      fc_userset[mapNum-1] = false;
      if (isModelString(phiStr.c_str())) {
        fcStr = phiStr;
      } else {
        fcStr = findPreferredFcItem(foStr, phiStr, fcItems, _models);
      }
    } else {
      phi_userset[mapNum-1] = false;
      if (isModelString(fcStr.c_str())) {
        phiStr = fcStr;
      } else {
        phiStr = findPreferredPhiItem(foStr, fcStr, phiItems, _models);
      }
    }
  }

  // Set values if not set by user
  if (!mapType_userset[mapNum-1] && !mapTypeStr.empty()) {
    setStringSelection(mapTypeChoice,mapTypeStr);
  }
  if (!fo_userset[mapNum-1] && !foStr.empty()) {
    setStringSelection(foChoice,foStr);
  }
  if (!fc_userset[mapNum-1] && !fcStr.empty()) {
    setStringSelection(fcChoice,fcStr);
  }
  if (!phi_userset[mapNum-1] && !phiStr.empty()) {
    setStringSelection(phiChoice,phiStr);
  }

  int i = modelStringToNumber(phiStr.c_str());
  if (i >= 0) {
    std::string s=::format("%s = %s", phiStr.c_str(), _models[i].c_str());
    modelText->setText(s.c_str());
  } else {
    modelText->setText("");    
  }
  
  if (mapTypeStr == std::string("Direct FFT")) {
    // Must call ShowExtras() after setting label because label size may change
    // and need new layout
    if (mapNum == 1) {
      map1F1Label->setText(F_TEXT);
    } else {
      map2F1Label->setText(F_TEXT);
    }
    ShowExtras(mapNum, false);
  } else {
    // Must call ShowExtras() after setting label because label size may change
    // and need new layout
    if (mapNum == 1) {
      map1F1Label->setText(FO_TEXT);
    } else {
      map2F1Label->setText(FO_TEXT);
    }
    ShowExtras(mapNum, true);
  }
}


void PhaseFileLoadDialog::UpdateResolution()
{
  QComboBox* fcChoice=map1FcChoice;
  float resmin_value, resmax_value;

  int modelnum=modelStringToNumber(fcChoice->currentText().toStdString().c_str());
  if (modelnum == -1) {
    resmin_value=_init_resmin;
    resmax_value=_init_resmax;
  } else {
    
    CMapHeaderBase hdr;
    sscanf(_cells[modelnum].c_str(),"%f %f %f %f %f %f",
           &hdr.a, &hdr.b, &hdr.c,
           &hdr.alpha, &hdr.beta, &hdr.gamma);
    *_res_tmp_map->mapheader = hdr;
    _res_tmp_map->RecalcResolution();
    // Note map header uses the numerical min/max for resmin/resmax
    // rather than the crystallography resolution min/max 
    resmin_value=_res_tmp_map->mapheader->resmax;
    resmax_value=_res_tmp_map->mapheader->resmin;
  }

  //NOTE: I don't know how to detect that the user has typed text, so
  //      for now _userset will stay false
  if (!resmin_userset)
  {
    resmin->setValue(resmin_value);
  }

  if (!resmax_userset)
  {
    resmax->setValue(resmax_value);
  }
}



bool PhaseFileLoadDialog::SetFile(const std::string& fname,
                                  const std::vector<std::string>& models,
                                  const std::vector<std::string>& cells) {
  _filename = fname;
  _models = models;
  _cells = cells;
  mapType_userset[0] = false;
  mapType_userset[1] = false;
  fo_userset[0] = false;
  fo_userset[1] = false;
  fc_userset[0] = false;
  fc_userset[1] = false;
  phi_userset[0] = false;
  phi_userset[1] = false;
  fom_userset[0] = false;
  fom_userset[1] = false;
  resmin_userset = false;
  resmax_userset = false;

  _showFom = true;

  _res_tmp_map=new EMap;
  _res_tmp_map->UseColumnLabels("F","FC","FOM","PHI","","FreeR_flag"); // don't care what these are set to, the mere fact that they're set will cause mtz reader to work w/o going interactive
  if (cells.size()) {
    CMapHeaderBase hdr;
    sscanf(cells[0].c_str(),"%f %f %f %f %f %f",
           &hdr.a, &hdr.b, &hdr.c,
           &hdr.alpha, &hdr.beta, &hdr.gamma);
    *_res_tmp_map->mapheader = hdr;
  }
  _res_tmp_map->LoadMapPhaseFile(fname.c_str(), false); // don't require Fo!
  // Note map header uses the numerical min/max for resmin/resmax
  // rather than the crystallography resolution min/max 
  _init_resmin=_res_tmp_map->mapheader->resmax;
  _init_resmax=_res_tmp_map->mapheader->resmin;
  

  std::vector<unsigned int> without_fc_types;
  std::vector<unsigned int> with_fc_types;  

  if (!EMapBase::GetPossibleMapTypes(_filename,without_fc_types,with_fc_types))
    return false;
  
  // remove DirectFFT option for non-mtz files
  if (!EMapBase::IsCCP4MTZFile(fname.c_str())) {
    std::vector<unsigned int>::iterator i
        = std::find(without_fc_types.begin(), without_fc_types.end(), MIMapType::DirectFFT);
    if (i != without_fc_types.end()) {
      without_fc_types.erase(i);
    }
    i = std::find(with_fc_types.begin(), with_fc_types.end(), MIMapType::DirectFFT);
    if (i != with_fc_types.end()) {
      with_fc_types.erase(i);
    }
  }
  
  if (!_models.size() && !without_fc_types.size()) {
    MIMessageBox("No source of phases for file: " + fname +".  Load model first.");
    return false;
  }

  filename->setText(fname.c_str());
  gridChoice->setCurrentIndex(1);

  map1Type->clear();
  map2Type->clear();
  
  //note: if no model, use without_fc_types, not with_fc_types
  std::vector<unsigned int>* types=&with_fc_types;
  if (_models.size()==0)
    types=&without_fc_types;
  for (size_t i=0; i < types->size(); ++i) {
    map1Type->addItem(StringForMapType((*types)[i]));
    map2Type->addItem(StringForMapType((*types)[i]));
  }

  // now set column fields based on file type
  std::vector<std::string> file_fcs;
  std::vector<std::string> file_fos;
  std::vector<std::string> file_phis;
  std::vector<std::string> file_foms;
  if (EMapBase::IsCif(fname.c_str())) {
    //file possibly has fc: give user option of fc from file or model
    if (without_fc_types.size() == with_fc_types.size()) 
      file_fcs.push_back("Fc");
    _showFom = false;
  } else if (EMapBase::IsPhsFile(fname.c_str())) {
    //file has fc: give user option of fc from file or model
    file_fcs.push_back("Fc");
    _showFom = false;

  } else if (EMapBase::IsScaFile(fname.c_str()) || 
             EMapBase::IsFinFile(fname.c_str()) || 
             EMapBase::IsRefFile(fname.c_str())) {
    // no decisions to make except which model to use for fc
    _showFom = false;

  } else if (EMapBase::IsCCP4MTZFile(fname.c_str())) {
    // possibly lots of decisions to make, depending on map type?
    std::vector<std::string> labels;
    std::vector<char> field_types;
    if (!EMapBase::GetMTZColumnInfo(fname, labels, field_types)) {
      MIMessageBox("Can't get MTZ file column info for: " + fname);
      return false;
    }
    
    for (size_t i=0; i < labels.size(); ++i) {
      switch (field_types[i]) {
        case 'F':
        case 'G':
        case 'J':
        case 'K':
          file_fos.push_back(labels[i]);
          file_fcs.push_back(labels[i]);
          break;
        case 'W':
          file_foms.push_back(labels[i]);
          break;
        case 'P':
          file_phis.push_back(labels[i]);
          break;
        default:
          break;
      }
    }

  } else {
    MIMessageBox("Unrecognized phase file type: " + fname);
    return false;
  }

  if (file_fos.size() == 0) {
    file_fos.push_back("From file");
  }
  if (file_phis.size() == 0) {
    file_phis.push_back("From file");
  }
  std::vector<std::string> empty;
  SetChoiceOptions(map1FoChoice, file_fos, empty);
  SetChoiceOptions(map2FoChoice, file_fos, empty);
  SetChoiceOptions(map1FcChoice, file_fcs, _models);
  SetChoiceOptions(map2FcChoice, file_fcs, _models);
  SetChoiceOptions(map1PhiChoice, file_phis, _models);
  SetChoiceOptions(map2PhiChoice, file_phis, _models);
  SetChoiceOptions(map1FomChoice, file_foms, empty);
  SetChoiceOptions(map2FomChoice, file_foms, empty);

  

  updateDefaults(1);
  updateDefaults(2);
  UpdateResolution();

  //FIXME: don't enable Ok button this until dialog checks out (all columns set correctly)
  return true;
}


void PhaseFileLoadDialog::GetData(MIData& data) {
  data["grid"].radio = gridChoice->currentIndex();
  // Note map header uses the numerical min/max for resmin/resmax
  // rather than the crystallography resolution min/max 
  data["resmin"].f = resmax->value();
  data["resmax"].f = resmin->value();

  data["enabled1"].b = map1Check->isChecked();
  data["type1"].str = map1Type->currentText().toStdString();
  data["fo1"].str = map1FoChoice->currentText().toStdString();
  data["fc1"].str = modelStringToModelName(map1FcChoice->currentText().toStdString().c_str());
  data["phi1"].str = modelStringToModelName(map1PhiChoice->currentText().toStdString().c_str());
  data["fom1"].str = map1FomChoice->currentText().toStdString();

  data["enabled2"].b = map2Check->isChecked();
  data["type2"].str = map2Type->currentText().toStdString();
  data["fo2"].str = map2FoChoice->currentText().toStdString();
  data["fc2"].str = modelStringToModelName(map2FcChoice->currentText().toStdString().c_str());
  data["phi2"].str = modelStringToModelName(map2PhiChoice->currentText().toStdString().c_str());
  data["fom2"].str = map2FomChoice->currentText().toStdString();

#if 1
  Logger::debug("Dialog returning grid type %d, resmin %f resmax %f\n"
         "map1: enabled=%d type=%s fc=%s fo=%s fom=%s phi=%s\n"
         "map2: enabled=%d type=%s fc=%s fo=%s fom=%s phi=%s\n",
         data["grid"].radio, data["resmin"].f,data["resmax"].f,
         data["enabled1"].b,
         data["type1"].str.c_str(),
         data["fc1"].str.c_str(),
         data["fo1"].str.c_str(),
         data["fom1"].str.c_str(),
         data["phi1"].str.c_str(),
         data["enabled2"].b ,
         data["type2"].str.c_str(),
         data["fc2"].str.c_str(),
         data["fo2"].str.c_str(),
         data["fom2"].str.c_str(),
         data["phi2"].str.c_str());
#endif
}

std::string PhaseFileLoadDialog::modelStringToModelName(const char* str) {
  std::string result;
  int i = modelStringToNumber(str);
  if (i >= 0) {
    result = "Model:" + _models[i];
  } else {
    result = str;
  }
  return result;
}



void PhaseFileLoadDialog::on_map1Type_activated(int) { OnMapTypeChange(1); }
void PhaseFileLoadDialog::on_map2Type_activated(int) { OnMapTypeChange(2); }
void PhaseFileLoadDialog::on_map1FoChoice_activated(int) { OnMapFoChange(1); }
void PhaseFileLoadDialog::on_map2FoChoice_activated(int) { OnMapFoChange(2); }

void PhaseFileLoadDialog::on_map1FcChoice_activated(int) { OnMapFcChange(1); }
void PhaseFileLoadDialog::on_map2FcChoice_activated(int) { OnMapFcChange(2); }
void PhaseFileLoadDialog::on_map1PhiChoice_activated(int) { OnMapPhiChange(1); }
void PhaseFileLoadDialog::on_map2PhiChoice_activated(int) { OnMapPhiChange(2); }
void PhaseFileLoadDialog::on_map1FomChoice_activated(int) { OnMapFomChange(1); }
void PhaseFileLoadDialog::on_map2FomChoice_activated(int) { OnMapFomChange(2); }
