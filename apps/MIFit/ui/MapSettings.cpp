#include "core/corelib.h"
#include "MapSettings.h"

#include "core/MIConfig.h"
#include "core/Colors.h"
#include <util/utillib.h>


static bool GetHardcodedDefaults(MapSettingsBase &settings, int styleno)
{
  //these are styleno-independent
  settings.m_radiusmax=50;
  settings.Radius=8.0;

  settings.ContourMethod=MAP_CUBE;
  settings.BlobRadius=1.7f;
  settings.maplinewidth=1.0;

  //these are styleno-dependent
  switch (styleno) {
    case 0:  // Blue Map 1,2,3,4,5 sigma
      settings.MapLevelColor[0]= Colors::MAP1;
      settings.MapLevelColor[1]= Colors::MAP2;
      settings.MapLevelColor[2]= Colors::MAP3;
      settings.MapLevelColor[3]= Colors::MAP4;
      settings.MapLevelColor[4]= Colors::MAP5;
      settings.MapLevelOn[0]= true;
      settings.MapLevelOn[1]= true;
      settings.MapLevelOn[2]= true;
      settings.MapLevelOn[3]= true;
      settings.MapLevelOn[4]= true;
      settings.MapLevel[0]= 50.0f;
      settings.MapLevel[1]= 100.0f;
      settings.MapLevel[2]= 150.0f;
      settings.MapLevel[3]= 200.0f;
      settings.MapLevel[4]= 250.0f;
      break;
    case 1:  // Single Level Blue 1 sigma
      settings.MapLevelColor[0]= Colors::MAP1;
      settings.MapLevelColor[1]= Colors::MAP2;
      settings.MapLevelColor[2]= Colors::MAP3;
      settings.MapLevelColor[3]= Colors::MAP4;
      settings.MapLevelColor[4]= Colors::MAP5;
      settings.MapLevelOn[0]= true;
      settings.MapLevelOn[1]= false;
      settings.MapLevelOn[2]= false;
      settings.MapLevelOn[3]= false;
      settings.MapLevelOn[4]= false;
      settings.MapLevel[0]= 50.0f;
      settings.MapLevel[1]= 100.0f;
      settings.MapLevel[2]= 150.0f;
      settings.MapLevel[3]= 200.0f;
      settings.MapLevel[4]= 250.0f;
      break;
    case 2: // Green Map 1,2,3,4,5 Level
      settings.MapLevelColor[0]= Colors::MAP6;
      settings.MapLevelColor[1]= Colors::MAP7;
      settings.MapLevelColor[2]= Colors::MAP8;
      settings.MapLevelColor[3]= Colors::MAP9;
      settings.MapLevelColor[4]= Colors::MAP10;
      settings.MapLevelOn[0]= true;
      settings.MapLevelOn[1]= true;
      settings.MapLevelOn[2]= true;
      settings.MapLevelOn[3]= true;
      settings.MapLevelOn[4]= true;
      settings.MapLevel[0]= 50.0f;
      settings.MapLevel[1]= 100.0f;
      settings.MapLevel[2]= 150.0f;
      settings.MapLevel[3]= 200.0f;
      settings.MapLevel[4]= 250.0f;
      break;
    case 3: // Single Level Green 1 sigma
      settings.MapLevelColor[0]= Colors::MAP6;
      settings.MapLevelColor[1]= Colors::MAP7;
      settings.MapLevelColor[2]= Colors::MAP8;
      settings.MapLevelColor[3]= Colors::MAP9;
      settings.MapLevelColor[4]= Colors::MAP10;
      settings.MapLevelOn[0]= true;
      settings.MapLevelOn[1]= false;
      settings.MapLevelOn[2]= false;
      settings.MapLevelOn[3]= false;
      settings.MapLevelOn[4]= false;
      settings.MapLevel[0]= 50.0f;
      settings.MapLevel[1]= 100.0f;
      settings.MapLevel[2]= 150.0f;
      settings.MapLevel[3]= 200.0f;
      settings.MapLevel[4]= 250.0f;
      break;
    case 4: // Difference Map -4,-3,3,4,5 sigma
      settings.MapLevelColor[0]= Colors::MAP4;
      settings.MapLevelColor[1]= Colors::MAP5;
      settings.MapLevelColor[2]= Colors::MAP1;
      settings.MapLevelColor[3]= Colors::MAP2;
      settings.MapLevelColor[4]= Colors::MAP3;
      settings.MapLevelOn[0]= true;
      settings.MapLevelOn[1]= true;
      settings.MapLevelOn[2]= true;
      settings.MapLevelOn[3]= true;
      settings.MapLevelOn[4]= true;
      settings.MapLevel[0]= -200.0f;
      settings.MapLevel[1]= -150.0f;
      settings.MapLevel[2]= 150.0f;
      settings.MapLevel[3]= 200.0f;
      settings.MapLevel[4]= 250.0f;
      break;
    case 5: // Difference Map -3 red,+3 blue sigma
      settings.MapLevelColor[0]= Colors::MAP4;
      settings.MapLevelColor[1]= Colors::MAP5;
      settings.MapLevelColor[2]= Colors::MAP1;
      settings.MapLevelColor[3]= Colors::MAP2;
      settings.MapLevelColor[4]= Colors::MAP3;
      settings.MapLevelOn[0]= false;
      settings.MapLevelOn[1]= true;
      settings.MapLevelOn[2]= true;
      settings.MapLevelOn[3]= false;
      settings.MapLevelOn[4]= false;
      settings.MapLevel[0]= -200.0f;
      settings.MapLevel[1]= -150.0f;
      settings.MapLevel[2]= 150.0f;
      settings.MapLevel[3]= 200.0f;
      settings.MapLevel[4]= 250.0f;
      break;
    default:
      return false;
      break;
  }
  return true;
}



void MapSettings::saveStyle(int i) {
  std::string key;
  std::string stylestr("");
  if (i!=-1)
    stylestr=format("%d",i);
  const char *sty=stylestr.c_str();

  key = format("MapSettings%s/m_radiusmax", sty);
  MIConfig::Instance()->Write(key, (long)m_radiusmax);
  key=format("MapSettings%s/Radius", sty);
  MIConfig::Instance()->Write(key, Radius);
  key=format("MapSettings%s/ContourMethod", sty);
  MIConfig::Instance()->Write(key, (long)ContourMethod);
  key=format("MapSettings%s/BlobRadius", sty);
  MIConfig::Instance()->Write(key, BlobRadius);
  key = format("MapSettings%s/maplinewidth", sty);
  MIConfig::Instance()->Write(key, maplinewidth);

  for (int i = 0; i < 5; ++i) {
    key = format("MapSettings%s/MapLevelColor%d", sty, i);
    MIConfig::Instance()->Write(key, (long)MapLevelColor[i]);
  }

  for (int i = 0; i < 5; ++i) {
    key = format("MapSettings%s/MapLevelOn%d", sty, i);
    MIConfig::Instance()->Write(key, MapLevelOn[i]);
  }

  for (int i = 0; i < 5; ++i) {
    key = format("MapSettings%s/MapLevel%d", sty, i);
    MIConfig::Instance()->Write(key, (double)MapLevel[i]);
  }
}

void MapSettings::loadStyle(int i) {

  // first get the hard-coded defaults
  GetHardcodedDefaults(*this,i);

  // now try in the user prefs file, using the hardcoded values as defaults
  // in case they're not in the settings file
  std::string key;
  std::string stylestr("");
  if (i!=-1)
    stylestr=format("%d",i);
  const char *sty=stylestr.c_str();

  long longValue;
  double d;

  key = format("MapSettings%s/m_radiusmax", sty);
  MIConfig::Instance()->Read(key, &longValue, (long)m_radiusmax);
  m_radiusmax = longValue;
  key=format("MapSettings%s/Radius",sty);
  MIConfig::Instance()->Read(key, &d, Radius);
  Radius = (float) d;
  key=format("MapSettings%s/ContourMethod",sty);
  MIConfig::Instance()->Read(key, &longValue, (long)ContourMethod);
  ContourMethod = longValue;
  key=format("MapSettings%s/BlobRadius",sty);
  MIConfig::Instance()->Read(key, &d, BlobRadius);
  BlobRadius = (float) d;
  key = format("MapSettings%s/maplinewidth", sty);
  MIConfig::Instance()->Read(key, &d, maplinewidth);
  maplinewidth = d;

  for (int i = 0; i < 5; ++i) {
    key = format("MapSettings%s/MapLevelColor%d", sty, i);
    MIConfig::Instance()->Read(key, &longValue, (long)MapLevelColor[i]);
    MapLevelColor[i] = longValue;
  }
  for (int i = 0; i < 5; ++i) {
    key = format("MapSettings%s/MapLevelOn%d", sty, i);
    MIConfig::Instance()->Read(key, &MapLevelOn[i], MapLevelOn[i]);
  }
  for (int i = 0; i < 5; ++i) {
    key = format("MapSettings%s/MapLevel%d", sty, i);
    MIConfig::Instance()->Read(key, &d, MapLevel[i]);
    MapLevel[i] = d;
  }
}

void MapSettings::save() {
  saveStyle(-1);
}

void MapSettings::load() {
  loadStyle(-1);
}

// copy every value in settings to MIData structure, *and* add mapmin, mapmax to data structure
void MapSettingsToData(const MapSettingsBase &settings, MIData &data, float mapmin, float mapmax)
{
  data["mapmin"].f=mapmin;
  data["mapmax"].f=mapmax;

  data["radiusMax"].i=settings.m_radiusmax;
  data["radius"].f=settings.Radius;
  data["contourMethod"].radio=settings.ContourMethod;
  data["contourMethod"].radio_count=3; // box, sphere, blob
  data["blobRadius"].f=settings.BlobRadius;
  data["lineWidth"].f=settings.maplinewidth;

  data["MapLevelColor0"].i=settings.MapLevelColor[0];
  data["MapLevelColor1"].i=settings.MapLevelColor[1];
  data["MapLevelColor2"].i=settings.MapLevelColor[2];
  data["MapLevelColor3"].i=settings.MapLevelColor[3];
  data["MapLevelColor4"].i=settings.MapLevelColor[4];

  data["MapLevelOn0"].b=settings.MapLevelOn[0];
  data["MapLevelOn1"].b=settings.MapLevelOn[1];
  data["MapLevelOn2"].b=settings.MapLevelOn[2];
  data["MapLevelOn3"].b=settings.MapLevelOn[3];
  data["MapLevelOn4"].b=settings.MapLevelOn[4];

  data["MapLevel0"].f=settings.MapLevel[0];
  data["MapLevel1"].f=settings.MapLevel[1];
  data["MapLevel2"].f=settings.MapLevel[2];
  data["MapLevel3"].f=settings.MapLevel[3];
  data["MapLevel4"].f=settings.MapLevel[4];
}

// copy every value in MIData structure to settings  *and* set mapmin, mapmax 
void DataToMapSettings(/* const */ MIData &data, MapSettingsBase &settings, float &mapmin, float &mapmax)
{
  mapmin=data["mapmin"].f;
  mapmax=data["mapmax"].f;

  settings.m_radiusmax = data["radiusMax"].i;
  settings.Radius = data["radius"].f;
  settings.ContourMethod = data["contourMethod"].radio;
  settings.BlobRadius = data["blobRadius"].f;
  settings.maplinewidth=data["lineWidth"].f;

  settings.MapLevelColor[0]=data["MapLevelColor0"].i;
  settings.MapLevelColor[1]=data["MapLevelColor1"].i;
  settings.MapLevelColor[2]=data["MapLevelColor2"].i;
  settings.MapLevelColor[3]=data["MapLevelColor3"].i;
  settings.MapLevelColor[4]=data["MapLevelColor4"].i;

  settings.MapLevelOn[0] = data["MapLevelOn0"].b;
  settings.MapLevelOn[1] = data["MapLevelOn1"].b;
  settings.MapLevelOn[2] = data["MapLevelOn2"].b;
  settings.MapLevelOn[3] = data["MapLevelOn3"].b;
  settings.MapLevelOn[4] = data["MapLevelOn4"].b;
  
  settings.MapLevel[0] = data["MapLevel0"].f;
  settings.MapLevel[1] = data["MapLevel1"].f;
  settings.MapLevel[2] = data["MapLevel2"].f;
  settings.MapLevel[3] = data["MapLevel3"].f;
  settings.MapLevel[4] = data["MapLevel4"].f;
}

bool GetMapSettingsForStyle(int styleno, MapSettingsBase &settings, bool hardcodedOnly)
{
  if (hardcodedOnly) {
    return GetHardcodedDefaults(settings,styleno);
  }
  
  MapSettings *p=dynamic_cast<MapSettings*>(&settings);
  if (!p)
    return false;
  p->loadStyle(styleno);
  return true;
}
