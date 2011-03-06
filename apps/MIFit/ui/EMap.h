#ifndef MI_EMap_H
#define MI_EMap_H

#include <QObject>
#include "core/corelib.h"
#include "CMapHeader.h"

class CArchive;
class XMLArchive;
class ViewPoint;

class EMap : public EMapBase
{
    Q_OBJECT
public:
    EMap();

    //@{
    // save the map to an archive.
    //@}
    void Save(CArchive &ar, int i);

    //@{
    // Get the map contouring levels with a dialog box.
    // Another function that would need to be changed to make this class non-interactive.
    //@}
    bool ContourLevels();

    // overloaded from base class to allow prompting user
    virtual long LoadMap(const char *pathname, int type = XtalView_phase);
    virtual long LoadCIFMap(const char *pathname, int datablock);

    bool PromptForCrystal();

    //@{
    // FFT the map.
    // Gets information from the user with a dialog box then calls FFTCalc.
    //@}
    virtual bool FFTMap(int maptype = -1, int gridlevel = -1,
                        float resMin = -1.0f, float resMax = -1.0f);

    void Export();

signals:
    void mapFftRecalculated(EMap*);

};


#endif // ifndef MI_EMap_H
