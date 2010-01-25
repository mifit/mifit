#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <map/maplib.h>
#include <QMessageBox>
#include "GenericDataDialog.h"
#include "core/corelib.h"
#include "EMap.h"
#include "macafxwin.h"
#include "Application.h"

#include "MapSettings.h"
#include "ui/ContourOptions.h"
#include "ui/SelectCrystal.h"
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>

struct mmtz_column_;

#ifdef _WIN32
#define _MVS
#ifndef i386
#define i386
#endif
#include <umtz/mmtzlib.h>
#undef _MVS
#else
#include <umtz/mmtzlib.h>
#endif






void EMap::Save(CArchive &ar, int i)
{
    char tmpbuf[1000];
    if (HasPhases())
    {
        sprintf(tmpbuf, "MapColumns %s%s %s%s %s%s %s%s %s%s %s%s\n",
                (_fostr.size() ? "FO=" : ""), (_fostr.size() ? _fostr.c_str() : ""),
                (_fcstr.size() ? "FC=" : ""), (_fcstr.size() ? _fcstr.c_str() : ""),
                (_fomstr.size() ? "FOM=" : ""), (_fomstr.size() ? _fomstr.c_str() : ""),
                (_phistr.size() ? "PHI=" : ""), (_phistr.size() ? _phistr.c_str() : ""),
                (_sigfstr.size() ? "SIGF=" : ""), (_sigfstr.size() ? _sigfstr.c_str() : ""),
                (_freeRstr.size() ? "FREE=" : ""), (_freeRstr.size() ? _freeRstr.c_str() : ""));
        ar.Write(tmpbuf, strlen(tmpbuf));
        sprintf(tmpbuf, "LoadMapPhase %d %s\n", i, pathName.c_str());
        ar.Write(tmpbuf, strlen(tmpbuf));

        sprintf(tmpbuf, "unitcell %f %f %f %f %f %f\n", mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta, mapheader->gamma);
        ar.Write(tmpbuf, strlen(tmpbuf));
        sprintf(tmpbuf, "spacegroupno %d\n", mapheader->spgpno);
        ar.Write(tmpbuf, strlen(tmpbuf));
        /*    sprintf(tmpbuf,"crystal %s\n", mapheader->crystal_name.c_str());
            ar.Write(tmpbuf, strlen(tmpbuf));
            sprintf(tmpbuf,"ctitle %s\n", mapheader->title.c_str());
            ar.Write(tmpbuf, strlen(tmpbuf));
         */
        sprintf(tmpbuf, "coefficients %s\n", StringForMapType(mapheader->maptype));
        ar.Write(tmpbuf, strlen(tmpbuf));
        sprintf(tmpbuf, "fftnx %d\nfftny %d\nfftnz %d\nresmin %0.2f\nresmax %0.2f\n", mapheader->nx,
                mapheader->ny, mapheader->nz, mapheader->resmin, mapheader->resmax);
        ar.Write(tmpbuf, strlen(tmpbuf));
        strcpy(tmpbuf, "fftapply\n");
        ar.Write(tmpbuf, strlen(tmpbuf));
    }
    else
    {
        sprintf(tmpbuf, "LoadMap %d %s\n", i, pathName.c_str());
        ar.Write(tmpbuf, strlen(tmpbuf));
    }
    sprintf(tmpbuf, "maptocont %d\n", i);
    ar.Write(tmpbuf, strlen(tmpbuf));
    sprintf(tmpbuf, "maplinewidth %f\n", settings->maplinewidth);
    ar.Write(tmpbuf, strlen(tmpbuf));
    sprintf(tmpbuf, "contourlevels %d\n",
            (settings->MapLevelOn[0] ? 1 : 0)
            +2*(settings->MapLevelOn[1] ? 1 : 0)
            +4*(settings->MapLevelOn[2] ? 1 : 0)
            +8*(settings->MapLevelOn[3] ? 1 : 0)
            +16*(settings->MapLevelOn[4] ? 1 : 0)   );
    ar.Write(tmpbuf, strlen(tmpbuf));
    sprintf(tmpbuf, "contourleveldefault %f %f %f %f %f\n",
            settings->MapLevel[0], settings->MapLevel[1], settings->MapLevel[2], settings->MapLevel[3], settings->MapLevel[4]);
    ar.Write(tmpbuf, strlen(tmpbuf));
    sprintf(tmpbuf, "contourradius %f\n", settings->Radius);
    ar.Write(tmpbuf, strlen(tmpbuf));
    for (int j = 0; j < 5; j++)
    {
        sprintf(tmpbuf, "color %d\n", settings->MapLevelColor[j]);
        ar.Write(tmpbuf, strlen(tmpbuf));
        sprintf(tmpbuf, "contourcolor %d\n", j+1);
        ar.Write(tmpbuf, strlen(tmpbuf));
    }
    sprintf(tmpbuf, "contourmap %d\n\n", i);
    ar.Write(tmpbuf, strlen(tmpbuf));
}

bool EMap::ContourLevels()
{
    MIData data;
    MapSettingsToData(*settings, data, mapmin, mapmax);

    QDialog dlg(0);

    dlg.setWindowTitle(MapID().c_str());
    dlg.setModal(true);
    dlg.setSizeGripEnabled(true);

    ContourOptions *co = new ContourOptions(&dlg, false);
    QDialogButtonBox *bb = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, &dlg);
    dlg.connect(bb, SIGNAL(accepted()), &dlg, SLOT(accept()));
    dlg.connect(bb, SIGNAL(rejected()), &dlg, SLOT(reject()));

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(co);
    mainLayout->addWidget(bb);
    dlg.setLayout(mainLayout);

    co->InitializeFromData(data);
    if (dlg.exec() != QDialog::Accepted)
    {
        return false;
    }
    co->GetData(data);

    float mapmin, mapmax;
    DataToMapSettings(data, *settings, mapmin, mapmax);
    SetMapLinewidth(settings->maplinewidth);
    //settings->save(); // No longer relevant.  only the preferences are saved
    mapContourLevelsChanged(this);
    return true;
}

long EMap::LoadMap(const char *pathname, int type)
{
    if (HasPhases())
    {
        if (QMessageBox::question(0, "Overwrite?",
                                  "This map already has a reflection list\nDo you want to overwrite it?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
        {
            return 0;
        }
    }
    return EMapBase::LoadMap(pathname, type);
}

long EMap::LoadCIFMap(const char *pathname, int datablock)
{
    if (refls.size() != 0)
    {
        if (QMessageBox::question(0, "Overwrite?",
                                  "This map already has a reflection list\nDo you want to overwrite it?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
        {
            return 0;
        }
    }
    return EMapBase::LoadCIFMap(pathname, datablock);
}

bool EMap::PromptForCrystal()
{
    MIData data;
    data["info"].str = CMapHeaderBase().Label();
    SelectCrystal::doSelectCrystal(data);
    CMapHeaderBase mh(data["info"].str);
    mapheader->updateSymmetryAndCell(mh);
    RecalcResolution();
    return true;
}

bool EMap::PromptForColumnLabels(unsigned int, mmtz_column_*,
                                 int&, int&, int&, int&, int&, int&)
{
    // no longer used, as column labels are now set in advance.
    return false;
}

EMap::EMap()
    : EMapBase()
{
    SetVersion(MIFit_version);

    // replace MapSettingsBase with MapSettings (which has GUI load & save support)
    // this will be deleted by EMapBase dtor
    delete settings;
    settings = new MapSettings();
    settings->load();

    // replace CMapHeaderBase with CMapHeader (which has GUI load & save support)
    // this will be deleted by EMapBase dtor
    delete mapheader;
    mapheader = new CMapHeader();
}

bool EMap::FFTMap(int mt, int gl, float rMin, float rMax)
{

    if (!HasPhic())
    {
        Logger::message("Before you can display the map, you need to use the\n"
                        "Calculate Structure Factors... command to\n"
                        "calculate phases from a model.");
        return false;
    }

    GenericDataDialog dlg;
    dlg.setWindowTitle("Fast Fourier Transform Map");

    QStringList mapTypeOptions;
    int mapTypeIndex = 0;
    QList<unsigned int> mapTypes;

    int currentMapType = (mt != -1) ? mt : mapheader->maptype;
    if (currentMapType == (int)MIMapType::DirectFFT)
    {
        // If DirectFFT, it can only be DirectFFT
        mapTypeOptions += StringForMapType(MIMapType::DirectFFT);
        mapTypes.push_back(MIMapType::DirectFFT);
    }
    else
    {
        // Limit map types to types which can be calculated given the fields loaded
        std::vector<unsigned int> types;
        std::vector<unsigned int> types_with_fc;
        GetMapTypes(types, types_with_fc);

        for (size_t i = 0; i < types.size(); ++i)
        {
            if (types[i] != MIMapType::DirectFFT)
            {
                // Exclude DirectFFT from list
                mapTypes += types[i];
                mapTypeOptions += StringForMapType(types[i]);
                if ((int)types[i] == currentMapType)
                {
                    mapTypeIndex = mapTypes.size()-1;
                }
            }
        }
    }

    QStringList gridTypes;
    gridTypes << "Coarse grid" << "Medium grid" << "Fine grid";
    int gridIndex = (gl == -1 || gl == -2) ? mapGrid : gl;

    dlg.addDoubleField("Min resolution:", rMin == -1.0f ? mapheader->resmin : rMin);
    dlg.addDoubleField("Max resolution:", rMax == -1.0f ? mapheader->resmax : rMax);
    dlg.addComboField("Map type:", mapTypeOptions, mapTypeIndex);
    dlg.addComboField("Grid:", gridTypes, gridIndex);


    if (mt!=-1 || (gl!=-1 && gl!=-2) || Application::instance()->GetSilentMode() || dlg.exec() == QDialog::Accepted)
    {
        bool result = EMapBase::FFTMap(mapTypes[dlg.value(2).toInt()],
                                       (gl==-2 ? -1 : dlg.value(3).toInt()),
                                       dlg.value(0).toDouble(),
                                       dlg.value(1).toDouble());
        mapFftRecalculated(this);
        return result;
    }

    return false;
}

void EMap::Export()
{
    if (HasPhases())
    {
        QString selectedFilter;
        QString path = QFileDialog::getSaveFileName(0, "Choose a name for the phase file", "",
                                                    "CCP4 phase file (*.mtz);;"
                                                    "mmCIF phase file (*.cif);;"
                                                    "Xtalview phase file (*.phs)", &selectedFilter);
        if (path.isEmpty())
            return;

        int format = CCP4_phase;
        if (selectedFilter == "CCP4 phase file (*.mtz)")
            format = CCP4_phase;
        else if (selectedFilter == "mmCIF phase file (*.cif)")
            format = mmCIF_phase;
        else if (selectedFilter == "Xtalview phase file (*.phs)")
            format = XtalView_phase;
        SavePhases(path.toAscii().constData(), format);
    }
}

int EMap::CenterVisibleEdges(float &x, float &y, float &z, ViewPoint *vp)
{
    // look thru the edges, find the on-scrren ones abd return their center
    x = y = z = 0;
    size_t n = 0;
    for (size_t i = 0; i < edges.size(); i++)
    {
        if (IsOnScreen(edges[i].p1, vp))
        {
            x += edges[i].p1.x;
            y += edges[i].p1.y;
            z += edges[i].p1.z;
            n++;
        }
        if (IsOnScreen(edges[i].p2, vp))
        {
            x += edges[i].p2.x;
            y += edges[i].p2.y;
            z += edges[i].p2.z;
            n++;
        }
    }
    x = x/(float)n;
    y = y/(float)n;
    z = z/(float)n;
    return n;
}

