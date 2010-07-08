#ifndef PHASE_FILE_LOAD_DIALOG_H
#define PHASE_FILE_LOAD_DIALOG_H

#include <string>
#include <vector>

#include "ui/uilib.h"

#include "ui_PhaseFileLoadDialog.h"

class PhaseFileLoadDialog : public QDialog, public Ui::PhaseFileLoadDialog
{
    Q_OBJECT

    std::string _filename;
    std::vector<std::string> _models;
    std::vector<std::string> _cells;
    bool mapType_userset[2];
    bool fo_userset[2];
    bool fc_userset[2];
    bool phi_userset[2];
    bool fom_userset[2];
    bool resmin_userset;
    bool resmax_userset;
    bool _showFom;

    EMap *_res_tmp_map;
    float _init_resmin, _init_resmax;

    void UpdateResolution();

    void OnMapTypeChange(unsigned int i);
    void OnMapFoChange(unsigned int i);
    void OnMapFcChange(unsigned int i);
    void OnMapPhiChange(unsigned int i);
    void OnMapFomChange(unsigned int i);

    void updateCheckBoxes(unsigned int num);
    void updateDefaults(unsigned int mapNum);
    void ShowExtras(unsigned int mapNum, bool show = true);
    std::string modelStringToModelName(const char *str);

private slots:
    void on_map1Type_activated(int);
    void on_map2Type_activated(int);
    void on_map1FoChoice_activated(int);
    void on_map2FoChoice_activated(int);

    void on_map1FcChoice_activated(int);
    void on_map2FcChoice_activated(int);
    void on_map1PhiChoice_activated(int);
    void on_map2PhiChoice_activated(int);
    void on_map1FomChoice_activated(int);
    void on_map2FomChoice_activated(int);

public:

    struct Data
    {
        int grid;

        // Note map header uses the numerical min/max for resmin/resmax
        // rather than the crystallography resolution min/max
        float resmin;
        float resmax;

        bool enabled1;
        std::string type1;
        std::string fo1;
        std::string fc1;
        std::string phi1;
        std::string fom1;

        bool enabled2;
        std::string type2;
        std::string fo2;
        std::string fc2;
        std::string phi2;
        std::string fom2;
    };

    PhaseFileLoadDialog(QWidget *parent);
    ~PhaseFileLoadDialog();

    bool SetFile(const std::string &fname,
                 const std::vector<std::string> &models,
                 const std::vector<std::string> &cells);
    void GetData(Data &data);
};


#endif // ifndef PHASE_FILE_LOAD_DIALOG_H
