#ifndef MI_DIALOG_H
#define MI_DIALOG_H

#include "core/corelib.h"
#include <cfloat>
#include <climits>
#include <vector>
#include <string>

class QWidget;
class QDialog;

class MIDialog
{
public:
    MIDialog(QWidget *parent, const std::string &name);
    virtual ~MIDialog();

    bool GetResults(MIData &data);

protected:
    QWidget *_qparent;
    std::string _name;

private:
    virtual bool PromptForResults(MIData &data) = 0;
};


int MIGetSingleChoiceIndex(const std::string &message,
                           const std::string &caption,
                           const std::vector<std::string> &choices,
                           QWidget *parent = 0);


class MIColorPromptDialog : public MIDialog
{
public:
    MIColorPromptDialog(QWidget *parent, const std::string &name);
private:
    bool PromptForResults(MIData &data);
    //variables red.u green.u blue.u
};

class MIColorPaletteDialog : public MIDialog
{
public:
    MIColorPaletteDialog(QWidget *parent, const std::string &name);
private:
    bool PromptForResults(MIData &data);
    //variables: color.radio
};

//NOTE: wraps MIColorPaletteDialog
int MIColorChooser(int start_color, const std::string name = "Choose color");


const unsigned int MI_OPEN_MODE = 0;
const unsigned int MI_SAVE_MODE = 2;

class MIFileDialog : public MIDialog
{
public:
    MIFileDialog(QWidget *parent, const std::string &message,
                 const std::string &deftDir = "",
                 const std::string &deftFile = "",
                 const std::string &filter = "",
                 unsigned int mode = MI_OPEN_MODE);

private:
    bool PromptForResults(MIData &data);
    std::string _deftDir, _deftFile, _filter;
    unsigned int _mode;
    //variables: filterIndex.i path.str pathList.strlist
};

std::string MIFileSelector(const std::string &message, const std::string &default_path = "",
                           const std::string &default_filename = "", const std::string &default_extension = "",
                           const std::string &wildcard = "*.*", unsigned int mode = MI_OPEN_MODE, QWidget *parent = NULL);


#endif // ifndef MI_DIALOG_H

