#include <nongui/nonguilib.h>
#include "LSQFitDialog.h"
#include <QFileDialog>
#include "core/corelib.h"
#include <math/mathlib.h>
#include "ui/uilib.h"
#include <chemlib/RESIDUE_.h>

// for getting current displaylist
#include "ui/MIMainWindow.h"
#include "ui/MIGLWidget.h"

using namespace chemlib;

LSQFitDialog::LSQFitDialog(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    targetres = NULL;
    sourceres = NULL;
    m_target = NULL;
    m_source = NULL;
    setSuccess(false);
}

void LSQFitDialog::updateButtons()
{
    applyButton->setEnabled(success);
    exportButton->setEnabled(success);
}

void LSQFitDialog::setSuccess(bool value)
{
    success = value;
    updateButtons();
    updateMatrix();
}

void LSQFitDialog::updateMatrix()
{
    if (success)
    {
        std::string mess = ::format("|%5.3f %5.3f %5.3f|   |%6.3f|\n"
                                    "|%5.3f %5.3f %5.3f| + |%6.3f|\n"
                                    "|%5.3f %5.3f %5.3f|   |%6.3f|\n",
                                    r[0][0], r[0][1],  r[0][2], v[0],
                                    r[1][0], r[1][1],  r[1][2], v[1],
                                    r[2][0], r[2][1],  r[2][2], v[2]);
        matrixText->setText(mess.c_str());
        mess = ::format("%0.2f", rms);
        rmsText->setText(mess.c_str());
    }
    else
    {
        matrixText->setText("");
        rmsText->setText("");
    }
}


void LSQFitDialog::on_sourceListWidget_currentTextChanged(const QString &str)
{
    if (m_source == NULL || str.size()==0)
    {
        sourceres = NULL;
        sourceTextCtrl->setText("");
        return;
    }

    std::string resstr = str.toStdString();
    RESIDUE *res = m_source->getResidues();
    while (res != NULL)
    {
        if (resid(res) == resstr)
        {
            sourceres = res;
        }
        res = res->next();
    }
    sourceTextCtrl->setText(resstr.c_str());
}

void LSQFitDialog::on_targetListWidget_currentTextChanged(const QString &str)
{
    if (m_target == NULL || str.size()==0)
    {
        targetres = NULL;
        targetTextCtrl->setText("");
        return;
    }

    std::string resstr = str.toStdString();
    RESIDUE *res = m_target->getResidues();
    while (res != NULL)
    {
        if (resid(res) == resstr)
        {
            targetres = res;
        }
        res = res->next();
    }
    targetTextCtrl->setText(resstr.c_str());
}



void LSQFitDialog::on_targetComboBox_currentIndexChanged(const QString &text)
{
    Molecule *mol = findMolecule(text.toStdString());
    if (mol)
    {
        m_target = mol;
        targetres = m_target->getResidues();
        ListTarget();
    }
    Matches.clear();
    matchListBox->clear();
    setSuccess(false);
}

void LSQFitDialog::on_sourceComboBox_currentIndexChanged(const QString &text)
{
    Molecule *mol = findMolecule(text.toStdString());
    if (mol)
    {
        m_source = mol;
        sourceres = m_source->getResidues();
        ListSource();
    }
    Matches.clear();
    matchListBox->clear();
    setSuccess(false);
}

static void copya(double a[3], MIAtom *atom)
{
    a[0] = atom->x();
    a[1] = atom->y();
    a[2] = atom->z();
}

void LSQFitDialog::on_calcButton_clicked()
{
    int nmatch = Matches.size();
    double (*a)[3], (*b)[3];
    double *w;
    MATCH *m;
    MIAtom *atom;
    RESIDUE *tres, *sres, **sourcefrom;
    MIAtom **sourceatoms;
    MIAtom **targetatoms;
    int n = 0, i;
    for (i = 0; i < nmatch; i++)
    {
        n += Matches[i].length;
    }
    if (n < 3)
    {
        Logger::message("Must have at least 3 matching pairs of atoms to calculate least-quares fit");
        return;
    }
    a = new double[n][3];
    b = new double[n][3];
    sourcefrom = new (RESIDUE(*[n]));
    targetatoms = new (MIAtom(*[n]));
    sourceatoms = new (MIAtom(*[n]));
    w = new double[n];
    if (!w)
    {
        Logger::message("Can't allocate memory - try fewer matches");
        if (a)
        {
            delete[] a;
        }
        if (b)
        {
            delete[] b;
        }
        if (sourcefrom)
        {
            delete[] sourcefrom;
        }
        if (targetatoms)
        {
            delete[] targetatoms;
        }
        if (sourceatoms)
        {
            delete[] sourceatoms;
        }
        return;
    }
    n = 0;
    for (i = 0; i < nmatch; i++)
    {
        m = &(Matches[i]);
        tres = m->target;
        sres = m->source;
        for (int j = 0; j < m->length; j++)
        {
            if ((tres != NULL) && (sres != NULL))
            {
                if ((atom = atom_from_name(m->atomtype.c_str(), *tres)) != NULL)
                {
                    copya(b[n], atom);
                    targetatoms[n] = atom;
                    if ((atom = atom_from_name(m->atomtype.c_str(), *sres)) != NULL)
                    {
                        copya(a[n], atom);
                        sourcefrom[n] = sres;
                        sourceatoms[n] = atom;
                        w[n] = 1.0;
                        n++;
                        tres = tres->next();
                        sres = sres->next();
                    }
                }
            }
        }
    }
    if (n < 3)
    {
        Logger::message("Must have at least 3 matching pairs of atoms to calculate least-quares fit");
        if (a)
        {
            delete[] a;
        }
        if (b)
        {
            delete[] b;
        }
        if (w)
        {
            delete[] w;
        }
        if (sourcefrom)
        {
            delete[] sourcefrom;
        }
        if (targetatoms)
        {
            delete[] targetatoms;
        }
        if (sourceatoms)
        {
            delete[] sourceatoms;
        }
        return;
    }
    rotlsqfit(a, b, w, n, r, v);

    double x, y, z, tx, ty, tz, d;

    rms = 0.0;
    for (i = 0; i < n; i++)
    {
        tx = sourceatoms[i]->x();
        ty = sourceatoms[i]->y();
        tz = sourceatoms[i]->z();
        x = r[0][0]*tx + r[0][1]*ty + r[0][2]*tz + v[0];
        y = r[1][0]*tx + r[1][1]*ty + r[1][2]*tz + v[1];
        z = r[2][0]*tx + r[2][1]*ty + r[2][2]*tz + v[2];
        tx = x - targetatoms[i]->x();
        ty = y - targetatoms[i]->y();
        tz = z - targetatoms[i]->z();
        d = tx*tx + ty*ty + tz*tz;
        rms += d;
    }

    rms = sqrt(rms/(double)n);

    setSuccess(true);

    delete[] a;
    delete[] b;
    delete[] w;
    delete[] sourcefrom;
    delete[] targetatoms;
    delete[] sourceatoms;
}

void LSQFitDialog::on_removeButton_clicked()
{
    int isel = matchListBox->currentRow();
    if (isel >= 0)
    {
        QListWidgetItem *item = matchListBox->takeItem(isel);
        delete item;
        Matches.erase(Matches.begin() + isel);
    }
    setSuccess(false);
}

void LSQFitDialog::on_addButton_clicked()
{
    MATCH m;
    if (targetres == NULL)
    {
        targetres = m_target->getResidues();
    }
    if (sourceres == NULL)
    {
        sourceres = m_source->getResidues();
    }
    m.target = targetres;
    m.source = sourceres;
    m.length = lengthSpinBox->value();
    if (m.length == 0)
    {
        m.length = m_source->getnresidues();
    }
    m.atomtype = atomTypeText->text().toStdString();
    if (MatchOK(&m))
    {
        Matches.push_back(m);
    }
    ListMatches();
    setSuccess(false);
}

void LSQFitDialog::ListMatches()
{
    matchListBox->clear();
    std::vector<MATCH>::iterator m = Matches.begin();
    while (m != Matches.end())
    {
        std::string s = ::format("%s = %s :%d %s", resid((*m).target).c_str(), resid((*m).source).c_str(), (*m).length, (*m).atomtype.c_str());
        matchListBox->addItem(s.c_str());
        m++;
    }
}

void LSQFitDialog::ListSource()
{
    if (m_source)
    {
        RESIDUE *res = m_source->getResidues();
        QListWidget *list = sourceListWidget;
        QComboBox *chains = chainsChoice;
        list->clear();
        chains->clear();
        int chainid = -9921;
        while (res != NULL)
        {
            if (res->chain_id() != chainid)
            {
                chainid = res->chain_id();
                char cid = (char) chainid;
                chains->addItem(QString(cid));
            }
            list->addItem(resid(res).c_str());
            res = res->next();
        }
    }
}

void LSQFitDialog::ListTarget()
{
    if (m_target)
    {
        RESIDUE *res = m_target->getResidues();
        QListWidget *list = targetListWidget;
        list->clear();
        while (res != NULL)
        {
            list->addItem(resid(res).c_str());
            res = res->next();
        }
    }
}

void LSQFitDialog::ListSourceChoices()
{
    sourceComboBox->clear();
    int isel = 0, i = 0;
    std::list<Molecule*>::iterator p = displaylist->begin();
    while (p != displaylist->end())
    {
        if (m_source == *p)
        {
            isel = i;
        }
        sourceComboBox->addItem((*p)->pathname.c_str());
        p++;
        i++;
    }
    sourceComboBox->setCurrentIndex(isel);
}

void LSQFitDialog::ListTargetChoices()
{
    targetComboBox->clear();
    int isel = 0, i = 0;
    std::list<Molecule*>::iterator p = displaylist->begin();
    while (p != displaylist->end())
    {
        if (m_target == *p)
        {
            isel = i;
        }
        targetComboBox->addItem((*p)->pathname.c_str());
        p++;
        i++;
    }
    targetComboBox->setCurrentIndex(isel);
}

void LSQFitDialog::SetAtomType(const char *t)
{
    atomTypeText->setText(t);
}

bool LSQFitDialog::MatchOK(MATCH *m)
{
    RESIDUE *res = m->target;
    int n = 0;
    while ((res != NULL) &&  n < m->length)
    {
        if (!atom_from_name(m->atomtype.c_str(), *res))
        {
            break;
        }
        n++;
        res = res->next();
    }
    if (n < m->length)
    {
        m->length = n;
    }
    if (m->length <= 0)
    {
        return false;
    }

    res = m->source;
    n = 0;
    while ((res != NULL) && n < m->length)
    {
        if (!atom_from_name(m->atomtype.c_str(), *res))
        {
            break;
        }
        n++;
        res = res->next();
    }
    if (n < m->length)
    {
        m->length = n;
    }

    if (m->length <= 0)
    {
        return false;
    }

    return true;
}

void LSQFitDialog::on_exportButton_clicked()
{
    QString pathname = QFileDialog::getSaveFileName(this, "Save LSQMatrix (.mtx) File", "",
                                                 "Matrix files (*.mtx);;All files (*.*)");
    if (!pathname.isEmpty())
    {
        LSQMatrix lsq_matrix;
        lsq_matrix.SetMatrix(r, v);
        lsq_matrix.Save(pathname.toAscii().constData());
    }
}

void LSQFitDialog::on_importButton_clicked()
{
    QString pathname = QFileDialog::getOpenFileName(this, "Load LSQMatrix (.mtx) File", "",
                                                 "Matrix files (*.mtx);;All files (*.*)");
    if (!pathname.isEmpty())
    {
        LSQMatrix lsq_matrix;
        lsq_matrix.Load(pathname.toAscii().constData());
        lsq_matrix.GetMatrix(r, v);
        setSuccess(true);
    }
}


void LSQFitDialog::InitializeFromData(const MIData&)
{
    MIGLWidget *glw = MIMainWindow::instance()->currentMIGLWidget();
    if (!glw)
        return;

    displaylist = glw->GetDisplaylist();
    std::list<Molecule*>::iterator psource = displaylist->end();
    std::list<Molecule*>::iterator ptarget = displaylist->begin();
    m_target = *ptarget;
    m_source = *psource;

    ListSource();
    ListTarget();
    ListSourceChoices();
    ListTargetChoices();
    SetAtomType("CA");
}

void LSQFitDialog::GetData(MIData &dat)
{
    dat["sourceModel"].str = sourceComboBox->currentText().toStdString();
    dat["targetModel"].str = targetComboBox->currentText().toStdString();
    dat["applyToChain"].b = chainOnlyCheckBox->isChecked();
    dat["chainId"].str = chainsChoice->currentText().toStdString();

    dat["r00"].f = r[0][0];
    dat["r01"].f = r[0][1];
    dat["r02"].f = r[0][2];

    dat["r10"].f = r[1][0];
    dat["r11"].f = r[1][1];
    dat["r12"].f = r[1][2];

    dat["r20"].f = r[2][0];
    dat["r21"].f = r[2][1];
    dat["r22"].f = r[2][2];

    dat["v0"].f = v[0];
    dat["v1"].f = v[1];
    dat["v2"].f = v[2];
}

Molecule*LSQFitDialog::findMolecule(const std::string &str)
{
    for (std::list<Molecule*>::iterator p = displaylist->begin(); p != displaylist->end(); ++p)
    {
        if ((*p)->pathname.c_str() == str)
        {
            return *p;
        }
    }
    return 0;
}
