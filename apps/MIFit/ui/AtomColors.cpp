#include "AtomColors.h"

#include <QInputDialog>
#include "ui/uilib.h"
#include "ui/MIColorPickerDlg.h"

AtomColors::AtomColors(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);
}

void AtomColors::setColor(int c)
{
    MIPalette *m_palette = Application::instance()->GetLPpal();
    int pi = PaletteIndex(c);
    QColor color(m_palette->colors[pi].red,
                 m_palette->colors[pi].green,
                 m_palette->colors[pi].blue);

    // color the background with the color
    QPalette palette = colorToolButton->palette();
    palette.setColor(QPalette::Window, color);
    colorToolButton->setPalette(palette);

    // for mac, we have to set a solid-color pixmap too, b/c the
    // background palette color is ignored however, if we *just* use a
    // pixmap, there's no (straightforward) way to retreive the color.
    QPixmap p(32, 32);
    p.fill(color);
    QIcon icon(p);
    colorToolButton->setIcon(icon);
}

void AtomColors::populateList(const std::vector<std::string> &atomNames,
                              const std::vector<std::string> &atomColors)
{
    if (atomNames.size() != atomColors.size())
        return;

    listWidget->clear();
    for (unsigned int i = 0; i < atomNames.size(); ++i)
    {
        std::string foo = atomNames[i] + "* " + atomColors[i];
        listWidget->addItem(foo.c_str());
    }

    this->atomNames = atomNames;
    this->atomColors = atomColors;
    currentAtomNames = atomNames;
    currentAtomColors = atomColors;
}


void AtomColors::GetData(std::vector<std::string> &atomNames,
                         std::vector<std::string> &atomColors, bool &save)
{
    atomNames = currentAtomNames;
    atomColors = currentAtomColors;
    save = saveCheckBox->isChecked();
}

void AtomColors::on_deleteTypePushButton_clicked()
{
    if (listWidget->currentItem())
    {
        listWidget->currentRow();
        currentAtomNames.erase(currentAtomNames.begin()
                                               +listWidget->currentRow());
        currentAtomColors.erase(currentAtomColors.begin()
                                                +listWidget->currentRow());
        delete listWidget->takeItem(listWidget->currentRow());
    }
}

void AtomColors::on_addTypePushButton_clicked()
{
    bool ok;
    QString foo = QInputDialog::getText(this, "Atom Type Name", "Enter atom type name", QLineEdit::Normal, "", &ok);
    if (!ok)
        return;
    int ci = MIColorPickerDlg::getColor(0, 1, "Choose atom color");
    foo += "* ";
    foo += Colors::colornames[ci];
    listWidget->addItem(foo);

    currentAtomNames.push_back(foo.toStdString());
    currentAtomColors.push_back(Colors::colornames[ci]);
}


void AtomColors::on_resetPushButton_clicked()
{
    populateList(atomNames, atomColors);
}

void AtomColors::on_editNamePushButton_clicked()
{
    if (!listWidget->currentItem())
        return;
    bool ok;
    QString foo = QInputDialog::getText(this, "Atom Type Name", "Enter atom type name", QLineEdit::Normal,
                                        currentAtomNames[listWidget->currentRow()].c_str(), &ok);
    if (!ok)
        return;

    currentAtomNames[listWidget->currentRow()] = foo.toStdString();
    foo += "* ";
    foo += currentAtomNames[listWidget->currentRow()].c_str();
    listWidget->currentItem()->setText(foo);
}

void AtomColors::on_colorToolButton_clicked()
{
    if (!listWidget->currentItem())
        return;

    int ci = MIColorPickerDlg::getColor(0, Colors::findColorNumber(currentAtomColors[listWidget->currentRow()]), "Choose atom color");
    std::string foo = currentAtomNames[listWidget->currentRow()];
    foo += "* ";
    foo += Colors::colornames[ci];
    currentAtomColors[listWidget->currentRow()] = Colors::colornames[ci];
    listWidget->currentItem()->setText(foo.c_str());

    setColor(ci);
}

void AtomColors::on_listWidget_currentRowChanged(int idx)
{
    if (idx < 0 || idx > (int)currentAtomColors.size())
        return;
    std::string name = currentAtomColors[idx];
    int i = Colors::findColorNumber(name);
    if (i!=-1)
    {
        setColor(i);
    }
}
