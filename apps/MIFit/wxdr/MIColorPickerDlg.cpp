#include "MIColorPickerDlg.h"
#include "core/corelib.h"
#include "ui/uilib.h"

#include <QButtonGroup>
#include <QPushButton>
#include <QToolButton>
#include <QLayout>
#include <QGridLayout>

MIColorPickerDlg::MIColorPickerDlg(QWidget *parent, 
                                   int selected) :
  QDialog(parent),result(selected) {
  setWindowTitle("Choose color");
  QColor colorTable[Colors_NUMBERCOLORS];
  for (int i=1;i<Colors_NUMBERCOLORS;++i) {
    MIPalette *m_palette=Application::instance()->GetLPpal();
    int color = PaletteIndex(i);
    colorTable[i]=QColor(m_palette->colors[color].red,
                         m_palette->colors[color].green,
                         m_palette->colors[color].blue);
  }

  QButtonGroup *bg=new QButtonGroup(this);
  bg->setExclusive(true);
  QGridLayout *layout=new QGridLayout(this);
  layout->setVerticalSpacing(1);
  layout->setHorizontalSpacing(1);

  for (int i=1;i<Colors_NUMBERCOLORS; ++i) {
    QToolButton *pb=new MIColorToolButton(colorTable[i],this);
    bg->addButton(pb,i);
    pb->setMaximumSize(25,25);
    pb->setMinimumSize(25,25);
    pb->setCheckable(true);
    pb->setChecked(i==selected);
    layout->addWidget(pb, (i-1)/10, (i-1)%10);
  }

  connect(bg, SIGNAL( buttonClicked(int) ), this, SLOT (colorPicked(int)));
}

void MIColorPickerDlg::colorPicked(int i) {
  result=i;
  accept();
}
