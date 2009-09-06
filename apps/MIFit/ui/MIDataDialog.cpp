#include "MIDialog.h"
#include "MIDataDialog.h"

#include <QApplication>
#include <QGridLayout>
#include <QLabel>
#include <QSpinBox>
#include <QLineEdit>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QToolButton>
#include <QColorDialog>
#include <QComboBox>

#include <nongui/nongui.h>
#include <util/utillib.h>
#include "core/corelib.h"
#include "ui/uilib.h"

MIDataDialog::MIDataDialog(QWidget* parent, Qt::WindowFlags f)
: QDialog(parent, f), data(NULL) {
  
  layout = new QGridLayout(this);
  setLayout(layout);
}

static void setButtonColor(QToolButton *valueControl, const QColor &color)
{
  // color the background with the color
  QPalette palette=valueControl->palette();
  palette.setColor(QPalette::Window, color);
  valueControl->setPalette(palette);
  
  // for mac, we have to set a solid-color pixmap too, b/c the
  // background palette color is ignored however, if we *just* use a
  // pixmap, there's no (straightforward) way to retreive the color.
  QPixmap p(32,32);
  p.fill(color);
  QIcon icon(p);
  valueControl->setIcon(icon);
}



void MIDataDialog::addControl(const std::string &key, const MIDatum &value, int row)
{
  QLabel* label = new QLabel(key.c_str());
  if (labels.find(key) != labels.end()) {
    label->setText(labels[key].c_str());
  }
  layout->addWidget(label, row, 0, Qt::AlignRight);

  QWidget* valueControl;
  ControlType valueControlType;
  if (value.str != MIDatum::INVALID_STRING) {
    valueControl = new QLineEdit(value.str.c_str());
    valueControlType = TEXT;
  } else if (value.radio != UINT_MAX) {
    if (value.radio_count != value.radio_labels.size()) {
      Logger::debug("Programmer error: radio_labels size != radio_count!");
      return;
    }
    valueControl = new QComboBox();
    valueControlType = RADIO;
    QStringList qsl;
    for (size_t i=0; i<value.radio_labels.size(); ++i) {
      qsl.append(value.radio_labels[i].c_str());
    }
    ((QComboBox*)valueControl)->addItems(qsl);
    ((QComboBox*)valueControl)->setCurrentIndex((int)value.radio);
  } else if (value.i != INT_MIN) {
    std::string str = format("%d", value.i);
    valueControl = new QLineEdit(str.c_str());
    valueControlType = INT;
  } else if (value.u != UINT_MAX) {
    std::string str = format("%u", value.u);
    valueControl = new QLineEdit(str.c_str());
    valueControlType = UNSIGNEDINT;
  } else if (value.s != SHRT_MIN) {
    std::string str = format("%d", value.s);
    valueControl = new QLineEdit(str.c_str());
    valueControlType = SHORT;
  } else if (value.f != FLT_MIN) {
    std::string str = format("%f", value.f);
    valueControl = new QLineEdit(str.c_str());
    valueControlType = FLOAT;
  } else if (value.d != DBL_MIN) {
    std::string str = format("%f", value.d);
    valueControl = new QLineEdit(str.c_str());
    valueControlType = DOUBLE;
  } else if (value.isColor) {
    valueControl = new QToolButton();
    valueControl->setAutoFillBackground(true);
    QColor color(value.color[0],value.color[1],value.color[2]);
    setButtonColor((QToolButton*)valueControl, color);
    valueControlType = COLOR;
    connect(valueControl, SIGNAL(clicked()), this, SLOT(colorButtonPressed()));
  } else if (value.isColorIndex) {
    valueControl = new QToolButton();
    colorIndexMap[valueControl] = value.color[0];
    valueControl->setAutoFillBackground(true);

    MIPalette *m_palette=Application::instance()->GetLPpal();
    int ci = PaletteIndex(value.color[0]);
    QColor color(m_palette->colors[ci].red,
                 m_palette->colors[ci].green,
                 m_palette->colors[ci].blue);
    setButtonColor((QToolButton*)valueControl, color);
    valueControlType = COLORINDEX;
    connect(valueControl, SIGNAL(clicked()), this, SLOT(colorIndexButtonPressed()));
  } else {
    QCheckBox* checkbox = new QCheckBox("");
    checkbox->setCheckState(value.b ? Qt::Checked : Qt::Unchecked);
    valueControl = checkbox;
    valueControlType = BOOLEAN;
  }
  layout->addWidget(valueControl, row, 1);
  dataControls[key] = std::make_pair(valueControlType, valueControl);
}

void MIDataDialog::setMIData(MIData* data) {
  this->data = data;
  dataControls.clear();

  QLayoutItem *child;
  while ((child = layout->takeAt(0)) != 0) {
    delete child;
  }
  
  int row = 0;
  if (order_.size()==0) {
    for (MIData::const_iterator iter = data->begin(); iter != data->end(); ++iter) {
      std::string key = iter->first;
      MIDatum value = iter->second;
      addControl(key,value,row++);
    }
  } else {
    for (unsigned int i=0; i< order_.size(); ++i) {
      if (data->find(order_[i]) != data->end()) {
        addControl(order_[i],(*data)[order_[i]],row++);
      }
    }
  }

  QDialogButtonBox* buttons = new QDialogButtonBox;
  buttons->addButton(QDialogButtonBox::Ok);
  buttons->addButton(QDialogButtonBox::Cancel);
  layout->addWidget(buttons, row, 0, 1, -1, Qt::AlignRight);
  QObject::connect(buttons, SIGNAL(accepted()), this, SLOT(accepted()));
  QObject::connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));

}

void MIDataDialog::accepted() {
  DataControlMap::iterator iter;
  for (iter = dataControls.begin(); iter != dataControls.end(); ++iter) {
    const std::string& key = iter->first;
    ControlType type = iter->second.first;
    QWidget* control = iter->second.second;
    switch (type) {
    case TEXT:
      {
        QLineEdit* text = (QLineEdit*) control;
        (*data)[key].str = text->text().toLatin1().data();
      }
      break;
    case RADIO:
      {
        QComboBox* combo = (QComboBox*) control;
        (*data)[key].radio=(unsigned int)combo->currentIndex();
      }
      break;
    case INT:
      {
        QLineEdit* text = (QLineEdit*) control;
        int i;
        int result = sscanf(text->text().toLatin1().data(), "%d", &i);
        if (result == 1) {
          (*data)[key].i = i;
        }
      }
      break;
    case UNSIGNEDINT:
      {
        QLineEdit* text = (QLineEdit*) control;
        unsigned int u;
        int result = sscanf(text->text().toLatin1().data(), "%u", &u);
        if (result == 1) {
          (*data)[key].u = u;
        }
      }
      break;
    case SHORT:
      {
        QLineEdit* text = (QLineEdit*) control;
        int i;
        int result = sscanf(text->text().toLatin1().data(), "%d", &i);
        if (result == 1) {
          (*data)[key].s = (short)i;
        }
      }
      break;
    case FLOAT:
      {
        QLineEdit* text = (QLineEdit*) control;
        float f;
        int result = sscanf(text->text().toLatin1().data(), "%f", &f);
        if (result == 1) {
          (*data)[key].f = f;
        }
      }
      break;
    case DOUBLE:
      {
        QLineEdit* text = (QLineEdit*) control;
        float f;
        int result = sscanf(text->text().toLatin1().data(), "%f", &f);
        if (result == 1) {
          (*data)[key].d = f;
        }
      }
      break;
    case BOOLEAN:
      {
        QCheckBox* checkbox = (QCheckBox*) control;
        (*data)[key].b = checkbox->checkState() != Qt::Unchecked;
      }
      break;
    case COLOR:
      {
        QToolButton *button=(QToolButton*) control;
        QColor color=button->palette().color(QPalette::Normal, QPalette::Window);
        (*data)[key].isColor = true;
        (*data)[key].color[0] = (unsigned char)color.red();
        (*data)[key].color[1] = (unsigned char)color.green();
        (*data)[key].color[2] = (unsigned char)color.blue();
      }
      break;
    case COLORINDEX:
      {
        QToolButton *button=(QToolButton*) control;
        (*data)[key].isColorIndex=true;
        (*data)[key].color[0]=colorIndexMap[button];
      }
      break;
    }
  }
  accept();
}


void MIDataDialog::colorIndexButtonPressed() {
  //find relevant control
  QWidget *control=(QToolButton*)sender();
  colorIndexMap[control]=(unsigned char)MIColorChooser(colorIndexMap[control]);

  MIPalette *m_palette=Application::instance()->GetLPpal();
  int ci = PaletteIndex(colorIndexMap[control]);
  QColor color(m_palette->colors[ci].red,
               m_palette->colors[ci].green,
               m_palette->colors[ci].blue);

  setButtonColor((QToolButton*)control, color);
}

void MIDataDialog::colorButtonPressed() {
  //find relevant control
  QWidget *w=(QWidget*)sender();

  DataControlMap::iterator iter;
  for (iter = dataControls.begin(); iter != dataControls.end(); ++iter) {
    QWidget* control = iter->second.second;

    if (control == w) {
      QColor initColor=control->palette().color(QPalette::Normal, QPalette::Window);
      QColor color=QColorDialog::getColor(initColor, this);
      if (color.isValid()) {
        setButtonColor((QToolButton*)control, color);
      }
      break;
    }
  }
}

void MIDataDialog::order(const std::string& key) {
  order_.push_back(key);
}

void MIDataDialog::label(const std::string& key, const std::string& label) {
  labels[key] = label;
}
