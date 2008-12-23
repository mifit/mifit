#include <QPushButton>
#include <QLineEdit>
#include <QFileDialog>

#include "MIBrowsePair.h"

#include <string>

MIBrowsePair::MIBrowsePair(QPushButton *button, QLineEdit *lineEdit, const QString &filter, bool isDir) : 
  QObject(button), _button(button), _lineEdit(lineEdit),  _isDir(isDir), _filter(filter)
{
  connect(button, SIGNAL(clicked()), this, SLOT(buttonClicked()));
}

void MIBrowsePair::buttonClicked() {
  QString str;

  if (_isDir) {
    str=QFileDialog::getExistingDirectory(0, "", _lineEdit->text());
  } else {
    std::string foo=_filter.toStdString();
    foo=foo.substr(0,foo.find_first_of("("));
    str=QFileDialog::getOpenFileName(0, foo.c_str(), _lineEdit->text(), _filter);
  }

  if (str.size()) {
      _lineEdit->setText(str);
  }
}
