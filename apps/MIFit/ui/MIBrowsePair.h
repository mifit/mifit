#ifndef MI_BROWSE_PAIR_H
#define MI_BROWSE_PAIR_H

class QPushButton;
class QLineEdit;

#include <string>

#include <QObject>
#include <QString>

class MIBrowsePair : public QObject
{
    Q_OBJECT
  public:
    MIBrowsePair(QPushButton *button, QLineEdit *lineEdit, 
      const QString & filter="",bool isDir=false);

  private slots:
    void buttonClicked();

  private:
    // these are *not* owned by this class
    QPushButton *_button;
    QLineEdit *_lineEdit;
    bool _isDir;
    QString _filter;
};

#endif
