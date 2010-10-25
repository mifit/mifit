#ifndef CursorArea_h
#define CursorArea_h

#include <QDeclarativeItem>

class CursorArea : public QDeclarativeItem
{
    Q_OBJECT
    Q_PROPERTY(QString shape READ cursor WRITE setCursor)

public:
    explicit CursorArea(QDeclarativeItem *parent = 0);

    QString cursor() const;
    void setCursor(const QString &shape);
};

#endif // CursorArea_h
