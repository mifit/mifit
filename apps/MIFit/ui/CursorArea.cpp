#include "CursorArea.h"

#include <QtGui/QCursor>

namespace
{
    QMap<Qt::CursorShape, QString> cursorShapeStrings;
    QMap<QString, Qt::CursorShape> cursorShapeMap;

    void initCursorShapeMaps()
    {
        static bool initialized = false;
        if (initialized)
            return;
        initialized = true;

        cursorShapeStrings[Qt::ArrowCursor] = "arrow";
        cursorShapeStrings[Qt::UpArrowCursor] = "upArrow";
        cursorShapeStrings[Qt::CrossCursor] = "cross";
        cursorShapeStrings[Qt::WaitCursor] = "wait";
        cursorShapeStrings[Qt::IBeamCursor] = "iBeam";
        cursorShapeStrings[Qt::SizeVerCursor] = "sizeVer";
        cursorShapeStrings[Qt::SizeHorCursor] = "sizeHor";
        cursorShapeStrings[Qt::SizeBDiagCursor] = "sizeBDiag";
        cursorShapeStrings[Qt::SizeFDiagCursor] = "sizeFDiag";
        cursorShapeStrings[Qt::SizeAllCursor] = "sizeAll";
        cursorShapeStrings[Qt::BlankCursor] = "blank";
        cursorShapeStrings[Qt::SplitVCursor] = "splitV";
        cursorShapeStrings[Qt::SplitHCursor] = "splitH";
        cursorShapeStrings[Qt::PointingHandCursor] = "pointingHand";
        cursorShapeStrings[Qt::ForbiddenCursor] = "forbidden";
        cursorShapeStrings[Qt::OpenHandCursor] = "openHand";
        cursorShapeStrings[Qt::ClosedHandCursor] = "closedHand";
        cursorShapeStrings[Qt::WhatsThisCursor] = "whatsThis";
        cursorShapeStrings[Qt::BusyCursor] = "busy";
        cursorShapeStrings[Qt::DragMoveCursor] = "dragMove";
        cursorShapeStrings[Qt::DragCopyCursor] = "dragCopy";
        cursorShapeStrings[Qt::DragLinkCursor] = "dragLink";
        cursorShapeStrings[Qt::BitmapCursor] = "bitmap";

        QMap<Qt::CursorShape, QString>::const_iterator i = cursorShapeStrings.constBegin();
        for (; i != cursorShapeStrings.constEnd(); ++i) {
            cursorShapeMap[i.value().toLower()] = i.key();
        }
    }
};

CursorArea::CursorArea(QDeclarativeItem *parent)
    : QDeclarativeItem(parent)
{
    initCursorShapeMaps();
}

QString CursorArea::cursor() const
{
    return cursorShapeStrings[QGraphicsItem::cursor().shape()];
}

void CursorArea::setCursor(const QString &shape)
{
    if (cursorShapeMap.contains(shape.toLower()))
        QGraphicsItem::setCursor(cursorShapeMap[shape.toLower()]);
}
