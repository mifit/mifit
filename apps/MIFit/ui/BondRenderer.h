#ifndef BONDRENDERER_H
#define BONDRENDERER_H

#include <QtGui/QVector3D>
#include <QtGui/QVector4D>
#include <QtOpenGL/qgl.h>

class BondRenderer
{
public:
    BondRenderer();

    void loadArrays() const;
    void addBond(const QVector3D &pos1, const QVector4D &color1, const QVector3D &pos2, const QVector4D &color2);
    void draw() const;

private:
    QVector<QVector3D> vertices;
    QVector<QVector4D> colors;
};

#endif // BONDRENDERER_H
