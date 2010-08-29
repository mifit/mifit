#include "BondRenderer.h"

BondRenderer::BondRenderer()
{
}

void BondRenderer::addBond(const QVector3D &pos1, const QVector4D &color1, const QVector3D &pos2, const QVector4D &color2)
{
    QVector3D midPoint = 0.5 * (pos1 + pos2);

    vertices.append(pos1);
    colors.append(color1);
    vertices.append(midPoint);
    colors.append(color1);
    vertices.append(midPoint);
    colors.append(color2);
    vertices.append(pos2);
    colors.append(color2);
}

void BondRenderer::draw() const
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, vertices.constData());
    glColorPointer(4, GL_FLOAT, 0, colors.constData());
    glDrawArrays(GL_LINES, 0, vertices.size());
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
}
