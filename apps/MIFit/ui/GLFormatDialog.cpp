#include "GLFormatDialog.h"
#include "ui_GLFormatDialog.h"
#include <QGLFormat>
#include <QFlags>
#include "GLFormatEdit.h"

namespace
{
const QString CURRENT_TYPE("Current");
const QString DEFAULT_TYPE("Default");
const QString DEFAULT_OVERLAY_TYPE("Default overlay");
}

GLFormatDialog::GLFormatDialog(QWidget *parent)
    : QDialog(parent),
      m_ui(new Ui::GLFormatDialog),
      currentFormat()
{
    m_ui->setupUi(this);

    // Hide initially and only show if format type set to Default.
    m_ui->changeDefault->hide();
    connect(m_ui->changeDefault, SIGNAL(clicked()),
            this, SLOT(changeDefault()));

    // Hide initially and only show if format type set to Default.
    m_ui->resetToQtDefault->hide();
    connect(m_ui->resetToQtDefault, SIGNAL(clicked()),
            this, SLOT(resetToQtDefault()));

    connect(m_ui->formatType, SIGNAL(currentIndexChanged(QString)),
            this, SLOT(formatTypeChanged(QString)));

    formatTypeChanged(DEFAULT_TYPE);
}

GLFormatDialog::~GLFormatDialog()
{
    delete m_ui;
}

void GLFormatDialog::setCurrentFormat(const QGLFormat &format)
{
    currentFormat = format;
    formatTypeChanged(CURRENT_TYPE);
}

namespace
{
typedef QPair<QString, QString> StringPair;

QString toTable(const QList<StringPair> &data)
{
    QString info("<table>");

    foreach (StringPair p, data)
    {
        info += QString("<tr><td align=right>%1:</td>"
                        "<td width=10>&nbsp;</td><td>%2</td></tr>")
                .arg(p.first).arg(p.second);
    }
    info += "</table>";
    return info;
}

QString formatToString(const QGLFormat &format, const QString &formatName)
{
    QList<StringPair> data;
    data += qMakePair<QString, QString>("Format", formatName);
    data += qMakePair<QString, QString>("Double buffer", format.doubleBuffer() ? "true" : "false");
    data += qMakePair<QString, QString>("Direct render", format.directRendering() ? "true" : "false");
    QStringList rgbSize;
    rgbSize << QString::number(format.redBufferSize())
            << QString::number(format.greenBufferSize())
            << QString::number(format.blueBufferSize());
    data += qMakePair<QString, QString>("RGB size", rgbSize.join(" "));
    data += qMakePair<QString, QString>("RGBA", format.rgba() ? "true" : "false");
    data += qMakePair<QString, QString>("Alpha", format.alpha() ? "true" : "false");
    data += qMakePair<QString, QString>("Alpha size", QString::number(format.alphaBufferSize()));
    data += qMakePair<QString, QString>("Depth", format.depth() ? "true" : "false");
    data += qMakePair<QString, QString>("Depth size", QString::number(format.depthBufferSize()));
    data += qMakePair<QString, QString>("Accum", format.accum() ? "true" : "false");
    data += qMakePair<QString, QString>("Accum size", QString::number(format.accumBufferSize()));
    data += qMakePair<QString, QString>("Stencil", format.stencil() ? "true" : "false");
    data += qMakePair<QString, QString>("Stencil size", QString::number(format.stencilBufferSize()));
    data += qMakePair<QString, QString>("Sample buffers", format.sampleBuffers() ? "true" : "false");
    data += qMakePair<QString, QString>("Samples", QString::number(format.samples()));
    data += qMakePair<QString, QString>("Overlay", format.hasOverlay() ? "true" : "false");
    data += qMakePair<QString, QString>("Stereo", format.stereo() ? "true" : "false");
    data += qMakePair<QString, QString>("Plane", QString::number(format.plane()));
    data += qMakePair<QString, QString>("Swap interval", QString::number(format.swapInterval()));

    return toTable(data);
}

}

void GLFormatDialog::formatTypeChanged(const QString &type)
{
    // Ensure format type UI set to type argument
    int typeIndex = m_ui->formatType->findText(type);
    if (typeIndex >= 0)
        m_ui->formatType->setCurrentIndex(typeIndex);

    QString info;
    QStringList flags;
    QGLFormat::OpenGLVersionFlags glFlags = QGLFormat::openGLVersionFlags();
    if (glFlags == QGLFormat::OpenGL_Version_None)
    {
        flags += "None";
    }
    else
    {
        if (glFlags.testFlag(QGLFormat::OpenGL_Version_3_0))
        {
            flags += "OpenGL_Version_3_0";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_2_1))
        {
            flags += "OpenGL_Version_2_1";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_2_0))
        {
            flags += "OpenGL_Version_2_0";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_1_5))
        {
            flags += "OpenGL_Version_1_5";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_1_4))
        {
            flags += "OpenGL_Version_1_4";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_1_3))
        {
            flags += "OpenGL_Version_1_3";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_1_2))
        {
            flags += "OpenGL_Version_1_2";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_Version_1_1))
        {
            flags += "OpenGL_Version_1_1";
        }
        if (glFlags.testFlag(QGLFormat::OpenGL_ES_Version_2_0))
        {
            flags += "OpenGL_ES_Version_2_0";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_ES_Common_Version_1_1))
        {
            flags += "OpenGL_ES_Common_Version_1_1";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_ES_Common_Version_1_0))
        {
            flags += "OpenGL_ES_Common_Version_1_0";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_ES_CommonLite_Version_1_1))
        {
            flags += "OpenGL_ES_CommonLite_Version_1_1";
        }
        else if (glFlags.testFlag(QGLFormat::OpenGL_ES_CommonLite_Version_1_0))
        {
            flags += "OpenGL_ES_CommonLite_Version_1_0";
        }
    }

    info += "<p>";
    info += flags.join(", ");
    info += "</p>\n";

    bool isDefault = type == DEFAULT_TYPE;
    m_ui->changeDefault->setVisible(isDefault);
    m_ui->resetToQtDefault->setVisible(isDefault);

    if (type == CURRENT_TYPE)
    {
        info += formatToString(currentFormat, type);
    }
    else if (type == DEFAULT_TYPE)
    {
        info += formatToString(QGLFormat::defaultFormat(), type);
    }
    else if (type == DEFAULT_OVERLAY_TYPE)
    {
        info += formatToString(QGLFormat::defaultOverlayFormat(), type);
    }
    m_ui->formatInfo->setText(info);
}

void GLFormatDialog::changeDefault()
{
    GLFormatEdit edit;
    edit.setFormat(QGLFormat::defaultFormat());
    if (edit.exec() == QDialog::Accepted)
    {
        QGLFormat::setDefaultFormat(edit.format());
        formatTypeChanged(DEFAULT_TYPE);
    }
}

void GLFormatDialog::resetToQtDefault()
{
    QGLFormat::setDefaultFormat(QGLFormat());
    formatTypeChanged(DEFAULT_TYPE);
}
