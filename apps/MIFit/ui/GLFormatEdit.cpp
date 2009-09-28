#include "GLFormatEdit.h"
#include "ui_GLFormatEdit.h"
#include <QPushButton>
#include <QSettings>

GLFormatEdit::GLFormatEdit(QWidget *parent) :
    QDialog(parent),
    m_ui(new Ui::GLFormatEdit)
{
    m_ui->setupUi(this);

    connect(m_ui->buttonBox, SIGNAL(clicked(QAbstractButton*)),
            this, SLOT(buttonClicked(QAbstractButton*)));
}

GLFormatEdit::~GLFormatEdit()
{
    delete m_ui;
}

QGLFormat GLFormatEdit::format() const
{
    QGLFormat result;
    result.setDoubleBuffer(m_ui->doubleBuffer->isChecked());
    result.setDirectRendering(m_ui->direct->isChecked());
    int rgbSize = m_ui->rgbSize->value();
    if (rgbSize >= 0) {
        result.setRedBufferSize(rgbSize);
        result.setGreenBufferSize(rgbSize);
        result.setBlueBufferSize(rgbSize);
    }
    result.setRgba(m_ui->rgba->isChecked());
    result.setAlpha(m_ui->alpha->isChecked());
    int alphaSize = m_ui->alphaSize->value();
    if (alphaSize >= 0)
        result.setAlphaBufferSize(alphaSize);
    result.setDepth(m_ui->depth->isChecked());
    int depthSize = m_ui->depthSize->value();
    if (depthSize >= 0)
        result.setDepthBufferSize(depthSize);
    result.setAccum(m_ui->accum->isChecked());
    int accumSize = m_ui->accumSize->value();
    if (accumSize >= 0)
        result.setAccumBufferSize(accumSize);
    result.setStencil(m_ui->stencil->isChecked());
    int stencilSize = m_ui->stencilSize->value();
    if (stencilSize >= 0)
        result.setStencilBufferSize(stencilSize);
    result.setSampleBuffers(m_ui->sampleBuffers->isChecked());
    int sampleBuffersSize = m_ui->sampleBuffersSize->value();
    if (sampleBuffersSize >= 0)
        result.setSamples(sampleBuffersSize);
    result.setStereo(m_ui->stereo->isChecked());
    result.setPlane(m_ui->plane->value());
    int swapInterval = m_ui->swapInterval->value();
    if (swapInterval >= 0)
        result.setSwapInterval(swapInterval);
    return result;
}

void GLFormatEdit::setFormat(const QGLFormat& format)
{
    format_ = format;
    reset();
}

void GLFormatEdit::reset()
{
    m_ui->direct->setChecked(format_.directRendering());
    m_ui->doubleBuffer->setChecked(format_.doubleBuffer());
    int rgbSize = qMin(format_.redBufferSize(),
                       qMin(format_.greenBufferSize(), format_.blueBufferSize()));
    m_ui->rgbSize->setValue(rgbSize);
    m_ui->rgba->setChecked(format_.rgba());
    m_ui->alpha->setChecked(format_.alpha());
    m_ui->alphaSize->setValue(format_.alphaBufferSize());
    m_ui->depth->setChecked(format_.depth());
    m_ui->depthSize->setValue(format_.depthBufferSize());
    m_ui->accum->setChecked(format_.accum());
    m_ui->accumSize->setValue(format_.accumBufferSize());
    m_ui->stencil->setChecked(format_.stencil());
    m_ui->stencilSize->setValue(format_.stencilBufferSize());
    m_ui->sampleBuffers->setChecked(format_.sampleBuffers());
    m_ui->sampleBuffersSize->setValue(format_.samples());
    m_ui->stereo->setChecked(format_.stereo());
    m_ui->plane->setValue(format_.plane());
    m_ui->swapInterval->setValue(format_.swapInterval());
}

void GLFormatEdit::buttonClicked(QAbstractButton* button)
{
    if (button == m_ui->buttonBox->button(QDialogButtonBox::Reset))
        reset();
}

namespace {
    const QString DOUBLE_BUFFER("doubleBuffer");
    const QString DIRECT_RENDER("directRender");
    const QString RGB_SIZE("rgbSize");
    const QString RGBA("rgba");
    const QString ALPHA("alpha");
    const QString ALPHA_SIZE("alphaSize");
    const QString DEPTH("depth");
    const QString DEPTH_SIZE("depthSize");
    const QString ACCUM("accum");
    const QString ACCUM_SIZE("accumSize");
    const QString STENCIL("stencil");
    const QString STENCIL_SIZE("stencilSize");
    const QString SAMPLE_BUFFERS("sampleBuffers");
    const QString SAMPLE_BUFFERS_SIZE("sampleBuffersSize");
    const QString STEREO("stereo");
    const QString PLANE("plane");
    const QString SWAP_INTERVAL("swapInterval");
}

QGLFormat GLFormatEdit::readSettings(const QSettings& settings)
{
    QGLFormat result;

    result.setDoubleBuffer(settings.value(DOUBLE_BUFFER).toBool());
    result.setDirectRendering(settings.value(DIRECT_RENDER).toBool());

    bool intOk;
    int rgbSize = settings.value(RGB_SIZE).toInt(&intOk);
    if (intOk && rgbSize >= 0) {
        result.setRedBufferSize(rgbSize);
        result.setGreenBufferSize(rgbSize);
        result.setBlueBufferSize(rgbSize);
    }

    result.setRgba(settings.value(RGBA).toBool());

    result.setAlpha(settings.value(ALPHA).toBool());
    int alphaSize = settings.value(ALPHA_SIZE).toInt(&intOk);
    if (intOk && alphaSize >= 0)
        result.setAlphaBufferSize(alphaSize);

    result.setDepth(settings.value(DEPTH).toBool());
    int depthSize = settings.value(DEPTH_SIZE).toInt(&intOk);
    if (intOk && depthSize >= 0)
        result.setDepthBufferSize(depthSize);

    result.setAccum(settings.value(ACCUM).toBool());
    int accumSize = settings.value(ACCUM_SIZE).toInt(&intOk);
    if (intOk && accumSize >= 0)
        result.setAccumBufferSize(accumSize);

    result.setStencil(settings.value(STENCIL).toBool());
    int stencilSize = settings.value(STENCIL_SIZE).toInt(&intOk);
    if (intOk && stencilSize >= 0)
        result.setStencilBufferSize(stencilSize);

    result.setSampleBuffers(settings.value(SAMPLE_BUFFERS).toBool());
    int sampleBuffersSize = settings.value(SAMPLE_BUFFERS_SIZE).toInt(&intOk);
    if (intOk && sampleBuffersSize >= 0)
        result.setSamples(sampleBuffersSize);

    result.setStereo(settings.value(STEREO).toBool());

    int plane = settings.value(PLANE).toInt(&intOk);
    if (intOk)
        result.setPlane(plane);

    int swapInterval = settings.value(SWAP_INTERVAL).toInt(&intOk);
    if (intOk && swapInterval >= 0)
        result.setSwapInterval(swapInterval);

    return result;
}

void GLFormatEdit::writeSettings(QSettings& settings, const QGLFormat& format)
{
    settings.setValue(DOUBLE_BUFFER, format.doubleBuffer());
    settings.setValue(DIRECT_RENDER, format.directRendering());
    int rgbSize = qMin(format.redBufferSize(),
                       qMin(format.greenBufferSize(), format.blueBufferSize()));
    settings.setValue(RGB_SIZE, rgbSize);
    settings.setValue(RGBA, format.rgba());
    settings.setValue(ALPHA, format.alpha());
    settings.setValue(ALPHA_SIZE, format.alphaBufferSize());
    settings.setValue(DEPTH, format.depth());
    settings.setValue(DEPTH_SIZE, format.depthBufferSize());
    settings.setValue(ACCUM, format.accum());
    settings.setValue(ACCUM_SIZE, format.accumBufferSize());
    settings.setValue(STENCIL, format.stencil());
    settings.setValue(STENCIL_SIZE, format.stencilBufferSize());
    settings.setValue(SAMPLE_BUFFERS, format.sampleBuffers());
    settings.setValue(SAMPLE_BUFFERS_SIZE, format.samples());
    settings.setValue(STEREO, format.stereo());
    settings.setValue(PLANE, format.plane());
    settings.setValue(SWAP_INTERVAL, format.swapInterval());
}

