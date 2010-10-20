include(../libs.pri)

TEMPLATE = lib
TARGET = miopengl
DEPENDPATH += interact zlib
INCLUDEPATH += zlib 

# Input
HEADERS += Arc.h \
           Axes.h \
           Axis.h \
           Camera.h \
           Circle.h \
           Frustum.h \
           Light.h \
           Object3d.h \
           OpenGL.h \
           QuatUtil.h \
           RelativeViewpoint.h \
           Renderable.h \
           Scene.h \
           Sphere.h \
           StereoView.h \
           Viewpoint.h \
           Viewport.h \
           ViewportRelativeViewpoint.h \
           interact/ArcBallFeedback.h \
           interact/FieldOfViewZoomCommand.h \
           interact/MouseArcBall.h \
           interact/MouseArcBallOrbitor.h \
           interact/MousePicker.h \
           interact/MouseTranslator.h \
           interact/MouseZoomer.h \
           interact/PropertyCommand.h \
           interact/SimpleMouseRotator.h \
           interact/TargetFeedback.h \
           interact/TranslationFeedback.h \
           zlib/infblock.h \
           zlib/infcodes.h \
           zlib/inffast.h \
           zlib/inffixed.h \
           zlib/inftrees.h \
           zlib/infutil.h \
           zlib/zconf.h \
           zlib/zlib.h \
           zlib/zutil.h \
           Text.h
SOURCES += Arc.cpp \
           Axes.cpp \
           Axis.cpp \
           Camera.cpp \
           Circle.cpp \
           Frustum.cpp \
           Light.cpp \
           Object3d.cpp \
           QuatUtil.cpp \
           RelativeViewpoint.cpp \
           Sphere.cpp \
           StereoView.cpp \
           Viewpoint.cpp \
           Viewport.cpp \
           ViewportRelativeViewpoint.cpp \
           interact/ArcBallFeedback.cpp \
           interact/FieldOfViewZoomCommand.cpp \
           interact/MouseArcBall.cpp \
           interact/MouseArcBallOrbitor.cpp \
           interact/MousePicker.cpp \
           interact/MouseTranslator.cpp \
           interact/MouseZoomer.cpp \
           interact/SimpleMouseRotator.cpp \
           interact/TargetFeedback.cpp \
           interact/TranslationFeedback.cpp \
           zlib/adler32.cpp \
           zlib/crc32.cpp \
           zlib/infblock.cpp \
           zlib/infcodes.cpp \
           zlib/inffast.cpp \
           zlib/inflate.cpp \
           zlib/inftrees.cpp \
           zlib/infutil.cpp \
           zlib/zutil.cpp \
           Text.cpp
