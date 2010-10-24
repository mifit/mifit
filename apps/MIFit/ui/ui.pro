include(../MIFit.pri)
TEMPLATE = lib
SOURCES = $$files(*.cpp)
HEADERS = $$files(*.h)
FORMS = $$files(*.ui)

RESOURCES += \
    ui.qrc

OTHER_FILES += \
    root.qml \
    stack.qml \
    message.qml
