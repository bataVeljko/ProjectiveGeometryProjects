TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lmgl

SOURCES += \
        main.cpp

unix: CONFIG += link_pkgconfig
