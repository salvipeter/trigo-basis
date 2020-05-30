# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp trigo-basis.hh
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp trigo-basis.cc

INCLUDEPATH += /usr/include/eigen3
LIBS *= -lQGLViewer-qt5 -L/usr/lib/OpenMesh -lOpenMeshCore -lGL -lGLU

RESOURCES = sample-framework.qrc
