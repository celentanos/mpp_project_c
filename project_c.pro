TEMPLATE = app
CONFIG += c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += terminal

TARGET = project_c

DESTDIR = $$PWD

# MPI Settings #################################################################
QMAKE_CXX = mpicxx
QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_CC = mpicc

#QMAKE_CFLAGS = $$system(mpicc --showme:compile)
#QMAKE_LFLAGS = $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
################################################################################

SOURCES += \
    main.cpp

unix:!macx: LIBS += -L/usr/lib/openmpi/lib/ -lmpi

INCLUDEPATH += /usr/lib/openmpi/include
DEPENDPATH += /usr/lib/openmpi/include
