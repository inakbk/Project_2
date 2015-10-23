TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES +=

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

LIBS += -larmadillo -llapack -lblas

HEADERS += \
    ../project_2_test_jacobi/jacobi.h \
    ../project_2_test_jacobi/jacobisolver.h \
    ../project_2_test_jacobi/writetofile.h
