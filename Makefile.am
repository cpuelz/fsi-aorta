bin_PROGRAMS = main3d

main3ddir = $(top_builddir)

main3d_SOURCES = source/BoundaryConditions.cpp \
source/CirculationModel.cpp                    \
source/FeedbackForcer.cpp                      \
source/main.cpp                                \
source/MechanicsModel.cpp                      \
source/MeshInfo.cpp                            \
source/ModelInitialization.cpp                 \
source/Source.cpp                              \
source/SourceDistributer.cpp                   \
source/VelocityBcCoefs.cpp

main3d_HEADERS = include/four_chambered_heart/CirculationModel.h \
include/four_chambered_heart/BoundaryConditions.h                \
include/four_chambered_heart/ModelParameters.h                   \
include/four_chambered_heart/VelocityBcCoefs.h                   \
include/four_chambered_heart/MeshInfo.h                          \
include/four_chambered_heart/FeedbackForcer.h                    \
include/four_chambered_heart/MechanicsModel.h                    \
include/four_chambered_heart/ModelInitialization.h               \
include/four_chambered_heart/SourceDistributer.h                 \
include/four_chambered_heart/Source.h

main3d_CPPFLAGS = $(FCH_IBAMR_CPPFLAGS) -DNDIM=3
main3d_CXXFLAGS = $(FCH_IBAMR_CXXFLAGS) -I$(top_srcdir)/include/ -I$(top_builddir)/include/
main3d_LDFLAGS = $(FCH_IBAMR_LDFLAGS) $(FCH_IBAMR_LIB_3D) $(FCH_IBAMR_LIBS)
