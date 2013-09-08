AM_CPPFLAGS = -I$(top_srcdir)/include/

bin_PROGRAMS = \
	chomp-cubical \
	chomp-matrix \
	chomp-simplicial \
	greyscale-to-cubical


chompincludedir = $(includedir)/chomp
chompinclude_HEADERS = \
	../include/chomp/BitmapSubcomplex.h       ../include/chomp/MorseComplex.h          \
	../include/chomp/Chain.h                  ../include/chomp/PolyRing.h              \
	../include/chomp/Closure.h                ../include/chomp/Preboundary.h           \
	../include/chomp/Complex.h                ../include/chomp/Prism.h                 \
	../include/chomp/ConleyIndex.h            ../include/chomp/Rect.h                  \
	../include/chomp/CoreductionDecomposer.h  ../include/chomp/RelativeMapHomology.h   \
	../include/chomp/CubicalComplex.h         ../include/chomp/RelativePair.h          \
	../include/chomp/Decomposer.h             ../include/chomp/Ring.h                  \
	../include/chomp/Draw.h                   ../include/chomp/SimplicialComplex.h     \
	../include/chomp/FiberComplex.h           ../include/chomp/SmithNormalForm.h       \
	../include/chomp/Field.h                  ../include/chomp/SparseMatrix.h          \
	../include/chomp/Generators.h             ../include/chomp/Subcomplex.h            \
	../include/chomp/GraphComplex.h           ../include/chomp/ToplexDetail.h          \
	../include/chomp/HashIndexer.h            ../include/chomp/Toplex.h                \
	../include/chomp/MatrixComplex.h          ../include/chomp/VectorIndexer.h         \
	../include/chomp/Matrix.h


# TODO: I think these are now part of boost, but in differently-named files
boostincludedir = $(includedir)/chomp/boost
boostinclude_HEADERS = \
	../include/boost/serialization/unordered_set.hpp \
	../include/boost/serialization/unordered_map.hpp \
	../include/boost/serialization/unordered_collections_load_imp.hpp \
	../include/boost/serialization/unordered_collections_save_imp.hpp


chomp_cubical_SOURCES = \
	chomp-cubical.cpp


chomp_matrix_SOURCES = \
	chomp-matrix.cpp


chomp_simplicial_SOURCES = \
	chomp-simplicial.cpp


greyscale_to_cubical_SOURCES = \
	greyscale-to-cubical.cpp \
	../include/CImg.h

greyscale_to_cubical_LDADD = \
	-lm -lpthread -lX11

