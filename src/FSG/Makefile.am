noinst_LIBRARIES = libfsg.a

libfsg_a_CXXFLAGS = -std=c++11

libfsg_a_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/StringGraph \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/SGA

libfsg_a_SOURCES = \
	parallel_util.hpp \
	compressor.hpp compressor.cpp \
	vector_locker.hpp \
	simple_bitvector.hpp\
	fsg.h fsg.cpp \
	bwt_interval.hpp \
	suff_mgr.hpp \
	ov_mgr.hpp ov_mgr.cpp \
	segment.hpp \
	FSGCommon.h \
	hits.hpp \
	hits_util.hpp hits_util.cpp \
	delta_interval.hpp
