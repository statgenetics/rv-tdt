# makefile for rareTdt project #

VERSION:= 2.0.0
SRCDIR:= ./src/
CC:= gcc		# The C compiler.
CFLAGS:= -O3 -Wall
			# C compilation options which relate to
			# optimization or debugging (usually
    			# just -g or -O).  Usually this wouldn't
    			# include -I options to specify the
    			# include directories, because then you
    			# couldn't override it on the command line
    			# easily as in the above example.
#CXX:= icpc
CXXPG:= g++ -pg
CXX:= g++ -g
			# The C++ compiler.
			# (Sometimes "CPP" instead of CXX.)
#CXXFLAGS:= -O3 -Wall -wd981	-wd1572 -wd383 -wd1418
CXXFLAGS:= -O3 -Wall 
			# C++ compilation options related to
    			# optimization or debugging (-O or -g).

INCLUDEDIR:= -I/usr/include -I./include/ -I./src/ -I./src/gsl -I./src/src_old -I./src/fb-skat -I./src/test
LIBDIR:= -L/usr/lib -L./libs/
			# Contains libraries we need to link in.
LIBS:= -lm  #-lgsl -lgslcblas # -ljson_linux-gcc-4.5.2_libmt

JSON:= src/lib_json/json_reader.cpp src/lib_json/json_value.cpp  src/lib_json/json_writer.cpp
JOBJS:= $(JSON:.cpp=.o)

GSL:=	src/gsl/error.c \
   src/gsl/sys/infnan.c \
   src/gsl/sys/coerce.c \
   src/gsl/sys/fdiv.c \
   src/gsl/sys/pow_int.c \
   src/gsl/sys/log1p.c \
   src/gsl/sys/fcmp.c \
   src/gsl/sys/invhyp.c \
   src/gsl/complex/math.c \
   src/gsl/specfunc/beta.c \
   src/gsl/specfunc/psi.c \
   src/gsl/specfunc/trig.c \
   src/gsl/specfunc/exp.c \
   src/gsl/specfunc/expint.c \
   src/gsl/specfunc/log.c \
   src/gsl/specfunc/erfc.c \
   src/gsl/specfunc/zeta.c \
   src/gsl/specfunc/elementary.c \
   src/gsl/specfunc/gamma.c \
   src/gsl/specfunc/gamma_inc.c \
   src/gsl/rng/rng.c \
   src/gsl/rng/default.c \
   src/gsl/rng/mt.c \
   src/gsl/rng/types.c \
   src/gsl/randist/binomial.c \
   src/gsl/randist/binomial_tpe.c \
   src/gsl/randist/beta.c \
   src/gsl/randist/exponential.c \
   src/gsl/randist/geometric.c \
   src/gsl/randist/nbinomial.c \
   src/gsl/randist/poisson.c \
   src/gsl/randist/multinomial.c \
   src/gsl/randist/chisq.c \
   src/gsl/randist/gauss.c \
   src/gsl/randist/tdist.c \
   src/gsl/randist/gausszig.c \
   src/gsl/randist/gamma.c \
   src/gsl/randist/hyperg.c \
   src/gsl/cdf/binomial.c \
   src/gsl/cdf/beta.c \
   src/gsl/cdf/betainv.c \
   src/gsl/cdf/gauss.c \
   src/gsl/cdf/gaussinv.c \
   src/gsl/cdf/tdist.c \
   src/gsl/cdf/tdistinv.c \
   src/gsl/cdf/chisq.c \
   src/gsl/cdf/chisqinv.c \
   src/gsl/cdf/gamma.c \
   src/gsl/cdf/gammainv.c \
   src/gsl/cdf/hypergeometric.c \
   src/gsl/cdf/poisson.c \
   src/gsl/blas/blas.c \
   src/gsl/cblas/caxpy.c \
   src/gsl/cblas/ccopy.c \
   src/gsl/cblas/cdotc_sub.c \
   src/gsl/cblas/cdotu_sub.c \
   src/gsl/cblas/cgbmv.c \
   src/gsl/cblas/cgemm.c \
   src/gsl/cblas/cgemv.c \
   src/gsl/cblas/cgerc.c \
   src/gsl/cblas/cgeru.c \
   src/gsl/cblas/chbmv.c \
   src/gsl/cblas/chemm.c \
   src/gsl/cblas/chemv.c \
   src/gsl/cblas/cher2.c \
   src/gsl/cblas/cher2k.c \
   src/gsl/cblas/cher.c \
   src/gsl/cblas/cherk.c \
   src/gsl/cblas/chpmv.c \
   src/gsl/cblas/chpr2.c \
   src/gsl/cblas/chpr.c \
   src/gsl/cblas/cscal.c \
   src/gsl/cblas/csscal.c \
   src/gsl/cblas/cswap.c \
   src/gsl/cblas/csymm.c \
   src/gsl/cblas/csyr2k.c \
   src/gsl/cblas/csyrk.c \
   src/gsl/cblas/ctbmv.c \
   src/gsl/cblas/ctbsv.c \
   src/gsl/cblas/ctpmv.c \
   src/gsl/cblas/ctpsv.c \
   src/gsl/cblas/ctrmm.c \
   src/gsl/cblas/ctrmv.c \
   src/gsl/cblas/ctrsm.c \
   src/gsl/cblas/ctrsv.c \
   src/gsl/cblas/dasum.c \
   src/gsl/cblas/daxpy.c \
   src/gsl/cblas/dcopy.c \
   src/gsl/cblas/ddot.c \
   src/gsl/cblas/dgbmv.c \
   src/gsl/cblas/dgemm.c \
   src/gsl/cblas/dgemv.c \
   src/gsl/cblas/dger.c \
   src/gsl/cblas/dnrm2.c \
   src/gsl/cblas/drot.c \
   src/gsl/cblas/drotg.c \
   src/gsl/cblas/drotm.c \
   src/gsl/cblas/drotmg.c \
   src/gsl/cblas/dsbmv.c \
   src/gsl/cblas/dscal.c \
   src/gsl/cblas/dsdot.c \
   src/gsl/cblas/dspmv.c \
   src/gsl/cblas/dspr2.c \
   src/gsl/cblas/dspr.c \
   src/gsl/cblas/dswap.c \
   src/gsl/cblas/dsymm.c \
   src/gsl/cblas/dsymv.c \
   src/gsl/cblas/dsyr2.c \
   src/gsl/cblas/dsyr2k.c \
   src/gsl/cblas/dsyr.c \
   src/gsl/cblas/dsyrk.c \
   src/gsl/cblas/dtbmv.c \
   src/gsl/cblas/dtbsv.c \
   src/gsl/cblas/dtpmv.c \
   src/gsl/cblas/dtpsv.c \
   src/gsl/cblas/dtrmm.c \
   src/gsl/cblas/dtrmv.c \
   src/gsl/cblas/dtrsm.c \
   src/gsl/cblas/dtrsv.c \
   src/gsl/cblas/dzasum.c \
   src/gsl/cblas/dznrm2.c \
   src/gsl/cblas/hypot.c \
   src/gsl/cblas/icamax.c \
   src/gsl/cblas/idamax.c \
   src/gsl/cblas/isamax.c \
   src/gsl/cblas/izamax.c \
   src/gsl/cblas/sasum.c \
   src/gsl/cblas/saxpy.c \
   src/gsl/cblas/scasum.c \
   src/gsl/cblas/scnrm2.c \
   src/gsl/cblas/scopy.c \
   src/gsl/cblas/sdot.c \
   src/gsl/cblas/sdsdot.c \
   src/gsl/cblas/sgbmv.c \
   src/gsl/cblas/sgemm.c \
   src/gsl/cblas/sgemv.c \
   src/gsl/cblas/sger.c \
   src/gsl/cblas/snrm2.c \
   src/gsl/cblas/srot.c \
   src/gsl/cblas/srotg.c \
   src/gsl/cblas/srotm.c \
   src/gsl/cblas/srotmg.c \
   src/gsl/cblas/ssbmv.c \
   src/gsl/cblas/sscal.c \
   src/gsl/cblas/sspmv.c \
   src/gsl/cblas/sspr2.c \
   src/gsl/cblas/sspr.c \
   src/gsl/cblas/sswap.c \
   src/gsl/cblas/ssymm.c \
   src/gsl/cblas/ssymv.c \
   src/gsl/cblas/ssyr2.c \
   src/gsl/cblas/ssyr2k.c \
   src/gsl/cblas/ssyr.c \
   src/gsl/cblas/ssyrk.c \
   src/gsl/cblas/stbmv.c \
   src/gsl/cblas/stbsv.c \
   src/gsl/cblas/stpmv.c \
   src/gsl/cblas/stpsv.c \
   src/gsl/cblas/strmm.c \
   src/gsl/cblas/strmv.c \
   src/gsl/cblas/strsm.c \
   src/gsl/cblas/strsv.c \
   src/gsl/cblas/xerbla.c \
   src/gsl/cblas/zaxpy.c \
   src/gsl/cblas/zcopy.c \
   src/gsl/cblas/zdotc_sub.c \
   src/gsl/cblas/zdotu_sub.c \
   src/gsl/cblas/zdscal.c \
   src/gsl/cblas/zgbmv.c \
   src/gsl/cblas/zgemm.c \
   src/gsl/cblas/zgemv.c \
   src/gsl/cblas/zgerc.c \
   src/gsl/cblas/zgeru.c \
   src/gsl/cblas/zhbmv.c \
   src/gsl/cblas/zhemm.c \
   src/gsl/cblas/zhemv.c \
   src/gsl/cblas/zher2.c \
   src/gsl/cblas/zher2k.c \
   src/gsl/cblas/zher.c \
   src/gsl/cblas/zherk.c \
   src/gsl/cblas/zhpmv.c \
   src/gsl/cblas/zhpr2.c \
   src/gsl/cblas/zhpr.c \
   src/gsl/cblas/zscal.c \
   src/gsl/cblas/zswap.c \
   src/gsl/cblas/zsymm.c \
   src/gsl/cblas/zsyr2k.c \
   src/gsl/cblas/zsyrk.c \
   src/gsl/cblas/ztbmv.c \
   src/gsl/cblas/ztbsv.c \
   src/gsl/cblas/ztpmv.c \
   src/gsl/cblas/ztpsv.c \
   src/gsl/cblas/ztrmm.c \
   src/gsl/cblas/ztrmv.c \
   src/gsl/cblas/ztrsm.c \
   src/gsl/cblas/ztrsv.c \
   src/gsl/linalg/svd.c \
   src/gsl/linalg/bidiag.c \
   src/gsl/linalg/householder.c \
   src/gsl/matrix/matrix.c \
   src/gsl/matrix/submatrix.c \
   src/gsl/matrix/rowcol.c \
   src/gsl/matrix/getset.c \
   src/gsl/matrix/init.c \
   src/gsl/vector/init.c \
   src/gsl/matrix/swap.c \
   src/gsl/vector/vector.c \
   src/gsl/vector/copy.c \
   src/gsl/vector/swap.c \
   src/gsl/vector/subvector.c \
   src/gsl/vector/oper.c \
   src/gsl/matrix/oper.c \
   src/gsl/matrix/copy.c \
   src/gsl/block/init.c

GOBJS:= $(GSL:.c=.o)

help :
	@echo " "
	@echo "rareTDT version $(VERSION) source code"
	@echo " "
	@echo "Type ...          To ..."
	@echo "make rvTDT        Compile the rvTDT program"
	@echo "make test         run a quick test"
	@echo "make clean        Delete temporary files and ALL .o objects"
	@echo "make clear        Deep clean"
	@echo " "

EXE_B:= rvTDT
SRC_B:= ./src/tped.cpp ./src/tdtTrio.cpp ./src/Argument_helper.cpp ./src/rareTdt.cpp ./src/mainReal.cpp ./src/fisher2.cpp 
OBJS_B:= $(SRC_B:.cpp=.o)
rvTDT: $(EXE_B)
$(EXE_B): $(OBJS_B) $(JOBJS) $(GOBJS)
	$(CXX) $(INCLUDEDIR) $(OBJS_B) $(JOBJS) $(GOBJS) -o $(EXE_B) $(LIBDIR) $(LIBS) 

.cpp.o: $*.cpp
	$(CXX)  $(CXXFLAGS) $(INCLUDEDIR) -c -o $@  $<
.c.o: $*.c
	$(CXX)  $(CXXFLAGS) $(INCLUDEDIR) -c -o $@  $<

test: rvTDT
	./rvTDT rvTDT-test -G ./example/ABCA7.tped -P ./example/autism_phased.phen  -M ./example/ABCA7.map -p 2000	

.SUFFIXES : .cpp .c .o $(SUFFIXES)

clean:
	rm -rf $(SRCDIR)*.o $(SRCDIR)*~
clear:
	find $(SRCDIR) -type f -name "*.o" -delete
