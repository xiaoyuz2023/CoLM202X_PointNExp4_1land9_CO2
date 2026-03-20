# =======================================================
# mpif90 - ifort 
# 

 FF = mpif90

 NETCDF_LIB = /stu01/lianghb21/CoLM202X/2023student/lib/netcdf/netcdf/lib
 NETCDF_INC = /stu01/lianghb21/CoLM202X/2023student/lib/netcdf/netcdf/include

 LAPACK_LIB = /stu01/lianghb21/CoLM202X/2023student/lib/lapack/lapack-3.11/lapack-311/lib

 MOD_CMD = -J

 GCC_VERSION := "`gcc -dumpversion`"
  IS_GCC_ABOVE_10 := $(shell expr "$(GCC_VERSION)" ">=" "10")
  ifeq "$(IS_GCC_ABOVE_10)" "1" 
     FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
           -ffpe-trap=invalid,zero,overflow -fbounds-check \
           -mcmodel=medium -fbacktrace -fdump-core -cpp \
           -ffree-line-length-0 -fallow-argument-mismatch 
  else
     FOPTS = -fdefault-real-8 -ffree-form -C -g -u -xcheck=stkovf \
           -ffpe-trap=invalid,zero,overflow -fbounds-check \
           -mcmodel=medium -fbacktrace -fdump-core -cpp \
           -ffree-line-length-0
  endif

 LDFLAGS = -L${NETCDF_LIB} -lnetcdff -L${LAPACK_LIB} -llapack -lblas