#!/bin/bash


## start with:  svn checkout url foldername
## cd foldername

# VERSION SPECIFIC PARAMETER!
export XBEACH_PROJECT_ID=xbx

export MODULEPATH=$MODULEPATH:/opt/apps/modules


# start with anaconda since it overwrites gcc !!!
module load anaconda3/2021.05

module load gcc/7.3.0
module load hdf5/1.12.0_gcc7.3.0
module load netcdf/v4.7.4_v4.5.3_gcc7.3.0
module load openmpi/4.0.4_gcc7.3.0

make distclean

./autogen.sh

mkdir -p "/opt/apps/xbeach/"$XBEACH_PROJECT_ID"_gcc_7.3.0_openmpi_4.0.4_HEAD"

FCFLAGS="-mtune=corei7-avx -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -O3 -ffast-math" ./configure  --with-netcdf --with-mpi --prefix="/opt/apps/xbeach/"$XBEACH_PROJECT_ID"_gcc_7.3.0_openmpi_4.0.4_HEAD"

make
make install

#cd $SVNPATH"/install"

#tar -cf $SVNPATH"/artifacts.tar" bin lib

cd /config/

/usr/share/Modules/bin/createmodule.py -p "/opt/apps/xbeach/"$XBEACH_PROJECT_ID"_gcc_7.3.0_openmpi_4.0.4_HEAD" ./teamcity-env.sh > "/opt/apps/modules/xbeach/xbeach-"$XBEACH_PROJECT_ID"_gcc_7.3.0_openmpi_4.0.4_HEAD"

chmod 777 "/opt/apps/modules/xbeach/xbeach-"$XBEACH_PROJECT_ID"_gcc_7.3.0_openmpi_4.0.4_HEAD"