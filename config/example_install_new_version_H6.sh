#!/bin/bash
. /etc/profile

# VERSION SPECIFIC PARAMETER!
export XBEACH_PROJECT_ID=xbx
export SVNNAME=xb_xbx


export MODULEPATH=$MODULEPATH:/opt/modules

export SVNPATH=$(pwd)
cd $SVNNAME

# start with anaconda since it overwrites gcc !!!
module load anaconda/lnx64_conda

module load gcc/4.9.2
module load hdf5/1.8.14_gcc_4.9.2
module load netcdf/v4.3.2_v4.4.0_gcc_4.9.2
module load openmpi/1.8.3_gcc_4.9.2

make distclean

./autogen.sh

mkdir -p "/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

FCFLAGS="-mtune=corei7-avx -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -O3 -ffast-math" ./configure  --with-netcdf --with-mpi --prefix="/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

make
make install

cd $SVNPATH"/install"

tar -cf $SVNPATH"/artifacts.tar" bin lib

cd $SVNPATH"/"$SVNNAME"/config/"

/usr/share/Modules/bin/createmodule.py -p "/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD" ./teamcity-env.sh > "/opt/xbeach/modules/xbeach-"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

chmod 777 "/opt/xbeach/modules/xbeach-"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"