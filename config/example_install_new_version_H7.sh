#!/bin/bash


## start with:  svn checkout url foldername
## cd foldername

# VERSION SPECIFIC PARAMETER!
export XBEACH_PROJECT_ID=xbx

export MODULEPATH=$MODULEPATH:/opt/apps/modules


# start with anaconda since it overwrites gcc !!!
module load anaconda3/py39_23.1.0-1

module load gcc/12.2.0_gcc12.2.0
module load hdf5/1.14.0_gcc12.2.0
module load netcdf/4.9.2_4.6.1_gcc12.2.0
module load openmpi/4.1.5_gcc12.2.0
module load curl

#module load intel/2023.1.0
#module load hdf5/1.12.0_intel2023.1.0
#module load netcdf/4.9.2_4.6.1_intel2023.1.0
#module load openmpi/5.0.0_intel2023.1.0


make distclean

./autogen.sh

mkdir -p "/opt/apps/xbeach/"$XBEACH_PROJECT_ID"_gcc_12.2.0_openmpi_4.1.5_HEAD"

FCFLAGS="-I/opt/apps/netcdf/4.9.2_4.6.1_gcc12.2.0/include -mtune=corei7-avx -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -O3 -ffast-math" ./configure  --with-netcdf --with-mpi --prefix="/opt/apps/xbeach/"$XBEACH_PROJECT_ID"_gcc_12.2.0_openmpi_4.1.5_HEAD"

make
make install

#cd $SVNPATH"/install"

#tar -cf $SVNPATH"/artifacts.tar" bin lib

cd /config/

/usr/share/Modules/bin/createmodule.py -p "/opt/apps/xbeach/"$XBEACH_PROJECT_ID"_gcc_12.2.0_openmpi_4.1.5_HEAD" ./teamcity-env.sh > "/opt/apps/modules/xbeach/xbeach-"$XBEACH_PROJECT_ID"_gcc_12.2.0_openmpi_4.1.5_HEAD"

chmod 777 "/opt/apps/modules/xbeach/xbeach-"$XBEACH_PROJECT_ID"_gcc_12.2.0_openmpi_4.1.5_HEAD"