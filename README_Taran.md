Xbeach
svn clone xbeach repository path

./autogen.sh
./configure --prefix=/usr/include
make

Path to run
/home/tarandeepk/Desktop/Xbeach_2/trunk/src/xbeach

./configure --prefix=/usr/include CPPFLAGS="-I/usr/include" LDFLAFS="-L/usr/include" NETCDF_INCIDIR="/usr/include" NETCDF_LIBDIR="/usr/include" --with-netcdf --with-mpi
mpirun -np 2 ../xbeach
