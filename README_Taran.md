Xbeach
aa) svn clone xbeach repository path

1) ./autogen.sh
2) ./configure 
3) --prefix=/usr/include
4) make

Path to run
/home/tarandeepk/Desktop/Xbeach_2/trunk/src/xbeach

./configure --prefix=/usr/include CPPFLAGS="-I/usr/include" LDFLAFS="-L/usr/include" NETCDF_INCIDIR="/usr/include" NETCDF_LIBDIR="/usr/include" --with-netcdf --with-mpi


5) mpirun -np 2 ../xbeach
