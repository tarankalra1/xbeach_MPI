Xbeach
A) get this repository 

1) ./autogen.sh
2) ./configure 
3) --prefix=/usr/include
4) make

Path to run
(Change this to your local path)
/home/tarandeepk/Desktop/Xbeach_2/trunk/src/xbeach

./configure --prefix=/usr/include CPPFLAGS="-I/usr/include" LDFLAFS="-L/usr/include" NETCDF_INCIDIR="/usr/include" NETCDF_LIBDIR="/usr/include" --with-netcdf --with-mpi


Take a set of xbeach runs 
B)  mpirun -np 2 ../xbeach
