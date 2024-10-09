! netcdf_test.f90
program netcdf_test
  use netcdf
  implicit none

  integer :: ncid, status

  ! Create a new NetCDF file
  status = nf90_create('test.nc', NF90_CLOBBER, ncid)
  if (status /= NF90_NOERR) then
    print *, 'Error creating NetCDF file'
    stop
  end if

  ! Close the NetCDF file
  status = nf90_close(ncid)
  if (status /= NF90_NOERR) then
    print *, 'Error closing NetCDF file'
    stop
  end if

  print *, 'NetCDF installation seems OK!'
end program netcdf_test
