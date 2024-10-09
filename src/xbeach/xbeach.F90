program xbeach
   use libxbeach_module,     only: executestep, init, final, outputext
   use introspection_module, only: getdoubleparameter
   use iso_c_binding,        only: c_double
   use process_input,        only: readinput
   use xmpi_module,          only: halt_program
   implicit none

   integer                                             :: rc
   real(c_double)                                      :: t, tstop

   ! Initialize program                                                          !
   rc = 0
   rc = readinput()
   if (rc.eq.1) then
      call halt_program
   endif

   rc = init()

   ! Start simulation                                                            !
   rc = getdoubleparameter("t", t)
   rc = getdoubleparameter("tstop", tstop)
   do while (t<tstop)
      rc = executestep()
      rc = getdoubleparameter("t", t)
      ! output
      rc = outputext()
   enddo

   ! Cleanup
   rc = final()
end program
