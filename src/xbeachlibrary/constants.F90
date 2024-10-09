!****h* General/Constants
!
!  SYNOPSIS
!     Module defining common constants
!
!  COPYRIGHT
!     (C) 2004 by WL|Delft Hydraulics
!
!  AUTHOR
!     Arjen Markus
!
!  DESCRIPTION
!  This module contains integer and real parameters:
!  - kind parameters for single (sp) and double (dp) precision
!  - definitions for pi (3.14...) and compi (imaginary unit)
!
! SOURCE
!
! wwvv never used
!      this would be a nice place to hold things as pi, twopi
!      and so on

module constants
   implicit none
   save
   private
   public spkind,dpkind,pi,compi,iFill,sFill,dFill

   integer,              parameter :: spkind = kind(1.0)
   integer,              parameter :: dpkind = kind(1.0d0)
   real(kind=dpkind),    parameter :: pi     = 4*atan(1.0_dpkind)
   complex(kind=dpkind), parameter :: compi  = (0.0_dpkind,1.0_dpkind)

   ! Fill numbers
   integer,parameter                     :: iFill = -huge(0)
   real,parameter                        :: sFill = -huge(0.0)
   real*8,parameter                      :: dFill = -dble(huge(0.0))  ! Robert: easier this way than to catch dFill<sFill later on

end module constants
!****
