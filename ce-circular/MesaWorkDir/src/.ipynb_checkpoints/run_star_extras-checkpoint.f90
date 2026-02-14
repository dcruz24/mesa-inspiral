! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      implicit none
      
      include 'test_suite_extras_def.inc'

        ! these routines are called by the standard run_star check_model
      contains

      !include 'standard_run_star_extras.inc'

      include "test_suite_extras.inc"

     subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_momentum => my_other_momentum_routine
         s% other_wind => my_other_wind_routine
      end subroutine extras_controls

     subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

 integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step

    subroutine my_other_momentum_routine(id, ierr)
         ! Adds potential to donor star from the compact star...
         real :: compactM, compactA, compactU
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         compactM = s% x_ctrl(1) ! [g]
         compactA = s% x_ctrl(2) ! [cm], Position of Accretor
         ! Notes in './star_data/public/star_data_step_work.inc'
         ! extra gravity (can be set by user)  added to -G*m/r^2 in momentum equation
         !type(auto_diff_real_star_order1), pointer, dimension(:) :: extra_grav
         ! type(auto_diff_real_star_orde1): Contains real(dp):: val & real(dp) :: 1dArray(33)
         ! Assume accessing val?
         
         do k = 1, s% nz
             s% extra_grav(k)%val = 1       ! Default hook had s% extra_grav(k) = 0d0 but get error with this
             if ( exp( s% lnR(k)) > compactA ) then
                 compactU = - standard_cgrav * compactM/(compactA**2)
                 s% extra_grav(k)%val = compactU       ! Default hook value has this, but cant match real auto_diff_real_star_order1
             end if
         end do
         
    end subroutine my_other_momentum_routine

    subroutine my_other_wind_routine(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), dimension(:), allocatable :: Ms, Rs, Vs, Rhos
         real :: vesc
         integer :: i, k, limitCell

         call star_ptr(id, s, ierr)         
         allocate(Ms(s% nz), Rs(s% nz), Vs(s% nz), Rhos(s% nz))
         ierr = 0
         limitCell = 0         
         Rs = exp(s% xh(s% i_lnR,:))
         Vs = s% xh(s% i_v, :)
         Rhos = exp(s% xh(s% i_lnD, :))
         
         !print *, 'Outer Radius: ', s% xh(s% i_lnR, 1)
         !print *, 'Inner Radius: ', s% xh(s% i_lnR, s% nz)
         
         ! Calculate the Enclosed Mass....
         Ms(1) = 0.0
         do i =2, s% nz
             Ms(i) = Ms(i-1) + 0.5 * (Rs(i) - Rs(i-1)) * (4.0 * pi * Rs(i)**2 * Rhos(i) + 4.0 * pi * Rs(i-1)**2 * Rhos(i-1))
         end do

         ! Find the last index where cell's velocity greater than escape velocity...
         do k = 1, s% nz
             vesc = sqrt(2*standard_cgrav*Ms(k)/Rs(k))
             if (Vs(k) > vesc) then
                 limitCell = k
             else
                 exit
             end if
         end do

         if (limitCell == 0) then
             w = 0
         else
             w = (Ms(1) - Ms(limitCell+1)) / Msun * 100  ! Remove this much mass..., should it be 100? Secds to years?
         end if
         deallocate(Rs, rhos, Ms, Vs)

      end subroutine my_other_wind_routine
      
    end module run_star_extras
