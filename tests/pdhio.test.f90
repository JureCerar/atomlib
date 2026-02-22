! This file is part of atomlib
! https://github.com/JureCerar/atomlib
!
! Copyright (C) 2022 Jure Cerar
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

program main
  use atomlib_pdhio
  implicit none
  integer, parameter :: NP = 100
  real, parameter :: DELTA = 0.001
  type(pdh_t) :: pdh
  character(128) :: file, message
  integer :: i, status

  ! Get file name form command line argument
  call get_command_argument (1, file, STATUS=status)
  if (status /= 0) error stop

  ! Allocate reallocate and check size
  call pdh%allocate (NP)
  call pdh%reallocate (2*NP)
  call pdh%reallocate (NP/2)
  call pdh%deallocate ()

  ! Allocate reallocate and check size
  call pdh%allocate (NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call pdh%reallocate (2*NP, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call pdh%reallocate (NP/2, STAT=status, ERRMSG=message)
  if (status /= 0) error stop
  call pdh%deallocate (STAT=status, ERRMSG=message)
  if (status /= 0) error stop

  ! Construct a dummy configuration
  call pdh%allocate(NP)
  pdh%text = "This is a template."
  pdh%key_words = [("XXXX", i = 1, PDH_NUM_CHAR_CONST)]
  pdh%int_const = [(i, i = 1, PDH_NUM_INT_CONST)]
  pdh%real_const = [(real(i), i = 1, PDH_NUM_REAL_CONST)]
  pdh%x = [(real(i), i = 1, NP)]
  pdh%y = [(real(i), i = 1, NP)]
  pdh%y_error = [(real(i), i = 1, NP)]

  ! Construct a dummy configuration
  pdh = pdh_t(NP, TEXT="This is a template.", KEY_WORDS=[("XXXX", i=1,PDH_NUM_CHAR_CONST)], &
  & INT_CONST=[(i, i=1,PDH_NUM_INT_CONST)], REAL_CONST=[(real(i), i=1,PDH_NUM_REAL_CONST)], &
  & X=[(real(i), i=1,NP)], Y=[(real(i), i=1,NP)], Y_ERROR=[(real(i), i=1,NP)])

  ! Try derived type write functions
  write (*, *, IOSTAT=status, IOMSG=message) pdh
  if (status /= 0) error stop message
  write (*,"(DT)", IOSTAT=status, IOMSG=message) pdh
  if (status /= 0) error stop message
  write (*,"(DT(14,6))", IOSTAT=status, IOMSG=message) pdh
  if (status /= 0) error stop message

  ! Process file
  write (*,*) "Processing file: '", trim(file), "'"

  ! Load (read) file and check values
  call pdh%load (file)
  if (pdh%num_points /= 100) error stop
  if (pdh%key_words(1) /= "SAXS") error stop
  if (pdh%int_const(2) /= 0) error stop
  if (abs(pdh%real_const(1) - 0.000) > DELTA) error stop
  if (abs(pdh%x(1) - 1.000) > DELTA) error stop
  if (abs(pdh%y(1) - 1.000) > DELTA) error stop
  if (abs(pdh%y_error(1) + 1.000) > DELTA) error stop

  ! Load (read)  file w/ error handling
  call pdh%load (file, STAT=status, ERRMSG=message)
  if (status /= 0) error stop  

  ! Write (save) to file.
  call pdh%save ("out.pdh")
    
  ! Write file w/ error handling
  call pdh%save ("out.pdh", STAT=status, ERRMSG=message)
  if (status /= 0) error stop

end program main