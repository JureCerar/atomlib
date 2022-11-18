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

module atomlib
  use xdrfor
  use atomlib_xyzio
  use atomlib_groio
  use atomlib_pdbio
  use atomlib_crdio
  use atomlib_trrio
  use atomlib_xtcio
  use atomlib_dcdio
  use atomlib_confio
  use atomlib_trajio
  use atomlib_ndxio
  use atomlib_pdhio
  implicit none

  character(*), parameter :: atomlib_version = "@PROJECT_VERSION@"

end module atomlib