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

module atomlib_trrio
  use iso_fortran_env, only: REAL64
  use xdrfor
  implicit none
  private
  public :: trnheader, trr_open, trr_close, trr_do_header, trr_do_real, trr_do_double, trr_do_data, trr_do_skip

  ! NOTE: xrdfile read is same as write. He he he.

  ! This is a fortran implementation of original GROMACS's code.
  ! * https://github.com/gromacs/gromacs/blob/master/src/gromacs/fileio/trrio.cpp

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: MAGIC_NUMBER = 1993 ! Magic constant
  character(*), parameter :: MAGIC_STRING = "GMX_trn_file" ! Magic string
  integer, parameter :: ERR_LEN = 128 ! Length of error string

  ! TRR header information 
  type trnheader
    integer :: ir_size = 0, e_size = 0, vir_size = 0, pres_size = 0, top_size = 0, sym_size = 0, nre = 0 ! Backwards compatibility
    integer :: box_size = 0 ! Non-zero if simulation box is present
    integer :: x_size = 0 ! Non-zero if particle coordinates are present
    integer :: v_size = 0 ! Non-zero if particle velocities are present
    integer :: f_size = 0 ! Non-zero if particle forces are present
    logical :: is_double = .false. ! Is trajectory double precision
    integer :: natoms = 0 ! Number of atoms
    integer :: step = 0 ! Number of frame
    real :: time = 0.0 ! Current time 
    real :: lambda = 0.0 ! Lambda value
  end type trnheader

contains

! Util: Determine if data is type real or double. 
integer function get_float_size (sh, stat, errmsg)
  use iso_fortran_env, only: REAL32, REAL64
  implicit none
  type(trnheader), intent(inout) :: sh
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  integer :: status

  get_float_size = 0

  catch: block

    status = exdrOK
  
    ! Calculate size depending on what we have
    if (sh%box_size /= 0) then
      get_float_size = sh%box_size / (DIM * DIM)
    else if (sh%x_size /= 0) then
      get_float_size = sh%x_size / (sh%natoms * DIM)
    else if (sh%v_size /= 0) then
      get_float_size = sh%v_size / (sh%natoms * DIM)
    else if (sh%f_size /= 0) then
      get_float_size = sh%f_size / (sh%natoms * DIM)
    else
      message = "Can not determine precision of TRR file"
      status = 1
      exit catch
    end if

    ! Is calculated size correct?
    if (get_float_size /= storage_size(0.0_REAL32) / 8 .and. get_float_size /= storage_size(0.0_REAL64) / 8) then
      message = "Unknown float size. Maybe different CPU?"
      status = 1
      exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end function get_float_size

! Open and check TRR file. Construct frame offset table. 
subroutine trr_open (file, xd, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: INT64
  implicit none
  character(*), intent(in) :: file ! Name of file.
  type(xdrfile), pointer, intent(inout) :: xd ! XDR file pointer.
  integer, intent(out) :: nframes ! Num. of frames in file.
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  integer, parameter :: SEEK_SET = 0
  character(ERR_LEN) :: message
  type(trnheader) :: sh
  integer(INT64) :: pos
  integer(INT64), allocatable :: temp(:)
  integer :: status

  catch: block

    ! Open XDR file in read mode.
    call xdr_open (xd, file, bRead=.true., STAT=status, ERRMSG=message)
    if (status /= exdrOK) exit catch

    ! Start at offset table of size 128. If needed, we will expand it as we go along.
    if (allocated(offset)) deallocate (offset, STAT=status, ERRMSG=message)
    allocate (offset(128), SOURCE=0_INT64, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    nframes = 0
    
    while: do while (.true.)
      ! Store current position
      pos = xdr_tell(xd)

      ! Read header and if EOF then there are (probably) no more frames in trajectory.
      call trr_do_header (xd, .true., sh, STAT=status, ERRMSG=message)
      if (status == exdrENDOFFILE) exit ! while
      if (status /= exdrOK) exit catch

      ! Expand offset table if needed
      if (nframes > size(offset)) then
        allocate (temp(2 * size(offset)), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        temp(1:nframes) = offset(1:nframes)
        call move_alloc (temp, offset)
      end if
      
      nframes = nframes + 1
      offset(nframes) = pos

      ! Skip remaining frame and only check for errors.
      call trr_do_skip (xd, sh, STAT=status, ERRMSG=message)
      if (status /= exdrOK) exit catch

    end do while

    ! Rewind trajectory file
    status = xdr_seek(xd, 0_INT64, SEEK_SET)
    if (status /= exdrOK) exit catch

    ! Rellocate offset table to actual size.
    allocate (temp(nframes), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    temp = offset(1:nframes)
    call move_alloc (temp, offset)

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_open

! Close TRR file.
subroutine trr_close (xd, stat, errmsg)
  implicit none
  type(xdrfile), pointer, intent(inout) :: xd ! XDR file pointer.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERR_LEN) :: message
  integer :: status

  call xdr_close (xd, STAT=status, ERRMSG=message)
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_close

! Read/write TRR file header.
subroutine trr_do_header (xd, bRead, sh, stat, errmsg)
  use iso_fortran_env, only: REAL64
  implicit none
  type(xdrfile), pointer, intent(in) :: xd
  logical, intent(in) :: bRead
  type(trnheader), intent(inout) :: sh
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  character(80) :: buffer
  real(REAL64) :: time, lambda
  logical :: bOK
  integer :: status, magic, slen

  catch: block

    status = exdrOK

    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    magic = MAGIC_NUMBER
    if (xdrfile_read_int(magic, 1, xd) /= 1) then
      message = "Can not read/write XDR file header"
      status = exdrENDOFFILE
    else if (magic /= MAGIC_NUMBER) then
      message = "Input file is not in TRR format"
      status = exdrMAGIC
    end if
    if (status /= 0) exit catch

    slen = len(MAGIC_STRING) + 1
    if (xdrfile_read_int(slen, 1, xd) /= 1) then
      message = "Can not read/write TRR magic string"
      status = exdrMAGIC
    else if (slen /= len(MAGIC_STRING) + 1) then
      message = "Can not read/write TRR magic string"
      status = exdrMAGIC
    end if
    if (status /= 0) exit catch

    if (bRead) then
      if (xdrfile_read_string(buffer, len(buffer), xd) <= 0) then
        message = "Can not read TRR version string"
        status = exdrMAGIC
      end if
    else
      if (xdrfile_write_string(MAGIC_STRING, xd) /= len(MAGIC_STRING) + 1) then
        message = "Can not write TRR version string"
        status = exdrMAGIC
      end if
    end if
    if (status /= 0) exit catch

    bOK = .true.
    bOK = bOK .and. (xdrfile_read_int(sh%ir_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%e_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%box_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%vir_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%pres_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%top_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%sym_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%x_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%v_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%f_size, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%natoms, 1, xd) == 1)
    if (.not. bOK) then
      message = "Can not read/write TRR header information"
      status = exdrHEADER
      exit catch
    end if

    ! Check trajectory precision
    sh%is_double = (get_float_size(sh, status, message) == storage_size(0.0_REAL64) / 8)
    if (status /= 0) exit catch

    bOK = bOK .and. (xdrfile_read_int(sh%step, 1, xd ) == 1)
    bOK = bOK .and. (xdrfile_read_int(sh%nre, 1, xd) == 1)
    if (sh%is_double) then
      time = real(sh%time, REAL64)
      lambda = real(sh%lambda, REAL64)
      bOK = bOK .and. (xdrfile_read_double(time, 1, xd) == 1)
      bOK = bOK .and. (xdrfile_read_double(lambda, 1, xd) == 1)
      sh%time = real(time)
      sh%lambda = real(lambda)
    else
      bOK = bOK .and. (xdrfile_read_float(sh%time, 1, xd) == 1)
      bOK = bOK .and. (xdrfile_read_float(sh%lambda, 1, xd) == 1)
    end if
    if (.not. bOK) then
      message = "Can not read/write TRR header information"
      status = exdrINT
      exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_do_header

! Read/write data in TRR float format.
subroutine trr_do_real (xd, sh, x, v, f, box, stat, errmsg)
  use iso_fortran_env, only: REAL32
  implicit none
  type(xdrfile), pointer, intent(in) :: xd
  type(trnheader), intent(in) :: sh
  real(REAL32), intent(inout) :: box(DIM,DIM)
  real(REAL32), intent(inout) :: x(DIM,*)
  real(REAL32), intent(inout) :: v(DIM,*)
  real(REAL32), intent(inout) :: f(DIM,*)
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  real(REAL32) :: dummy(DIM,DIM) = 0.0
  integer :: status
  logical :: bOK
  
  catch: block
    
    status = exdrOK

    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    if (sh%box_size /= 0) then
      if (xdrfile_read_float(box, DIM * DIM, xd) /= DIM * DIM) then
        message = "Can not read/write XDR simulation box data"
        status = exdrFLOAT
        exit catch
      end if
    end if
    
    bOK = .true.
    if (sh%vir_size /= 0) bOK = bOK .and. (xdrfile_read_float(dummy, DIM * DIM, xd) /= DIM * DIM)
    if (sh%pres_size /= 0) bOK = bOK .and. (xdrfile_read_float(dummy, DIM * DIM, xd) /= DIM * DIM)
    if (.not. bOK) then
      message = "Can not read/write legacy data"
      status = exdrFLOAT
      exit catch
    end if

    ! Read/write coordinates 
    if (sh%x_size /= 0) then
      if (xdrfile_read_float(x, sh%natoms * DIM, xd) /= sh%natoms * DIM) then
        message = "Can not read/write coordinate data"
        status = exdr3DX
        exit catch
      end if
    end if

    ! Read/write velocities
    if (sh%v_size /= 0) then
      if (xdrfile_read_float(v, sh%natoms * DIM, xd) /= sh%natoms * DIM) then
        message = "Can not read/write velocity data"
        status = exdr3DX
        exit catch
      end if
    end if
      
    ! Read/write forces
    if (sh%f_size /= 0) then
      if (xdrfile_read_float(f, sh%natoms * DIM, xd) /= sh%natoms * DIM) then
        message = "Can not read/write force data"
        status = exdr3DX
        exit catch
      end if
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_do_real

! Read/write data in TRR double format.
subroutine trr_do_double (xd, sh, x, v, f, box, stat, errmsg)
  use iso_fortran_env, only: REAL64
  implicit none
  type(xdrfile), pointer, intent(in) :: xd
  type(trnheader), intent(in) :: sh
  real(REAL64), intent(inout) :: x(DIM,*)
  real(REAL64), intent(inout) :: v(DIM,*)
  real(REAL64), intent(inout) :: f(DIM,*)
  real(REAL64), intent(inout) :: box(DIM,DIM)
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  real(REAL64) :: dummy(DIM,DIM) = 0d0
  character(128) :: message
  logical :: bOK
  integer :: status
  
  catch: block

    status = exdrOK

    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    if (sh%box_size /= 0) then
      if (xdrfile_read_double(box, DIM * DIM, xd) /= DIM * DIM) then
        message = "Can not read/write XDR simulation box data"
        status = exdrDOUBLE
        exit catch
      end if
    end if
  
    bOK = .true.
    if (sh%vir_size /= 0) bOK = bOK .and. (xdrfile_read_double(dummy, DIM * DIM, xd) /= DIM * DIM)
    if (sh%pres_size /= 0) bOK = bOK .and. (xdrfile_read_double(dummy, DIM * DIM, xd) /= DIM * DIM)
    if (.not. bOK) then
      message = "Can not read/write legacy data"
      status = exdrDOUBLE
      exit catch
    end if

    if (sh%x_size /= 0) then
      if (xdrfile_read_double(x, sh%natoms * DIM, xd) /= sh%natoms * DIM) then
        message = "Can not read/write coordinate data"
        status = exdr3DX
        exit catch
      end if
    end if
    
    if (sh%v_size /= 0) then
      if (xdrfile_read_double(v, sh%natoms * DIM, xd) /= sh%natoms * DIM) then
        message = "Can not read/write velocities data"
        status = exdr3DX
        exit catch
      end if
    end if

    if (sh%f_size /= 0) then
      if (xdrfile_read_double(f, sh%natoms * DIM, xd) /= sh%natoms * DIM) then
        message = "Can not read/write forces data"
        status = exdr3DX
        exit catch
      end if
    end if 

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_do_double

! Wrapper for reading/writing data in TRR format.
! NOTE: Loss of data if used on double precision trajectory. Use 'trr_do_double'.
subroutine trr_do_data (xd, sh, coor, box, stat, errmsg)
  use iso_fortran_env, only: REAL32, REAL64
  implicit none
  type(xdrfile), pointer, intent(in) :: xd
  type(trnheader), intent(in) :: sh
  real, intent(inout) :: coor(DIM,sh%natoms)
  real, intent(inout) :: box(DIM,DIM)
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  real(REAL64), allocatable :: x_double(:,:), dmy_double(:,:)
  real(REAL64) :: box_double(DIM,DIM)
  real, allocatable :: dmy_float(:,:)
  character(ERR_LEN) :: message
  integer :: status, natoms

  catch: block

    status = exdrOK

    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    natoms = sh%natoms

    if (sh%is_double) then
      allocate (x_double(DIM,natoms), SOURCE=real(coor,REAL64), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      allocate (dmy_double(DIM,natoms), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

      call trr_do_double (xd, sh, x_double, dmy_double, dmy_double, box_double, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

      coor = real(x_double, REAL32)
      box = real(box_double, REAL32)

    else ! is float
      allocate (dmy_float(DIM,natoms), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      call trr_do_real (xd, sh, coor, dmy_float, dmy_float, box, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    end if 

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_do_data 

! Read through TRR file but do not store any data.
subroutine trr_do_skip (xd, sh, stat, errmsg)
  use iso_fortran_env, only: REAL64, REAL32, INT64
  implicit none
  type(xdrfile), pointer, intent(in) :: xd
  type(trnheader), intent(in) :: sh
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  integer(INT64) :: framebytes
  integer :: status

  catch: block
  
    status = exdrOK

    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    framebytes = 0
    if (sh%box_size /= 0) framebytes = framebytes + DIM * DIM
    if (sh%vir_size /= 0) framebytes = framebytes + DIM * DIM
    if (sh%pres_size /= 0) framebytes = framebytes + DIM * DIM
    if (sh%x_size /= 0) framebytes = framebytes + sh%natoms * DIM
    if (sh%v_size /= 0) framebytes = framebytes + sh%natoms * DIM
    if (sh%f_size /= 0) framebytes = framebytes + sh%natoms * DIM
    if (sh%is_double) then
      framebytes = framebytes * storage_size(0.0_REAL64) / 8
    else
      framebytes = framebytes * storage_size(0.0_REAL32) / 8
    end if

    ! Skip number of bytes from current position (SEEK_CUR = 1).
    if (xdr_seek(xd, framebytes, 1) /= exdrOK) then
      message = "Can not skip TRR frame data"
      status = exdr3DX
      exit catch
    end if
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine trr_do_skip

end module atomlib_trrio
