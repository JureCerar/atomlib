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

module atomlib_xtcio
  use xdrfor
  implicit none
  private
  public :: xtc_open, xtc_close, xtc_do_header, xtc_do_data, xtc_do_skip

  ! NOTE: xrdfile read is same as write. He he he.

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: MAGIC_NUMBER = 1995 ! Magic constant
  integer, parameter :: XTC_SHORTHEADER_SIZE = (20 + DIM * DIM * 4) ! Num of bytes in header: magic, natoms, step, time, and box
  integer, parameter :: XTC_SHORT_BYPESPERATOM = 12 ! Short XTCs store each coordinate as a 32-bit float.
  integer, parameter :: XTC_HEADER_SIZE = (DIM * DIM * 4 + DIM * 2 + 46) ! Short XTCs store each coordinate as a 32-bit float.
  integer, parameter :: ERR_LEN = 128 ! Length of error string

contains

! Check XTC file and construct frame offset table.
subroutine xtc_open (file, xd, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: INT64
  implicit none
  character(*), intent(in) :: file ! Name of file.
  type(xdrfile), pointer, intent(inout) :: xd ! XDR file pointer.
  integer, intent(out) :: nframes ! Num. of frames in file
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  integer, parameter :: SEEK_SET = 0
  character(ERR_LEN) :: message
  integer(INT64), allocatable :: temp(:)
  integer(INT64) :: pos
  integer :: status, natoms, step
  real :: time

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
      call xtc_do_header (xd, natoms, step, time, STAT=status, ERRMSG=message)
      if (status == exdrENDOFFILE) exit while
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
      call xtc_do_skip (xd, natoms, STAT=status, ERRMSG=message)
      if (status /= exdrOK) exit catch

    end do while

    ! Rewind trajectory file
    status = xdr_seek(xd, 0_INT64, SEEK_SET)
    if (status /= exdrOK) exit catch

    ! Relocate offset table to actual size.
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

end subroutine xtc_open

! Close XTC file.
subroutine xtc_close (xd, stat, errmsg)
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

end subroutine xtc_close

! Read/write XTC file header
subroutine xtc_do_header (xd, natoms, step, time, stat, errmsg)
  implicit none
  type(xdrfile), pointer, intent(in) :: xd ! XDR file pointer
  integer, intent(inout) :: natoms ! Num. of atoms in frame
  integer, intent(inout) :: step ! Current step number
  real,  intent(inout) :: time ! Simulation time in [ps]
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  integer :: status, magic
  logical :: bOK

  catch: block

    status = 0
    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    ! Read/write and check magic number
    magic = MAGIC_NUMBER
    if (xdrfile_read_int(magic, 1, xd) /= 1) then
      message = "Can not read/write XDR file header"
      status = exdrENDOFFILE
    else if (magic /= MAGIC_NUMBER) then
      message = "Input file is not in XTC format"
      status = exdrMAGIC
    end if
    if (status /= 0) exit catch

    ! Read/write number of atoms, step, and time
    bOK = .true.
    bOK = bOK .and. (xdrfile_read_int(natoms, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_int(step, 1, xd) == 1)
    bOK = bOK .and. (xdrfile_read_float(time, 1, xd) == 1)
    if (.not. bOK) then
      message = "Can not read/write XTC header information"
      stat = exdrHEADER
      exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine xtc_do_header

! Read/write XTC file data
subroutine xtc_do_data (xd, bRead, natoms, prec, coor, box, stat, errmsg)
  implicit none
  type(xdrfile), intent(in), pointer :: xd ! XDR file pointer
  logical, intent(in) :: bRead ! Read or write to file
  integer, intent(in) :: natoms ! Number of atoms in file
  real, intent(inout) :: prec ! Coordinate precission
  real, intent(inout) :: coor(DIM,natoms) ! Particle coordinates
  real, intent(inout) :: box(DIM,DIM) ! Simulation box size
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  integer :: status

  catch: block

    status = 0
    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if
  
    if (xdrfile_read_float(box, DIM * DIM, xd) /= DIM * DIM) then
      message = "Can not read/write XDR simulation box data"
      status = exdrFLOAT
      exit catch
    end if
  
    ! XTC coordinates are compressed: read = decompression and write = compression
    if (bRead) then
      if (xdrfile_decompress_coord_float(coor, natoms, prec, xd) /= natoms) then
        message = "Can not read decompressed XDR coordinates"
        status = exdr3DX
        exit catch
      end if
    else
      if (xdrfile_compress_coord_float(coor, natoms, prec, xd) /= natoms) then
        message = "Can not write compressed XDR coordinates"
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

end subroutine xtc_do_data

! Skip XTC file data
subroutine xtc_do_skip (xd, natoms, stat, errmsg)
  use iso_fortran_env, only: INT64
  implicit none
  type(xdrfile), intent(in), pointer :: xd ! XDR file pointer
  integer, intent(in) :: natoms ! Number of particles in frame
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  integer, parameter :: SEEK_CUR = 1
  character(ERR_LEN) :: message
  integer :: nbytes
  integer(INT64) :: framebytes
  integer :: status

  catch: block

    status = 0
    if (.not. associated(xd)) then
      message = "File XDR pointer not associated"
      status = exdrFILENOTFOUND
      exit catch
    end if

    ! If less than 10 atoms XTC file is not compressed
    if (natoms < 10) then
      ! Calculate num. of bytes in frame: box, dummy, and coor
      ! * I do not know what "dummy" variable is, only that it's 4 bytes long.
      framebytes = DIM * DIM * 4 + 4 + natoms * XTC_SHORT_BYPESPERATOM
      if (xdr_seek(xd, framebytes, 1) /= exdrOK) then
        message = "Can not skip XDR frame data"
        status = exdr3DX
        exit catch
      end if
  
    else
      ! Calculate num. of bytes in frame: box, and dummy
      ! * I do not know what "dummy" data is only that it's 36 bytes long.
      framebytes = DIM * DIM * 4 + 36
      if (xdr_seek(xd, framebytes, SEEK_CUR) /= exdrOK) then
        message = "Can not skip XDR frame data"
        status = exdr3DX
        exit catch
      end if
  
      ! Read actual byte size of compressed coordinates
      if (xdrfile_read_int(nbytes, 1, xd) /= 1) then
        message = "Can not skip XDR frame data"
        status = exdrINT
        exit catch
      end if
  
      ! Skip specified number of bytes but fitst round to next 32-bit boundry.
      ! C source: int64_t framebytes = (nbytes + 3) & ~0x03
      framebytes = int(and((nbytes + 3), not(int(z'03'))), INT64)
      if (xdr_seek(xd, framebytes, SEEK_CUR) /= exdrOK) then
        message = "Can not skip XDR frame data"
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

end subroutine xtc_do_skip

end module atomlib_xtcio
