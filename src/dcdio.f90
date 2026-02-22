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

module atomlib_dcdio
  use, intrinsic :: iso_fortran_env, only: INT64, REAL64
  implicit none
  private
  public :: dcd_open, dcd_close, dcd_read_header, dcd_write_header  
  public :: dcd_read_natoms, dcd_read_data, dcd_skip_data, dcd_write_data

  ! TODO:
  ! - [ ] Handle remarks fixed length read.
  ! - [ ] Imporve error messages (see notes).
  ! - [ ] Implement proper box vector transform. 

  ! Global constants
  integer, parameter :: DIM = 3 ! Default dimension
  integer, parameter :: MAGIC_NUMBER = 84 ! Magic number
  character(*), parameter :: MAGIC_STRING = "CORD" ! Magic string
  integer, parameter :: ERR_LEN = 128 ! length of error string
  integer, parameter :: REMARK_LEN = 80 ! Length of remark string

contains

! Open and check DCD file (DCD files have some convoluted file check procedures). Construct frame offset table.
subroutine dcd_open (file, unit, nframes, offset, stat, errmsg)
  use iso_fortran_env, only: IOSTAT_END
  implicit none
  character(*), intent(in) :: file ! File name
  integer, intent(out) :: unit ! Fortran file unit
  integer, intent(out) :: nframes ! Num. of frames in file.
  integer(INT64), allocatable, intent(out) :: offset(:) ! Stream position offsets.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(REMARK_LEN), allocatable :: remark(:)
  character(ERR_LEN) :: message
  integer :: start_time, every_time, end_time
  integer :: magic_num, charmm_version, has_extra_block, four_dimensions
  real :: timestep
  character(4) :: magic_str
  integer :: i, status, natoms

  catch: block

    ! Open file in native endinness.
    open (NEWUNIT=unit, FILE=trim(file), FORM="unformatted", ACCESS="stream", STATUS="old", &
    &     IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Read magic number and magic string
    read (unit, POS=1, IOSTAT=status, IOMSG=message) magic_num, magic_str
    if (status /= 0) exit catch

    ! If MAGIC number and string are incorrect, close file and try swapping endinness
    if (magic_num /= MAGIC_NUMBER .or. magic_str /= MAGIC_STRING) then
      close (unit, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch
      open (NEWUNIT=unit, FILE=trim(file), FORM="unformatted", ACCESS="stream", &
      & STATUS="old", CONVERT="swap", IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

      ! Read magic number and magic string
      read (unit, POS=1, IOSTAT=status, IOMSG=message) magic_num, magic_str
      if (status /= 0) exit catch

      ! We tried both native and reverse endiness and failed
      if (magic_num /= MAGIC_NUMBER .or. magic_str /= MAGIC_STRING) then
        message = "Can not read file header: File not in DCD format"
        status = 1
        exit catch
      end if

    end if

    ! Check if the file identifies as CHARMM file (LAMMPS files pretend to be CHARMM v24).
    read (unit, POS=85, IOSTAT=status, IOMSG=message) charmm_version
    if (charmm_version == 0) then
      message = "File version not supported"
      status = 1
    end if
    if (status /= 0) exit catch 

    ! We only support files with the extra unitcell block.
    read (unit, POS=49, IOSTAT=status, IOMSG=message) has_extra_block
    if (has_extra_block /= 1) then
      message = "File version not supported"
      status = 1
    end if
    if (status /= 0) exit catch 

    ! We don't support files with >3 dimensions
    read (unit, IOSTAT=status, IOMSG=message) four_dimensions
    if (four_dimensions == 1) then
      message = "Files with >3 dimensions are not currently supported"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Now we can do actual file check. Return to start of file
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! DCD header already has number of frames, so less work for us.
    call dcd_read_header (unit, nframes, remark, start_time, every_time, end_time, &
    & timestep, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) return

    ! We cal allocate exact size of offset table. No need for reallocation calls. 
    if (allocated(offset)) deallocate (offset, STAT=status, ERRMSG=message)
    allocate (offset(nframes), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    ! Skip remaining frame and only check for errors. Store frame positions.
    do i = 1, nframes
      inquire (UNIT=unit, POS=offset(i), IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch
      call dcd_skip_data (unit, natoms, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end do

    ! Return to start of file
    rewind (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

  end block catch

  ! Error handling
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine dcd_open

! Close DCD file.
subroutine dcd_close (unit, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERR_LEN) :: message
  integer :: status

  close (unit, IOSTAT=status, IOMSG=message)
  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine dcd_close

! Read DCD header. DCD files have only one header.
subroutine dcd_read_header (unit, nframes, remark, start_time, every_time, end_time, timestep, natoms, stat, errmsg)
  use iso_fortran_env, only: INT64
  use iso_c_binding, only: C_NULL_CHAR
  implicit none
  integer, intent(in) :: unit
  character(*), allocatable, intent(inout) :: remark(:)
  integer, intent(out) :: nframes
  integer, intent(out) :: start_time
  integer, intent(out) :: every_time
  integer, intent(out) :: end_time
  real, intent(out) :: timestep
  integer, intent(out) :: natoms
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  integer(INT64) :: filesize, framesize, nframes2
  character(ERR_LEN) :: message
  character(REMARK_LEN) :: buffer
  logical :: opened
  character(9) :: action, access
  integer :: i, nremarks, n, dummy, pos, status

  catch: block
  
    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit must be opened for reading"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Read number of frames
    ! ERROR: "Can not read number of frames or the file header is corrupt"
    read (unit, POS=9, IOSTAT=status, IOMSG=message) nframes
    if (status /= 0) exit catch

    ! Read simulation times
    ! DEBUG: Check if POS is OK (4 bytes)
    ! ERROR: "Can not read time data or the file header is corrupt"
    read (unit, POS=13, IOSTAT=status, IOMSG=message) start_time, every_time, end_time
    if (status /= 0) exit catch

    ! Read timestep
    ! ERROR: "Can not read timestep data or the file header is corrupt"
    read (unit, POS=45, IOSTAT=status, IOMSG=message) timestep
    if (status /= 0) exit catch

    ! Read number of remarks and remark(s)
    read (unit, POS=97, IOSTAT=status, IOMSG=message) nremarks
    if (status /= 0) exit catch

    ! Allocate and read remarks. And trim C_NULL_CHAR
    if (allocated(remark)) deallocate(remark)
    allocate (remark(nremarks), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    do i = 1, nremarks
      read (unit, IOSTAT=status, IOMSG=message) buffer
      if (status /= 0) exit catch
      n = index(buffer, C_NULL_CHAR)
      if (n == 0) then
        message = "Cannot read remark"
        status = 1
        exit catch
      end if
      remark(i) = buffer(1:n-1)
    end do

    ! Read two dummy arguments (NOTE: Honestly, I have no idea what they do)
    read (unit, IOSTAT=status, IOMSG=message) dummy, dummy
    if (status /= 0) exit catch

    ! Read number of atoms and dummy arguments
    read (unit, IOSTAT=status, IOMSG=message) natoms, dummy
    if (status /= 0) exit catch

    ! Inquire current position and file size
    inquire (UNIT=unit, POS=pos, SIZE=filesize, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch
    pos = pos - 1

    ! Each frame has natoms*3 (4 bytes each) = natoms*12
    ! plus 6 box dimensions (8 bytes each) = 48
    ! Additionally there are 32 bytes of file information in each frame
    framesize = natoms * 12 + 80

    ! Header is typically 276 bytes, but inquire gives us exact size
    ! ERROR: "Excpected '", nframes, "' number of frames got '", nframes2, "'"
    nframes2 = (filesize - pos) / framesize
    if (nframes2 /= nframes) then
      message = "Extected and actual number of frames do not match"
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

end subroutine dcd_read_header

! Write DCD header. DCD files have only one header.
subroutine dcd_write_header (unit, remark, start_time, every_time, end_time, timestep, natoms, stat, errmsg)
  use iso_c_binding, only: C_NULL_CHAR
  implicit none
  integer, intent(in) :: unit
  character(*), intent(in) :: remark(:)
  integer, intent(in) :: start_time
  integer, intent(in) :: every_time
  integer, intent(in) :: end_time
  real, intent(in) :: timestep
  integer, intent(in) :: natoms
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  character(9) :: action, access
  character(80) :: buffer
  logical :: opened
  integer :: i, status

  catch: block

    ! Check if unit is assigned and opened for reading/writing.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READWRITE") == 0) then
      message = "File unit must be opened as readwrite"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch
  
    ! Write magic number and string
    write (unit, IOSTAT=status, IOMSG=message) MAGIC_NUMBER, MAGIC_STRING
    if (status /= 0) exit catch

    ! Write number of frames; This value is updated with each new frame write.
    write (unit, IOSTAT=status, IOMSG=message) 0 ! nframes
    if (status /= 0) exit catch

    ! Timestep of first, last snapshot, and snapshot stride.
    write (unit, IOSTAT=status, IOMSG=message) start_time, every_time, end_time
    if (status /= 0) exit catch
    
    ! Write five dummy arguments (don't know what they are).
    write (unit, IOSTAT=status, IOMSG=message) (0, i = 1, 5)
    if (status /= 0) exit catch

    ! Write simulation timestep.
    write (unit, IOSTAT=status, IOMSG=message) timestep
    if (status /= 0) exit catch

    ! Write nine dummy arguments: "has unit cell" + 8 zeros.
    write (unit, IOSTAT=status, IOMSG=message) 1, (0, i = 1, 8)
    if (status /= 0) exit catch
  
    ! We pretend to be CHARMM version 24 (don't know where the values come form).
    write (unit, IOSTAT=status, IOMSG=message) 24, 84, 164
    if (status /= 0) exit catch
    
    ! Write remarks and terminate them with C_NULL_CHAR
    write (unit, IOSTAT=status, IOMSG=message) size(remark)
    if (status /= 0) exit catch
    do i = 1, size(remark)
      write (buffer, "(a79,a)") remark(i), C_NULL_CHAR
      write (unit, IOSTAT=status, IOMSG=message) buffer
      if (status /= 0) exit catch
    end do

    ! Write number of atoms in snapshot and some more magic numbers
    write (unit, IOSTAT=status, IOMSG=message) 164, 4, natoms, 4
    if (status /= 0) exit catch

    ! Flush output
    flush (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine dcd_write_header

! Read number of atoms in tajectory from DCD file header (DCD files have only one header) and
! * return back to original position in file.
subroutine dcd_read_natoms (unit, natoms, stat, errmsg)
  use iso_fortran_env, only: INT64
  implicit none
  integer, intent(in) :: unit
  integer, intent(out) :: natoms
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  character(ERR_LEN) :: message
  logical :: opened
  character(9) :: action, access
  character(REMARK_LEN) :: remark
  integer :: i, nremarks, dummy, status
  integer(INT64) :: pos

  catch: block
  
    ! Check if file unit is opened for reading and store current position.
    inquire (UNIT=unit, OPENED=opened, POS=pos, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit must be opened for reading"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Read number of remarks and remark(s)
    read (unit, POS=97, IOSTAT=status, IOMSG=message) nremarks
    if (status /= 0) exit catch
    do i = 1, nremarks
      read (unit, IOSTAT=status, IOMSG=message) remark
      if (status /= 0) exit catch
    end do

    ! NOTE: Honestly, I have no idea what dummy parameters do. 
    read (unit, IOSTAT=status, IOMSG=message) dummy, dummy, natoms, dummy
    if (status /= 0) exit catch

    ! If original position is 1, that means header was not yet processed.
    ! * This means current position is beginning of first frame.
    if (pos /= 1) then
      read (unit, POS=pos-4, IOSTAT=status, IOMSG=message) dummy
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine dcd_read_natoms

! Read data in DCD format.
subroutine dcd_read_data (unit, natoms, coor, box, stat, errmsg)
  use iso_fortran_env, only: IOSTAT_END, REAL32, REAL64
  implicit none
  integer, intent(in) :: unit
  integer, intent(in) :: natoms
  real, intent(out) :: coor(DIM,natoms)
  real, intent(out) :: box(DIM,DIM)
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  real(REAL64) :: a, b, c, alpha, beta, gamma, k
  real :: x(natoms), y(natoms), z(natoms)
  character(ERR_LEN) :: message
  logical :: opened
  character(9) :: action, access
  integer :: i, status, dummy(6), magic, nbytes

  catch: block

    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File unit must be opened for reading"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch  

    ! Expected size in bytes
    nbytes = natoms * 4
  
    read (unit, IOSTAT=status, IOMSG=message) magic
    if (magic /= 48) then
      message = "Cannot read magic number or snapshot is corrupt"
      status = 1
    end if
    if (status /= 0) exit catch 
    
    ! Simulation box size is in format: a, gamma, b, beta, alpha, c. Box dimensions are double (REAL64).
    ! ERROR: "Cannot read simulation box vectors or snapshot is corrupt"
    read (unit, IOSTAT=status, IOMSG=message) a, gamma, b, beta, alpha, c
    if (status /= 0) exit catch 
    
    ! TODO: Apply box transform (See: http://gisaxs.com/index.php/Unit_cell)
    k = (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)
    ! box(:,1)= real([ a, 0d0, 0d0 ])
    ! box(:,2)= real([ b * cos(gamma), b * sin(gamma), 0d0 ])
    ! box(:,3)= real([ c * cos(beta), c * k, c * sqrt(1d0 - cos(beta)**2 - k**2) ])
    
    ! NOTE:For now we only support cubic simulation box
    box = 0.
    box(1,1) = real(a) 
    box(2,2) = real(b) 
    box(3,3) = real(c) 

    ! Read x, y, and z coordinates as follows:
    ! * magic number
    ! * num. of bytes, x-coordinates, num. of bytes
    ! * num. of bytes, y-coordinates, num. of bytes
    ! * num. of bytes, z-coordinates, num. of bytes
    read (unit, IOSTAT=status, IOMSG=message) magic, dummy(1), x, dummy(2:3), y, dummy(4:5), z, dummy(6)
    if (magic /= 48) then
      message = "Cannot read coordinate values or snapshot is corrupt"
      status = 1
    else if (any(dummy(:) /= nbytes)) then
      message = "Number of bytes in snapshot is incorrect for size of coor array passed"
      status = 1
    end if
    if (status /= 0) exit catch 
    
    ! Write back to array
    do i = 1, natoms
      coor(:,i) = [x(i), y(i), z(i)] 
    end do

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine dcd_read_data

! Skip DCD file data
subroutine dcd_skip_data (unit, natoms, stat, errmsg)
  use iso_fortran_env, only: INT64
  implicit none
  integer, intent(in) :: unit
  integer, intent(in) :: natoms
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  integer(INT64) :: pos, newpos, nbytes
  integer :: status, magic, dummy
  character(ERR_LEN) :: message
  character(9) :: action, access
  logical :: opened
  real(REAL64) :: box(6)

  catch: block

    ! Check if file unit is opened for reading.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, POS=pos, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READ") == 0) then
      message = "File must be opened for reading"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch   
  
    read (unit, IOSTAT=status, IOMSG=message) magic
    if (magic /= 48) then
      message = "Cannot read magic number or snapshot is corrupt"
      status = 1
    end if
    if (status /= 0) exit catch

    read (unit,IOSTAT=status) box
    if ( status /= 0 ) error stop

    ! Each frame has natoms*DIM*(4 bytes each) = natoms*DIM*4
    ! * plus 6 box dimensions*(8 bytes each) = 48
    ! * Additionally there are 32 bytes of file information in each frame
    nbytes = natoms * DIM * 4 + 48 + 32

    ! We subtract 4 bytes so that the next read of the 4 byte integer will
    ! * line things up properly for the next read.
    newpos = pos + nbytes - 4
    read (unit, POS=newpos, IOSTAT=status, IOMSG=message) dummy
    if (status /= 0) exit catch
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

  return
end subroutine dcd_skip_data

! Write data in DCD format. Write also updates header.
subroutine dcd_write_data (unit, natoms, coor, box, stat, errmsg)
  use iso_fortran_env, only: INT64, REAL64
  implicit none
  integer, intent(in) :: unit
  integer, intent(in) :: natoms
  real, intent(in) :: coor(DIM,natoms)
  real, intent(in) :: box(DIM,DIM)
  integer, intent(inout), optional :: stat ! Error status code
  character(*), intent(inout), optional :: errmsg ! Error message
  real(REAL64) :: a, b, c, alpha, beta, gamma
  integer(INT64) :: pos
  integer :: status, nframes, nbytes, start_time, every_time, end_time
  logical :: opened
  character(9) :: action, access
  character(ERR_LEN) :: message

  catch: block

    ! Check if unit is assigned and opened for reading/writing.
    inquire (UNIT=unit, OPENED=opened, ACTION=action, ACCESS=access, IOSTAT=status, IOMSG=message)
    if (.not. opened .or. index(action, "READWRITE") == 0) then
      message = "File unit must be opened as readwrite"
      status = 1
    else if (index(access, "STREAM") == 0) then
      message = "File must be opened as stream"
      status = 1
    end if
    if (status /= 0) exit catch

    ! Write magic number
    write (unit, IOSTAT=status, IOMSG=message) 48
    if (status /= 0) exit catch

    ! Box transformation. TODO: implement triclinic box
    a = dble(box(1,1))
    b = dble(box(2,2))
    c = dble(box(3,3))
    alpha = 90.0d0
    beta = 90.0d0
    gamma = 90.0d0  
    
    ! Write simulation box
    write (unit, IOSTAT=status, IOMSG=message) a, gamma, b, beta, alpha, c
    if (status /= 0) exit catch

    ! Calculate number of bytes for coordinates
    nbytes = natoms * 4
  
    ! Write x, y, and z coordinates as follows:
    ! * magic number
    ! * num. bytes, x-coordinates, num. bytes
    ! * num. bytes, y-coordinates, num. bytes
    ! * num. bytes, z-coordinates, num. bytes
    write (unit, IOSTAT=status, IOMSG=message) 48
    if (status /= 0) exit catch
    write (unit, IOSTAT=status, IOMSG=message) nbytes, coor(1,:), nbytes
    if (status /= 0) exit catch
    write (unit, IOSTAT=status, IOMSG=message) nbytes, coor(2,:), nbytes
    if (status /= 0) exit catch
    write (unit, IOSTAT=status, IOMSG=message) nbytes, coor(3,:), nbytes
    if (status /= 0) exit catch

    ! Inquire current stream position
    inquire (UNIT=unit, POS=pos, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Update num. of frames in header (read->update->write)
    read (unit, POS=9, IOSTAT=status, IOMSG=message) nframes, start_time, every_time, end_time
    if (status /= 0) exit catch
    write (unit, POS=9, IOSTAT=status, IOMSG=message) nframes + 1, start_time, every_time, end_time + every_time
    if (status /= 0) exit catch

    ! Return back to end of stream
    write (unit, POS=pos, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

    ! Flush output
    flush (unit, IOSTAT=status, IOMSG=message)
    if (status /= 0) exit catch

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine dcd_write_data

end module atomlib_dcdio
