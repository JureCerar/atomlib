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

module atomlib_trajio
  use iso_fortran_env, only: INT64
  use xdrfor
  use atomlib_xyzio
  use atomlib_groio
  use atomlib_pdbio
  use atomlib_crdio
  use atomlib_xtcio
  use atomlib_trrio
  use atomlib_dcdio
  use atomlib_frame
  implicit none
  private
  public :: traj_t
  
  integer, parameter :: DIM = 3 ! Default dimension.
  integer, parameter :: ERRLEN = 128 ! Length of error string.
  integer, parameter :: SEEK_SET = 0, SEEK_CUR = 1, SEEK_END = 2 ! Positioning constants.

  type :: traj_t
    type(xdrfile), pointer, private :: xd => null() ! XDR file pointer
    integer, private :: unit = 0 ! Fortran file unit
    character(3), private :: file_type = "" ! Current file (opened) type.
    integer, private :: file_num_frames = 0 ! Total number of frames in (opened) file.
    integer, private :: file_current_frame = 0 ! Current frame in (opened) file.
    integer(INT64), allocatable, private :: offset(:) ! Frame offset table
    ! Public variables ------------
    integer :: nframes = 0 ! Number of frames
    type(frame_t), allocatable :: frame(:) ! Frame data
  contains
    ! Allocation functions
    procedure :: allocate => traj_allocate
    procedure :: reallocate => traj_reallocate
    procedure :: deallocate => traj_deallocate
    procedure :: allocated => traj_allocated
    ! Formatted I/O
    procedure, private :: traj_write_formatted
    generic :: write(formatted) => traj_write_formatted
    ! File I/O
    procedure :: open => traj_open
    procedure :: close => traj_close
    procedure :: fseek => traj_fseek
    procedure :: ftell => traj_ftell
    procedure :: read_next => traj_read_next
    procedure :: save => traj_save
    procedure :: load => traj_load
  end type traj_t

  ! Constructor interface
  interface traj_t
    module procedure :: traj_constructor
  end interface traj_t

contains 

! Util: Returns extension of file (always in upper case).
function extension (file) result (out)
  implicit none
  character(:), allocatable :: out
  character(*), intent(in) :: file
  integer, parameter :: offset = ichar("a") - ichar("A")
  integer :: i

  i = index(file, ".", BACK=.true.)
  if (i == 0) then
    out = ""
  else
    out = file(i+1:len_trim(file))
  end if

  do i = 1, len_trim(out)
    select case (out(i:i))
    case ("a" : "z")
      out(i:i) = char(ichar(out(i:i)) - offset)
    end select
  end do

end function extension

! Constructor for traj_t dferived type.
type(traj_t) function traj_constructor (nframes, frame) result (this)
  implicit none
  integer, intent(in) :: nframes
  type(frame_t), intent(in), optional :: frame(nframes)

  call this%allocate(nframes)
  if (present(frame)) this%frame = frame

end function traj_constructor

! Allocate and initialize traj_t data for specified number of frames.
subroutine traj_allocate (this, nframes, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  integer, intent(in) :: nframes ! Number of frames to allocate
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status
  
  catch: block 

    if (nframes < 0) then
      message = "Invalid number of frames"
      status = 1
      exit catch
    end if
    
    if (allocated(this%frame)) deallocate (this%frame, STAT=status)
    allocate (this%frame(nframes), STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    this%nframes = nframes

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_allocate

! Reallocate traj_t data for specified number of frames (bigger or smaller).
subroutine traj_reallocate (this, nframes, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  integer, intent(in) :: nframes ! New number of frames.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  type(frame_t), allocatable :: temp(:)
  character(ERRLEN) :: message
  integer :: i, status

  catch: block 

    if (nframes <= 0) then
      message = "Invalid number of frames"
      status = 1
      exit catch
    end if

    this%nframes = nframes

    if (allocated(this%frame)) then
      allocate (temp(nframes), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      i = min(nframes, size(this%frame))
      temp(:i) = this%frame(:i)
      call move_alloc (temp, this%frame)
    else
      allocate (this%frame(nframes), STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_reallocate 

! Deallocate traj_t data.
subroutine traj_deallocate (this, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  integer, intent(out), optional :: stat
  character(*), intent(out), optional :: errmsg
  character(ERRLEN) :: message
  integer :: status

  catch: block

    status = 0

    if (allocated(this%frame)) deallocate (this%frame, STAT=status, ERRMSG=message) 
    if (status /= 0) exit catch

    this%nframes = 0  

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)  
  
end subroutine traj_deallocate 

! Check if traj_t is (properly) allocated.
function traj_allocated (this) result (out)
  implicit none
  logical :: out
  class(traj_t), intent(in) :: this

  out = allocated(this%frame)

end function traj_allocated 

! Formatted traj_t data.
subroutine traj_write_formatted (this, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(traj_t), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype 
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  character(80) :: fmt
  integer :: i

  catch: block

    if (.not. this%allocated()) then
      iomsg = "Object not allocated"
      iostat = 1
      exit catch
    end if 

    if (iotype == "LISTDIRECTED") then
      write (unit, "(i0,/)", IOSTAT=iostat, IOMSG=iomsg) this%nframes
      if (iostat /= 0) exit catch

      do i = 1, size(this%frame)
        write (unit, *, IOSTAT=iostat, IOMSG=iomsg) this%frame(i)
        if (iostat /= 0) exit catch
      end do
    
    else if (iotype == "DT") then
      if (size(v_list) == 0) then
        fmt = "(DT)"
      else if (size(v_list) == 2) then
        write (fmt,"(a,i0,a,i0,a)") "(DT(", v_list(1), ",", v_list(2), "))"
      else
        iomsg = "Integer-list for DT descriptor must contain two integers"
        iostat = 1
        exit catch
      end if
      
      write (unit, "(i0,/)", IOSTAT=iostat, IOMSG=iomsg)  this%nframes
      if (iostat /= 0) exit catch

      do i = 1, size(this%frame)
        write (unit, fmt, IOSTAT=iostat, IOMSG=iomsg) this%frame(i)
        if (iostat /= 0) exit catch
      end do

    else
      iomsg = "Unsupported iotype"
      iostat = 1
      exit catch

    end if

  end block catch

end subroutine traj_write_formatted

! Util: Read from unit/pointer in specified format.
subroutine frame_read_unit (unit, type, frame, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file unit.
  character(*), intent(in) :: type ! File type (in uppercase).
  type(frame_t), intent(inout) :: frame ! Frame data.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status, natoms
  
  catch: block

    select case (type)
    case ("XYZ")
      ! NOTE: XYZ default units are Angstroms.
      xyz: block
        character(80) :: comment
        call xyz_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate(natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call xyz_read_coor (unit, natoms, frame%coor, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        frame%coor = frame%coor * 0.1 ! to [nm].
        frame%box = frame%box * 0.1 ! to [nm].
      end block xyz
    
    case ("GRO")
      gro: block
        character(128) :: comment
        call gro_read_header (unit, comment, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call gro_read_coor (unit, natoms, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block gro

    case ("PDB")
      ! NOTE: PDB default units are Angstroms.
      pdb: block
        call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch      
        call pdb_read_coor (unit, natoms, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        frame%coor = frame%coor * 0.1 ! to [nm].
        frame%box = frame%box * 0.1 ! to [nm].
      end block pdb

    case ("CRD", "COR")
      ! NOTE: CRD/COR default units are Angstroms.
      crd: block
        logical :: extended
        call crd_read_header (unit, extended, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch    
        call crd_read_coor (unit, extended, natoms, frame%coor, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        frame%coor = frame%coor * 0.1 ! to [nm].
        frame%box = frame%box * 0.1 ! to [nm].
      end block crd

    case ("DCD")
      ! NOTE: DCD file default units are Angstroms.
      dcd: block
        call dcd_read_natoms (unit, natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call dcd_read_data (unit, natoms, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        frame%coor = frame%coor * 0.1 ! to [nm].
        frame%box = frame%box * 0.1 ! to [nm].
      end block dcd

    case default
      message = "Unsupported file extension: '" // trim(type) // "'"
      status = 1
      exit catch

    end select
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)
 
end subroutine frame_read_unit

subroutine frame_read_xdr (xd, type, frame, stat, errmsg)
  implicit none
  type(xdrfile), pointer, intent(in) :: xd ! XDR file pointer.
  character(*), intent(in) :: type ! File type (in uppercase).
  type(frame_t), intent(inout) :: frame ! Frame data.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: natoms, status

  catch: block
    select case (type)
    case ("XTC")
      xtc: block
        integer :: step
        real :: prec
        call xtc_do_header (xd, natoms, step, frame%time, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate (natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call xtc_do_data (xd, .true., natoms, prec, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block xtc

    case ("TRR")
      trr: block
        type(trnheader) :: sh
        call trr_do_header (xd, .true., sh, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call frame%allocate (sh%natoms, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call trr_do_data (xd, sh, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        frame%time = sh%time
      end block trr

    case default
      message = "Unsupported file extension: '" // trim(type) // "'"
      status = 1
      exit catch

    end select
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)
  
end subroutine frame_read_xdr

! Util: Write to unit/pointer in specified format
subroutine frame_write_unit (unit, type, frame, stat, errmsg)
  implicit none
  integer, intent(in) :: unit ! Fortran file handle.
  character(*), intent(in) :: type ! File type (in uppercase).
  class(frame_t), intent(inout) :: frame ! Frame data.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: i, natoms, status
  
  catch: block

    if (.not. frame%allocated()) then
      message = "Object not allocated"
      status = 1
      exit catch
    end if

    natoms = size(frame%coor, DIM=2)
    if (natoms <= 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if 

    select case (type)  
    case ("XYZ")
      ! NOTE: XYZ default units are Angstroms.
      xyz: block
        character(80) :: comment = ""
        call xyz_write (unit, comment, natoms, [("", i = 1, natoms)], frame%coor * 10, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block xyz

    case ("GRO")
      gro: block
        character(128) :: comment = ""
        call gro_write (unit, comment, natoms, [(i, i = 1, natoms)], [("", i = 1, natoms)], &
        & [(0, i = 1, natoms)], [("UNK", i = 1, natoms)], frame%coor, frame%box,  STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block gro

    case ("PDB")
       message = "Trajectory write in PDB format is not supported"
      status = 1
      exit catch

    case ("CRD", "COR")
      error stop

    case ("DCD")
      ! NOTE: DCD default units are Angstroms.
      dcd: block
        call dcd_write_data (unit, frame%natoms, frame%coor * 10, frame%box * 10, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block dcd 
      
    case default
      message = "Unsupported file extension: '" // trim(type) // "'"
      status = 1
      exit catch

    end select
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)
 
end subroutine frame_write_unit

subroutine frame_write_xdr (xd, type, frame, stat, errmsg)
  implicit none
  type(xdrfile), pointer, intent(in) :: xd ! XDR file pointer.
  character(*), intent(in) :: type ! File type (in uppercase).
  class(frame_t), intent(inout) :: frame ! Frame data.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: natoms, status
  
  catch: block

    if (.not. frame%allocated()) then
      message = "Object not allocated"
      status = 1
      exit catch
    end if

    natoms = size(frame%coor, DIM=2)
    if (natoms <= 0) then
      message = "Invalid number of particles"
      status = 1
      exit catch
    end if 

    select case (type)  
    case ("XTC")
      xtc: block
        real :: prec = 0.0
        integer :: step = 0
        call xtc_do_header (xd, natoms, step, frame%time, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call xtc_do_data (xd, .false., natoms, prec, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block xtc

    case ("TRR")
      trr: block
        use iso_fortran_env, only: REAL32
        type(trnheader) :: sh
        ! Generate custom header for writing
        sh%is_double = .false.
        sh%natoms = natoms
        sh%x_size = DIM * sh%natoms * storage_size(0.0_REAL32) / 8
        sh%box_size = DIM * DIM * storage_size(0.0_REAL32) / 8
        sh%time = frame%time
        call trr_do_header (xd, .false., sh, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
        call trr_do_data (xd, sh, frame%coor, frame%box, STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end block trr

    case default
      message = "Unsupported file extension: '" // trim(type) // "'"
      status = 1
      exit catch

    end select
  
  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)
 
end subroutine frame_write_xdr

! Open trajectory for reading. Check file and construct offset table.
subroutine traj_open (this, file, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  character(*), intent(in) :: file ! File name
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block
  
    call this%close (STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    this%file_type = extension(file)
    this%file_current_frame = 1

    select case (this%file_type)
    case ("XYZ")
      call xyz_open (file, this%unit, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case ("GRO")
      call gro_open (file, this%unit, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case ("PDB")
      call pdb_open (file, this%unit, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case ("CRD", "COR")
      call crd_open (file, this%unit, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case ("XTC")
      call xtc_open (file, this%xd, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch 

    case ("TRR")
      call trr_open (file, this%xd, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case ("DCD")
      call dcd_open (file, this%unit, this%file_num_frames, this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case default
      message = "Unsupported file type: '" // extension(file) // "'"
      status = 1
      exit catch

    end select

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_open 

! Close trajectory.
subroutine traj_close (this, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status

  catch: block

    select case (this%file_type)
    case ("XYZ", "GRO", "PDB", "CRD", "COR", "DCD")
      close (this%unit, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

    case ("XTC", "TRR")
      call xdr_close (this%xd, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

    case default
      status = 0

    end select

    this%file_num_frames = 0
    this%file_current_frame = 0
    this%file_type = ""

    if (allocated(this%offset)) then
      deallocate (this%offset, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
    end if 

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_close 

! Move position to the specified frame from whence.
subroutine traj_fseek (this, offset, whence, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  integer, intent(in) :: offset ! Position offset
  integer, intent(in) :: whence ! Origin of offset: beginning, current, end 
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: status, frame

  catch: block

    select case (whence)  
    case (SEEK_SET) ! Absolute position
      frame = offset
    case (SEEK_CUR) ! From current frame
      frame = this%file_current_frame + offset
    case (SEEK_END) ! From last frame
      frame = this%file_num_frames + offset
    case default
      message = "Unsupported whence position"
      status = 1
      exit catch
    end select

    frame = max(frame, 1)
    frame = min(frame, this%file_num_frames)

    select case (this%file_type)
    case ("XYZ", "GRO", "PDB", "CRD", "COR")
      ! Return to original position using dummy stream read i.e. read 4 bytes before 
      ! * the selected position. If position is at beggining of the file (<4 bytes) 
      ! * then simply rewind the file.  
      if (this%offset(frame) < 4) then
        rewind (this%unit, IOSTAT=status, IOMSG=message)
      else
        read (this%unit, *, POS=this%offset(frame)-4, IOSTAT=status, IOMSG=message)
      end if
      if (status /= 0) exit catch
      
    case ("XTC", "TRR")
      status = xdr_seek(this%xd, this%offset(frame), SEEK_SET)
      if (status /= 0) then
        message = "Cannot set position in XDR file"
        exit catch
      end if

    case ("DCD")
      read (this%unit, POS=this%offset(frame), IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

    case default
      message = "No file currently opened"
      status = 1
      exit catch

    end select

    this%file_current_frame = frame 

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_fseek

! Returns current position in trajectory file.
function traj_ftell (this) result (out)
  implicit none
  integer :: out
  class(traj_t), intent(inout) :: this

  out = this%file_current_frame

end function traj_ftell

! Read next frame (or next N frames) in trajectory file. 
subroutine traj_read_next (this, nframes, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  integer, intent(in), optional :: nframes ! Number of frames to read. Default is 1.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: i, num, status
  
  catch: block

    ! Check how many frames are left i.e. how many can we actually read.
    num = merge(nframes, 1, present(nframes))
    num = min(num, this%file_num_frames - this%file_current_frame + 1)
    if (num == 0) then
      message = "No frames were read"
      status = -1
      exit catch

    else if (num < 0) then
      message = "Invalid number of frames"
      status = 1
      exit catch

    end if

    ! Allocate specified number of frames.
    call this%allocate (num, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch
    
    ! Read according to data type.
    select case (this%file_type)
    case ("XYZ", "GRO", "PDB", "CRD", "COR", "DCD")
      do i = 1, num
        call frame_read_unit (this%unit, this%file_type, this%frame(i), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end do

    case ("XTC", "TRR")
      do i = 1, num
        call frame_read_xdr (this%xd, this%file_type, this%frame(i), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end do

    case default
      message = "No file currently opened"
      status = 1
      exit catch

    end select

    ! Update current frame.
    this%file_current_frame = this%file_current_frame + num

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_read_next

! Save (write) trajectory file. Format is based on file extension.  
subroutine traj_save (this, file, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  character(*), intent(in) :: file ! Input file.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(128) :: message
  character(3) :: type
  integer :: i, natoms, status

  catch: block

    if (.not. this%allocated()) then
      message = "Object not allocated"
      status = 1
      exit catch
    end if

    type = extension(file)
    select case (type)
    case ("XYZ", "GRO", "PDB", "CRD", "COR")
      open (NEWUNIT=this%unit, FILE=file, ACCESS="stream", FORM="formatted", &
      & STATUS="unknown", ACTION="write", IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

      do i = 1, size(this%frame)
        call frame_write_unit (this%unit, type, this%frame(i), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end do

      close (this%unit, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

    case ("DCD")
      open (NEWUNIT=this%unit, FILE=file, FORM="unformatted", ACCESS="stream", & 
      & ACTION="readwrite", STATUS="unknown", IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch
      
      ! DCD trajectory has fixed frame size.
      natoms = maxval(this%frame(:)%natoms)
      if (any(natoms /= this%frame(:)%natoms)) then
        message = "DCD files do not support multiple frame sizes."
        status = 1
        exit catch
      end if

      ! DCD files have single header.
      call dcd_write_header (this%unit, [""], 0, 0, 0, 0., natoms, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch

      do i = 1, size(this%frame)
        call frame_write_unit (this%unit, type, this%frame(i), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end do

      close (this%unit, IOSTAT=status, IOMSG=message)
      if (status /= 0) exit catch

    case ("XTC", "TRR")
      call xdr_open(this%xd, file, .false., STAT=status, ERRMSG=message)
      if (status /= exdrOK) exit catch

      ! XTC/TRR trajectory have fixed frame sizes.
      natoms = maxval(this%frame(:)%natoms)
      if (any(natoms /= this%frame(:)%natoms)) then
        message = "XTC/TRR files do not support multiple frame sizes."
        status = 1
        exit catch
      end if

      do i = 1, size(this%frame)
        call frame_write_xdr (this%xd, type, this%frame(i), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      end do

      call xdr_close (this%xd, STAT=status, ERRMSG=message)
      if (status /= exdrOK) exit catch

    case default
      message = "Unsupported file extension: '" // trim(type) // "'"
      status = 1
      exit catch

    end select

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)
  
end subroutine traj_save 

! Load (read) trajectory file.
subroutine traj_load (this, file, first, last, stride, stat, errmsg)
  implicit none
  class(traj_t), intent(inout) :: this
  character(*), intent(in) :: file ! File name
  integer, intent(in), optional :: first ! First frame in trajectory to read. Default is 1.
  integer, intent(in), optional :: last ! Last frame in trajectory to read. default is -1 (last).
  integer, intent(in), optional :: stride ! Read only n-th frame.
  integer, intent(out), optional :: stat ! Error status code. Returns zero if no error.
  character(*), intent(out), optional :: errmsg ! Error message.
  character(ERRLEN) :: message
  integer :: fmin, fmax, step
  integer :: i, next, status, nframes

  catch: block

    call this%open (file, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    ! Collect optional parameters
    fmin = merge(first, 1, present(first))
    fmax = merge(last, -1, present(last))
    step = merge(stride, 1, present(stride))

    ! Correct values if they do not make sense.
    fmin = max(1, fmin)
    fmax = merge(fmax, this%file_num_frames, fmax /= -1)
    fmax = max(1, min(this%file_num_frames, fmax))

    ! Calculate number of frames that can be read and allocate data
    nframes = ceiling(real(fmax - fmin + 1) / step)
    call this%allocate(nframes, STAT=status, ERRMSG=message)
    if (status /= 0) exit catch

    next = 1
    do i = fmin, fmax, step
      ! Go to selected frame
      call this%fseek (i, SEEK_SET, STAT=status, ERRMSG=message)
      if (status /= 0) exit catch
      
      if (next > nframes) then
        message = "Something went wrong: number of frames missmatch."
        status = 1
        exit catch
      end if

      ! Read according to file type
      select case (this%file_type)
      case ("XYZ", "GRO", "PDB", "CRD", "COR", "DCD")
        call frame_read_unit (this%unit, this%file_type, this%frame(next), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      case ("XTC", "TRR")
        call frame_read_xdr (this%xd, this%file_type, this%frame(next), STAT=status, ERRMSG=message)
        if (status /= 0) exit catch
      case default
        message = "Unsupported file extension: '" // extension(this%file_type) // "'"
        status = 1
        exit catch
      end select

      ! Increment frame counter
      next = next + 1
 
    end do ! for i

  end block catch

  if (present(stat)) then
    stat = status
  else if (status /= 0) then
    error stop message
  end if
  if (present(errmsg)) errmsg = trim(message)

end subroutine traj_load

end module atomlib_trajio

! program main
!   use atomlib_traj
!   implicit none
!   type(traj_t) :: obj
!   type(frame_t) :: frame
!   integer :: i, unit

!   ! test object
!   ! open (NEWUNIT=unit, FILE="../test/conf.xyz", STATUS="old")
!   ! call frame%read (unit, "XYZ")
!   ! print *, frame
!   ! close (unit)


!   ! call obj%allocate(10)
!   ! obj%coor = reshape([(i, i = 1, 3*10)], shape(obj%coor))
!   ! obj%box = reshape([(i, i = 1, 3*3)], shape(obj%box))
!   ! write (*,*) obj%frame(1)
!   ! obj%frame(1)%coor = reshape([(i, i = 1, 3*10)], shape(obj%coor))
!   ! obj%frame(1)%box = reshape([(i, i = 1, 3*3)], shape(obj%box))
!   ! write (*,*) obj%frame(1)

!   ! call obj%open("../test/conf.xyz")
!   ! call obj%skip_next (1)
!   ! call obj%read_next (2)
!   ! write (*, *) obj

!   ! call obj%read("../test/conf.xyz")
!   call obj%load("../test/conf.gro")
!   ! call obj%skip_next (1)
!   ! call obj%read_next (2)
!   write (*, *) obj

!   write (*, *) obj%frame(1)
!   write (*, "(DT(20,10))") obj%frame(1)

! end program main
