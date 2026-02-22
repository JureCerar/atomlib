# xtcio

The xtc format is a portable format for trajectories. It uses the xdr routines for writing and reading data which was created for the Unix NFS system. The trajectories are written using a reduced precision algorithm which works in the following way: the coordinates (in nm) are multiplied by a scale factor, typically 1000, so that you have coordinates in pm. These are rounded to integer values. Then several other tricks are performed, for instance making use of the fact that atoms close in sequence are usually close in space too (e.g. a water molecule). To this end, the xdr library is extended with a special routine to write 3-D float coordinates. 

**NOTE:** xrdfile *read* is same as *write*. He he he.

---------------------------------------------------------------------

## `xtc_open`
Open file for reading and check for errors. Frame offset table is constructed (byte position where frame starts). Move to selected frame with [xdr_seek](xdrfor.md#xdr_seek). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call xtc_open (file, xd, nframes, offset, stat, errmsg)
```
#### *Parameters*
`character(*), intent(IN) :: file`
> Input file name.

` type(XDRFILE), pointer, intent(INOUT) :: xd`
> XDR file pointer.

`integer, intent(OUT) :: nframes`
> Number of frames in file.

`integer(INT64), allocatable, dimension(:), intent(INOUT) :: offset`
> Frame offset table. Array size same as nframes.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `xtc_close`
Close opened file.
```fortran
call xtc_close (xd, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(INOUT) :: xd`
> XDR file pointer.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `xtc_do_header`
Read/write file header information in XTC format.
```fortran
call xtc_do_header (xd, natoms, step, time, stat, errmsg)
```
`type(XDRFILE), pointer, intent(IN) :: xd`
> XDR file pointer.

`integer, intent(INOUT) :: natoms`
> Number of atoms in current frame.

`integer, intent(INOUT) :: step`
> Current simulation step number.

`integer, intent(INOUT) :: time`
> Time of current frame in ps.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `xtc_do_data`
Read/write data in XTC format.
```fortran
call xtc_do_data (xd, bRead, natoms, prec, coor, box, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(IN) :: xd`
> XDR file pointer.

`logical, intent(IN) :: bRead`
> Read data? Else write it.

`integer, intent(INOUT) :: natoms`
> Number of atoms in current frame.

`real, intent(INOUT) :: prec`
> Coordinate precission.

`real, dimension(DIM,natoms), intent(INOUT) :: coor`
> Atom cartesian coordinates in nm.

`real, dimension(DIM,DIM), intent(INOUT) :: box`
> Simulation box dimensions in nm.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `xtc_do_skip`
Read through file (skip) in XTC format but do not store any data.
```fortran
call xtc_do_skip (xd, natoms, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(IN) :: xd`
> XDR file pointer.

`integer, intent(IN) :: natoms`
> number of atoms in current frame.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## Example
```fortran
program main
  use xdrfor
  use atomlib_xtcio
  use iso_fortran_env, only: INT64
  implicit none
  integer, paramete
  type(xdrfile), pointer :: xd => null()
  integer(INT64), allocatable :: offset(:)
  real, allocatable :: coor(:,:)
  real :: box(DIM,DIM)
  real :: time, precision
  integer :: natoms, nframes, step, status

  ! Read
  call xtc_open ("traj.xtc", xd, nframes, offset, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xtc_do_header (xd, natoms, step, time, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  allocate (coor(DIM,natoms), STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xtc_do_data (xd, .True., natoms, prec, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xtc_close (unit, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message  

  ! Write
  call xdr_open (xd, "out.xtc", .False., STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xtc_do_header (xd, natoms, step, time, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xtc_do_data (xd, .False., natoms, prec, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call xdr_close (xd, STAT=status, ERRMSG=message)
  if (status /= 0) error stop messa

end program main
```

