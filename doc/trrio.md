# trrio

Gromacs portable binary trajectory format. File contains all the coordinates, velocities, forces and energies. 

**NOTE:** xrdfile *read* is same as *write*. He he he.

---------------------------------------------------------------------

## `trnheader`

```fotran
type trnheader
  integer :: ir_size, e_size, vir_size, pres_size, top_size, sym_size, nre
  integer :: box_size, x_size, v_size, v_size 
  logical :: is_double
  integer :: natoms, step
  real :: time, lambda
end type trnheader
```

---------------------------------------------------------------------

## `trr_open`
Open file for reading and check for errors. Frame offset table is constructed (byte position where frame starts). Move to selected frame with [xdr_seek](xdrfor.md#xdr_seek). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call trr_open (file, xd, nframes, offset, stat, errmsg)
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

## `trr_close`
Close opened file.
```fortran
call trr_close (xd, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(INOUT) :: xd`
> XDR file pointer.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `trr_do_header`
Read/write file header information in TRR format.
```fortran
call trr_do_header (xd, bRead, sh, stat, errmsg)
```
`type(XDRFILE), pointer, intent(IN) :: xd`
> XDR file pointer.

`logical, intent(IN) :: bRead`
> Read header information? Else write header.

`type(TRNHEADER), intent(INOUT) :: sh`
> TRR file header information.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `trr_do_data`
Read/write data in TRR format. Wrapper for `trr_do_real` and `trr_do_double`.
```fortran
call trr_do_data (xd, sh, coor, box, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(IN) :: xd`
> XDR file pointer.

`type(TRNHEADER), intent(INOUT) :: sh`
> TRR file header information.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in nm.

`real, dimension(DIM,DIM), intent(OUT) :: box`
> Simulation box dimensions in nm.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

#### *Notes*
Loss of data if used on double precision trajectory. Use `trr_do_double`. To get velocity, and force information use `trr_do_real` and `trr_do_double`.

---------------------------------------------------------------------

## `trr_do_skip`
Read through file (skip) in TRR format but do not store any data.
```fortran
call trr_do_skip (xd, sh, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(IN) :: xd`
> XDR file pointer.

`type(TRNHEADER), intent(INOUT) :: sh`
> TRR file header information.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## Example
```fortran
program main
  use xdrfor
  use atomlib_trrio
  use iso_fortran_env, only: INT64
  implicit none
  integer, paramete
  type(xdrfile), pointer :: xd => null()
  integer(INT64), allocatable :: offset(:)
  type(trnheader) :: sh
  real, allocatable :: coor(:,:)
  real :: box(DIM,DIM)
  integer :: nframes, status

  ! Read
  call trr_open ("traj.trr", xd, nframes, offset, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call trr_do_header (xd, .True., sh, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  allocate (coor(DIM,sh%natoms), STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call trr_do_data (xd, sh, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call trr_close (unit, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message  

  ! Write
  call xdr_open (xd, "out.trr", .False., STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call trr_do_header (xd, .False., sh, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call trr_do_data (xd, sh, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop messa

  call xdr_close (xd, STAT=status, ERRMSG=message)
  if (status /= 0) error stop messa

end program main
```
