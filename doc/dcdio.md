# dcdio

DCD is a binary trajectory file used by NAMD, CHARMM and LAMMPS as the default trajectory format.

---------------------------------------------------------------------

## `dcd_open`
Open file for reading and checked (DCD files have some convoluted file check procedures). Frame offset table is constructed (byte position where frame starts). Move to selected frame with [fseek](https://gcc.gnu.org/onlinedocs/gfortran/FSEEK.html#FSEEK). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call dcd_open (file, unit, nframes, offset, stat, errmsg)
```
#### *Parameters*
`character(*), intent(IN) :: file`
> Input file name.

`integer, intent(OUT) :: unit`
> Fortran file unit.

`integer, intent(OUT) :: nframes`
> Number of frames in file.

`integer(INT64), allocatable, dimension(:), intent(INOUT) :: offset`
> Frame offset table. Array size same as nframes.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_close`
Close opened file.
```fortran
call dcd_close (unit, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_read_header`
Read DCD file header. **NOTE:** DCD files have only one header.
```fortran
call dcd_read_header (unit, remark, start_time, every_time, end_time, timestep, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`character(*), allocatable, dimension(:), intent(OUT) :: remark`
> Header remark.

`integer, intent(OUT) :: nframes`
> Number of frames in trajectory.

`integer, intent(OUT) :: start_time, every_time, end_time`
> Trajectory staring time, time step, and end time 

`real, intent(OUT) :: timestep`
> Simulation timestepo in ps. Try start_time * timestep.

`integer, intent(OUT) :: natoms`
> Number of atoms in trajectory.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_write_header`
Write DCD file header. **NOTE:** DCD files have only one header.
```fortran
call dcd_write_header (unit, remark, start_time, every_time, end_time, timestep, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`character(*), dimension(:), intent(IN) :: remark`
> Header remark.

`integer, intent(IN) :: nframes`
> Number of frames in trajectory.

`integer, intent(IN) :: start_time, every_time, end_time`
> Trajectory staring time, time step, and end time 

`real, intent(IN) :: timestep`
> Simulation timestepo in ps. Try start_time * timestep.

`integer, intent(IN) :: natoms`
> Number of atoms in trajectory.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_read_natoms`
Read number of atoms in tajectory from DCD file header (DCD files have only one header) and return back to original position in file. Useful when header was already read (in another subroutine per say).
```fortran
call dcd_read_natoms (unit, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT) :: natoms`
> Number of atoms in trajectory.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_read_data`
Read data in DCD format.
```fortran
call dcd_read_data (unit, natoms, coor, box, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms in trajectory.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(DIM,DIM), intent(OUT) :: box`
> Simulation box dimensions in Å.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_skip_data`
Read through file (skip) in DCD format but do not store any data.
```fortran
call dcd_skip_data (unit, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms in trajectory.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `dcd_write_data`
Write data in DCD format. Write also updates header info (nframes, etc.).
```fortran
call dcd_write_data (unit, natoms, coor, box, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms in trajectory.

`real, dimension(DIM,natoms), intent(IN) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(DIM,DIM), intent(IN) :: box`
> Simulation box dimensions in Å.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## Example
```fortran
program main
  use atomlib_dcdio
  use iso_fortran_env, only: INT64
  implicit none
  integer, paramete
  integer(INT64), allocatable :: offset(:)
  real, allocatable :: coor(:,:)
  real :: timestep, box(DIM,DIM)
  integer :: start_time, every_time, end_time
  integer :: natoms, nframes, unit, status

  ! Read
  call dcd_open ("traj.dcd", unit, nframes, offset, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call dcd_read_header (unit, nframes, remark, start_time, every_time, end_time, &
  & timestep, natoms, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  allocate (coor(DIM,natoms), STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call dcd_read_data (unit, natoms, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call dcd_close (unit, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message  

  ! Write
  open (NEWUNIT=unit, FILE="out.dcd", FORM="unformatted", ACCESS="stream", & 
  & ACTION="readwrite", STATUS="unknown", IOSTAT=status, IOMSG=message)
  if (status /= 0) error stop message

  call dcd_write_header (unit, remark, start_time, every_time, end_time, &
  & timestep, natoms, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  call dcd_write_data (unit, natoms, coor, box, STAT=status, ERRMSG=message)
  if (status /= 0) error stop message

  close (unit, IOSTAT=status, IOMSG=message)
  if (status /= 0) error stop messa

end program main
```
