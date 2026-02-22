# pdbio

The Protein Data Bank (pdb) file format is a textual file format describing the three-dimensional structures of molecules held in the Protein Data Bank. The pdb format accordingly provides for description and annotation of protein and nucleic acid structures including atomic coordinates, secondary structure assignments, as well as atomic connectivity. In addition experimental metadata are stored. The units are in Ångströms.

```
TITLE     Two water molecules, t= 0.0
REMARK    THIS IS A SIMULATION BOX
CRYST1   18.206   18.206   18.206  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  OW1 WAT A   1       1.260  16.240  16.790  1.00  0.00           O  
ATOM      2  HW2 WAT A   1       1.900  16.610  17.470  1.00  0.00           H  
ATOM      3  HW3 WAT A   1       1.770  15.680  16.130  1.00  0.00           H  
ATOM      4  OW1 WAT A   2      12.750   0.530   6.220  1.00  0.00           O  
ATOM      5  HW2 WAT A   2      13.370   0.020   6.800  1.00  0.00           H  
ATOM      6  HW3 WAT A   2      13.260   1.200   5.680  1.00  0.00           H  
TER
ENDMDL
```

PDB files contais all sort of metadata we are only interested with keyword`ATOM` or `HETATM` tha contain atom infromation. Lines contain the following information:
- Record type, 
- Atom serial number,
- Atom name,
- Alternate location indicator,
- Residue name.
- Chain identifier,
- Residue sequence number, 
- Code for insertions of residues,
- Orthogonal coordinates for `x`, `y`, and `z`,
- Occupancy,
- Temperature factor,
- Segment identifier,
- Element symbol,
- and atoms charge

Fortran format: `(6x,i5,1x,a4,a1,a3,x,a1,i4,a1,3x,3f8.3,2f6.2,10x,a2,i2)`

For more information please see [reference](https://www.wwpdb.org/documentation/file-format).

---------------------------------------------------------------------
## `pdb_open`
Open file for reading and checked. Frame offset table is constructed (byte position where frame starts). Move to selected frame with [fseek](https://gcc.gnu.org/onlinedocs/gfortran/FSEEK.html#FSEEK). If the function fails stat is a positive number, otherwise it is set to zero. If stat is not present end error occurs error stop is raised.
```fortran
call pdb_open (file, unit, nframes, offset, stat, errmsg)
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

## `pdb_close`
Close opened file.
```fortran
call pdb_close (unit, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `pdb_get_natoms`
Count num. of atoms in current PDB frame i.e. count occurance of "ATOM" and "HETATOM" keywords. File must be opened as Fortran 'formatted stream'.
```
call pdb_get_natoms (unit, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(OUT) :: natoms`
> Number of atoms.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `pdb_read_data`
Read data in PDB format.
```fortran
call pdb_read_data (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
& altloc, qfact, bfact, element, charge, coor, box, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, dimension(natoms), intent(OUT) :: atomi, resi`
> Atom and residue index.

`character(*), dimension(natoms), intent(OUT) :: atomn, resn`
> Atom and residue name.

`character(*), intent(OUT) :: chain`
> Chain identifier.

`character(*), intent(OUT) :: resic`
> Code for insertions of residues

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(natoms), intent(OUT) :: qfact, bfact`
> Occupancy and temperature factor.

`character(*), dimension(natoms), intent(OUT) :: element`
> Element symbol.

`integer, dimension(natoms), intent(OUT) :: charge`
> Elemental charge.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `pdb_read_coor`
Read ONLY coordinate data in PDB format.
```fortran
call pdb_read_coor (unit, natoms, coor, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`real, dimension(DIM,natoms), intent(OUT) :: coor`
> Atom cartesian coordinates in Å.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `pdb_skip_data`
Read through file (skip) in PDB format but do not store any data.
```fortran
call pdb_skip_data (unit, natoms, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## `pdb_write`
Read data in PDB format.
```fortran
call pdb_write (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
& altloc, qfact, bfact, element, charge, coor, box, stat, errmsg)
```
`integer, intent(IN) :: unit`
> Fortran file unit.

`integer, intent(IN) :: natoms`
> Number of atoms.

`integer, dimension(natoms), intent(IN) :: atomi, resi`
> Atom and residue index.

`character(*), dimension(natoms), intent(IN) :: atomn, resn`
> Atom and residue name.

`character(*), intent(IN) :: chain`
> Chain identifier.

`character(*), intent(IN) :: resic`
> Code for insertions of residues

`real, dimension(DIM,natoms), intent(IN) :: coor`
> Atom cartesian coordinates in Å.

`real, dimension(natoms), intent(IN) :: qfact, bfact`
> Occupancy and temperature factor.

`character(*), dimension(natoms), intent(IN) :: element`
> Element symbol.

`integer, dimension(natoms), intent(IN) :: charge`
> Elemental charge.

`integer, intent(IN), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(IN), optional :: errmsg`
> Error status message. Empty if no error.
---------------------------------------------------------------------

## Example
Example of reading PDB file.
```fortran
program main
  use atomlib_pdbio
  use iso_fortran_env, only: INT64
  implicit none
  integer, parameter :: DIM = 3
  integer(INT64), allocatable :: offset(:)
  integer, allocatable :: atomi(:), resi(:), charge(:)
  character(6), allocatable :: atomn(:), resn(:), chain(:), resic(:), altloc(:), element(:)
  real, allocatable :: qfact(:), bfact(:), coor(:,:)
  real :: box(DIM,DIM)
  integer :: i, unit, nframes, natoms, status

    call pdb_open (arg, unit, nframes, offset, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

    call pdb_get_natoms (unit, natoms, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

    allocate (atomi(natoms), atomn(natoms), resi(natoms), resn(natoms), chain(natoms), resic(natoms), &
    & altloc(natoms), qfact(natoms), bfact(natoms), element(natoms), charge(natoms), coor(DIM,natoms))

    call pdb_read_data (unit, natoms, atomi, atomn, resi, resn, chain, resic, &
    & altloc, qfact, bfact, element, charge, coor, box, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message

    call pdb_close (unit, STAT=status, ERRMSG=message)
    if (status /= 0) error stop message 

end program main
```