# xdrfor

This is a fortran wrapper for `xdrfile.c`. **NOTE:** xrdfile *read* is same as *write*. He he he.

---------------------------------------------------------------------

## `xdrfile_read_int`
Read/write one or more integer type variable(s). Returns number of integers read/written. If this is negative, an error occured.
```fortran
out = xdrfile_read_int(ptr, ndata, xfp)
```
#### *Parameters*
`integer(C_INT), dimension(..), intent(INOUT) :: ptr`
> Pointer to memory where data should be read/written.

`integer(C_INT), value, intent(IN) :: ndata`
> Number of integers to write/read.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of values written/read. Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_read_float`
Read/write one or more float type variable(s). Returns number of integers read/written. If this is negative, an error occured.
```fortran
out = xdrfile_read_float(ptr, ndata, xfp)
```
#### *Parameters*
`real(C_FLOAT), dimension(..), intent(INOUT) :: ptr`
> Pointer to memory where data should be read/written.

`integer(C_INT), value, intent(IN) :: ndata`
> Number of integers to write/read.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of values written/read. Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_read_double`
Read/write one or more double type variable(s). Returns number of integers read/written. If this is negative, an error occured.
```fortran
out = xdrfile_read_double(ptr, ndata, xfp)
```
#### *Parameters*
`real(C_DOUBLE), dimension(..), intent(INOUT) :: ptr`
> Pointer to memory where data should be read/written.

`integer(C_INT), value, intent(IN) :: ndata`
> Number of integers to write/read.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of values written/read. Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_decompress_coord_float`
Decompress coordiates from XDR file to array of floats. Returns number of coordinate triplets read/written. If this is negative, an error occured.
```fortran
out = xdrfile_decompress_coord_float(ptr, ncoord, precision, xfp)
```
#### *Parameters*
`real(C_FLOAT), dimension(..), intent(INOUT) :: ptr`
> Pointer to coordinates to compress (length >= 3*ncoord).

`integer(C_INT), intent(IN) :: ncoord`
> Max number of coordinate triplets to read on input.

`real(C_FLOAT), intent(OUT) :: precision`
> The precision used in the compression.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of values read. Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_compress_coord_float`
Compresses coordiates to XDRFILE to array of floats. Returns number of coordinate triplets written. If this is negative, an error occured.
```fortran
out = xdrfile_compress_coord_float(ptr, ncoord, precision, xfp)
```
#### *Parameters*
`real(C_FLOAT), dimension(..), intent(IN) :: ptr`
> Pointer to coordinates to compress (length >= 3*ncoord).

`integer(C_INT), intent(IN) :: ncoord`
> Max number of coordinate triplets to read on input.

`real(C_FLOAT), intent(IN) :: precision`
> The precision used in the compression.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of values written. Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_read_string`
Read string type (array of characters) variable. If no end-of-string is encountered, one byte less than this is read and end-of-string appended. Returns number of characters read, including end-of-string.
```fortran
out = xdrfile_read_string(ptr, maxlen, xfp)
```
#### *Parameters*
`character(C_CHAR), dimension(*), intent(OUT) :: ptr`
> Pointer to memory where data should be written.

`integer(C_INT), value, intent(IN) :: maxlen`
> Maximum length of string.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of characters read. Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_write_string`
Write string type (array of characters) variable. Returns number of characters written, including end-of-string.
```fortran
out = xdrfile_write_string(ptr, xfp)
```
#### *Parameters*
`character(C_CHAR), dimension(*), intent(IN) :: ptr`
> Pointer to memory where data should be written.

`type(XDRFILE), intent(IN) :: xfp`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Number of characters written. Returns non-zero if error.

---------------------------------------------------------------------

## `xdr_tell`
Returns absolute current position in XDRFILE.
```fortran
out = xdr_tell(xd)
```
#### *Parameters*
`type(XDRFILE), intent(IN) :: xd`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT64_T) :: out`
> Absolute current position in XDRFILE.

---------------------------------------------------------------------

## `xdr_seek`
Moves XDRFILE position to the specified offset from whence. **NOTE:** This uses C based indexing (0, 1, 2, ...).
```fortran
out = xdr_seek(xd, offset, whence)
```
#### *Parameters*
`type(XDRFILE), intent(IN) :: xd`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT64_T), value, intent(IN) :: offset`
> Position offset

`integer(C_INT), value, intent(IN) :: whence`
> Origin of offset: `0 = SEEK_SET`, `1 = SEEK_CUR`, `2 = SEEK_END`

`integer(C_INT) :: out`
> Returns non-zero if error.

---------------------------------------------------------------------

## `xdrfile_open`
Open a portable binary file, just like `fopen`. Returns pointer to abstract xdr file datatype, or `NULL` if an error occurs.
```fortran
out = xdrfile_open(file, mode)
```
#### *Parameters*
`character(C_CHAR), dimension(*), intent(IN) :: file`
> Full or relative path (including name) of the file

`character(C_CHAR), dimension(*), intent(IN) :: mode`
> `r` for reading, `w` for writing, and `a` for append.

`type(C_PTR) :: out` 
> Pointer to an abstract XDRFILE datatype.

---------------------------------------------------------------------

## `xdrfile_close`
Close a previously opened portable binary file, just like `fclose`. Returns non-zero on error.
```fortran
out = xdrfile_close(xd)
```
#### *Parameters*
`type(XDRFILE), intent(IN) :: xd`
> Pointer to an abstract XDRFILE datatype.

`integer(C_INT) :: out`
> Returns non-zero if error.

---------------------------------------------------------------------

## `xdr_open`
Wrapper for opening XDR files in Fortran. Use instead of `xdrfile_open`.
```fortran
call xdr_open (xd, file, bRead, stat, errmsg)
```
`type(XDRFILE), pointer, intent(OUT) :: xd`
> Pointer to an abstract XDRFILE datatype.

`character(*), intent(IN) :: file`
> Full or relative path of the file.

`ogical, intent(IN) :: bRead`
> Open in read mode? Otherwise in write mode.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `xdr_close`
Wrapper for closing XDR files in Fortran.
```fortran
call xdr_close(xd, stat, errmsg)
```
#### *Parameters*
`type(XDRFILE), pointer, intent(IN) :: xd`
> Pointer to an abstract XDRFILE datatype.

`integer, intent(OUT), optional :: stat`
> Error status code. Returns zero if no error.

`character(*), intent(OUT), optional :: errmsg`
> Error status message. Empty if no error.

---------------------------------------------------------------------

## `exdr_message`
Get error message for xdrfile error code.
```fortran
out = exdr_message(status)
```
#### *Parameters*
`integer, intent(IN) :: status`
> Error status code.

`character(:), allocatable :: out`
> Decoded error message.