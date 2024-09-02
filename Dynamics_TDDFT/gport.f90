module gport
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char, c_int16_t
  implicit none
  private
  public :: MAKEDIRQQ

  interface
     function c_mkdir(path, mode) bind(C, name="mkdir")
       use, intrinsic :: iso_c_binding
       implicit none
       character(kind=c_char), dimension(*), intent(in) :: path
       integer(c_int16_t), value :: mode
       integer(c_int) :: c_mkdir
     end function c_mkdir
  end interface

  contains

    function MAKEDIRQQ(dirname) result(res)
      implicit none
      character(len=*), intent(in)  :: dirname
      logical :: res
      character(len=len(dirname) + 1, kind=c_char) :: c_dirname
      integer(c_int) :: status
      integer(c_int16_t) :: mode

      c_dirname = trim(dirname) // c_null_char
      mode = int(o'772', c_int16_t)

      ! Call the C function mkdir
      status = c_mkdir(c_dirname, mode)

      ! Check if the directory was successfully created
      res = (status == 0)

    end function MAKEDIRQQ

end module
