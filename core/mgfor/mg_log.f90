module mg_log
  use mg_types

  implicit none

contains

  subroutine logprint(mg, msg)
    type(MGinfos) :: mg
    character(len=*):: msg

    if(mg%param%debug>0)then
       print*, "[LOG] ", msg
    endif
  end subroutine logprint
end module mg_log
