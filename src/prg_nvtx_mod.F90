module prg_nvtx_mod
  use iso_c_binding
  implicit none


  integer(4),private :: col(7) = [ int(Z'0000ff00',4), int(Z'000000ff',4), &
  & int(Z'00ffff00',4), int(Z'00ff00ff',4), int(Z'0000ffff',4), &
  & int(Z'00ff0000',4), int(Z'00ffffff',4)]
  character(kind=c_char,len=256),private,target :: tempName

  type, bind(C):: nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1
     type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes

#if defined(USE_NVTX)
  interface gpmdRangePush
     ! push range with custom label and standard color
     subroutine gpmdRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR,len=*) :: name
     end subroutine gpmdRangePushA

     ! push range with custom label and custom color
     subroutine gpmdRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import:: nvtxEventAttributes
       type(nvtxEventAttributes):: event
     end subroutine gpmdRangePushEx
  end interface gpmdRangePush

  interface gpmdRangePop
     subroutine gpmdRangePop() bind(C, name='nvtxRangePop')
     end subroutine gpmdRangePop
  end interface gpmdRangePop

  public :: gpmdStartRange, gpmdEndRange

contains

  subroutine gpmdStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
    type(nvtxEventAttributes):: event

    tempName=trim(name)//c_null_char
!    write(*,*)"gpmdcov_nvtx"," Tag = "//name
    if ( .not. present(id)) then
       call gpmdRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call gpmdRangePushEx(event)
    end if
  end subroutine gpmdStartRange

  subroutine gpmdEndRange
    call gpmdRangePop
  end subroutine gpmdEndRange
#endif 
!END if defined USE_NVTX
  
end module prg_nvtx_mod
