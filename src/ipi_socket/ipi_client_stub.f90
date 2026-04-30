module ipi_client_stub
  implicit none
contains
  subroutine ipi_socket_not_implemented()
    error stop "i-PI socket support is not implemented yet; standalone validation comes first."
  end subroutine ipi_socket_not_implemented
end module ipi_client_stub

