module precision_mod

  implicit none

  !*** Single precision
  integer, parameter :: si = selected_int_kind(8)
  integer, parameter :: sp = selected_real_kind(4)

  !*** Double precision
  integer, parameter :: di = selected_int_kind(16)
  integer, parameter :: dp = selected_real_kind(8)

  !*** Custom precision
  integer, parameter :: cp = selected_real_kind(4)

  !*** Special precision
  integer, parameter :: hp = selected_real_kind(8)

end module precision_mod
