module constants_mod

  use precision_mod
  
  !*** Constants for inititalization
  real(kind=cp), parameter :: mytiny = 1.e-9_cp
  real(kind=cp), parameter :: myhuge = 1.e30_cp

  real(kind=cp),    parameter :: myverysmall = 1.e-25_cp
  real(kind=dp),    parameter :: mysmall     = 0.000001_dp
  integer(kind=si), parameter :: myhugeint   = 100000000

  real(kind=cp), parameter :: zero = 0._cp
  real(kind=cp), parameter :: one  = 1._cp
  real(kind=cp), parameter :: half = 0.5_cp

  real(kind=cp), parameter :: onethird  = 1._cp/3._cp
  real(kind=cp), parameter :: twothird  = 2._cp/3._cp
  real(kind=cp), parameter :: foorthird = 4._cp/3._cp

  !*** Constants 
  real(kind=cp), parameter :: pi  = 3.141592653589793_cp
  real(kind=cp), parameter :: twopi = 2._cp * pi

  !*** Constants for conversions
  real(kind=cp), parameter :: deg2rad = 0.017453292519943
  real(kind=cp), parameter :: rad2deg = 57.295779513082323 

  !*** Check stability in modeling
  real(kind=cp), parameter :: stability_criterion = 1.e28_cp

end module
