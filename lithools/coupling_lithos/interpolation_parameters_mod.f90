module interpolation_parameters_mod

  use precision_mod
  use constants_mod
  use mpi

  implicit none

  !*** Array for wavefields
  real(kind=cp), dimension(:,:,:), allocatable :: vel_inc, trac_inc
  real(kind=cp), dimension(:,:),   allocatable :: vel_inc2, stress_inc

  !*** Mapping arrays
  integer(kind=si), dimension(:),       allocatable :: abs_bnd_ielem
  integer(kind=si), dimension(:,:),     allocatable :: mapipt
  integer(kind=si), dimension(:,:,:),   allocatable :: abs_bnd_ijk
  integer(kind=si), dimension(:,:,:,:), allocatable :: loc2glob

  !*** Boundary informations
  real(kind=cp), dimension(:,:),   allocatable :: abs_bnd_jacobian2Dw
  real(kind=cp), dimension(:,:,:), allocatable :: abs_bnd_normal

  !*** Coordinates
  real(kind=cp), dimension(:), allocatable :: xcoord, ycoord, zcoord

  !*** Arrays with AxiSEM wavefields
  real(kind=cp), dimension(:,:), allocatable :: vxold,  vyold,  vzold
  real(kind=cp), dimension(:,:), allocatable :: sxxold, syyold, szzold
  real(kind=cp), dimension(:,:), allocatable :: syzold, sxzold, sxyold

  real(kind=cp), dimension(:,:), allocatable :: vxold1,  vyold1,  vzold1
  real(kind=cp), dimension(:,:), allocatable :: sxxold1, syyold1, szzold1
  real(kind=cp), dimension(:,:), allocatable :: syzold1, sxzold1, sxyold1

  real(kind=cp), dimension(:,:), allocatable :: vxold2,  vyold2,  vzold2
  real(kind=cp), dimension(:,:), allocatable :: sxxold2, syyold2, szzold2
  real(kind=cp), dimension(:,:), allocatable :: syzold2, sxzold2, sxyold2

  real(kind=cp), dimension(:), allocatable   :: vpow, stalta, taper
  real(kind=cp), dimension(:,:), allocatable :: data_tmp

  !*** Parameters
  integer(kind=si) :: ntold, ntnew, npts, num_bnd_faces, i, j, k, ipt, itold, itnew, ind, oldlen
  integer(kind=si) :: ngllsquare, iface, igll, ielem, iglob, iunit, isim, nsim, omp_proc, nanan
  real(kind=cp) :: vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, nx, ny, nz, alpha_mesh
  real(kind=cp) :: tx, ty, tz, x, y, z, feold, dtold, dtnew, lat_src, lon_src, lat_mesh, lon_mesh
  real(kind=cp) :: sumx,sumy,sumz
  character(len=250) :: input_field_name, output_field_name, input_point_file, fichier, dummy, rep
  character(len=3)   :: fwdtool
  character(len=250), dimension(:), allocatable :: working_axisem_dir
  character(len=250), dimension(3) :: output_veloc_name
  character(len=250), dimension(6) :: output_stress_name
  integer(kind=si), dimension(:), allocatable ::  ivx, ivy, ivz, isxx, isyy, iszz, isxy, isxz, isyz
  integer(kind=si) :: nbrec, ntime, istap, itime

  integer(kind=si) :: nptsa, nelx, nely, nelz, nelem, ngllx, nglly, ngllz, nlta, nsta, itbeg, itend
  real(kind=cp)    :: thres, alpha, tbeg, tend, fmax

  !** For filter
  real(kind=cp), dimension(:), allocatable :: Ker, conv, taptap
  integer(kind=si) :: ibeg, iend, nt
  integer(kind=si) :: n1
  real(kind=cp) :: dt, ttt2, ttt1
  real(kind=cp), dimension(:), allocatable :: tab_sinc 

  !*** From vadim
  integer(kind=si), dimension(MPI_STATUS_SIZE) :: statut
  integer(kind=si), parameter :: etq=100
  integer(kind=si) :: nb_rec_by_proc, nb_remain_proc, irecmin, irecmax, iproc, nrec_to_store, irec
  integer(kind=si), dimension(:), allocatable :: i_inf, i_sup, nb_received, shift, nb_received_sv, shift_sv


contains





end module interpolation_parameters_mod
