module global_parameters_mod

  use precision_mod
  use constants_mod

  implicit none

  !*** Folders and files names
  character(len=256) :: working_axisem_dir, input_field_name
  character(len=256) :: input_veloc_name(3), input_stress_name(6)
  character(len=256) :: output_veloc_name(3), output_stress_name(6)
  character(len=256) :: output_field_name, input_point_file
  character(len=3)   :: coup_tool

  !*** Coordinates of source and mesh
  real(kind=cp) :: lat_src,lon_src,lat_mesh,lon_mesh

  !*** Others
  integer(kind=si) :: nbproc, nel, ntime
  integer(kind=si) :: ibeg, iend
  integer(kind=si), allocatable, dimension(:)  :: nb_stored
  real(kind=cp), allocatable, dimension(:,:,:) :: scoor, zcoor
  
  !*** AxiSEM mesh specifications
  integer(kind=si), parameter :: NGNOD=8, NGLLX=5, NGLLY=5, NGLLS=5, NGLLZ=5
    
  !*** Field arrays
  real(kind=sp), allocatable, dimension(:,:,:) :: data_read
  real(kind=sp), allocatable, dimension(:,:)   :: data_rec, stress_rec, stress_to_write
  real(kind=sp), allocatable, dimension(:,:)   :: strain_rec

  !*** Work arrays 
  real(kind=sp), allocatable, dimension(:,:)   :: data_reduce, stress_reduce
        
  !*** Input box point 
  integer(kind=si) :: nbrec
  real(kind=cp), allocatable, dimension(:,:) :: reciever_cyl, reciever_geogr, reciever_sph
  real(kind=cp), allocatable, dimension(:)   :: reciever_interp_value, xi_rec, eta_rec

  integer(kind=si), allocatable :: rec2elm(:)
  
  !*** Rotation matrix with respect to the source
  real(kind=cp), dimension(3,3) :: rot_mat, trans_rot_mat

  !*** Rotation matrix with respect to the chosen coupling mesh
  real(kind=cp), dimension(3,3) :: rot_mat_mesh, trans_rot_mat_mesh
  real(kind=cp), dimension(3,3) :: mat, tmat

  !*** Info about simulation
  integer(kind=si) :: nsim
  character(len=100), allocatable, dimension(:)   :: simdir
  character(len=10),  allocatable, dimension(:,:) :: src_type

  !*** Post process
  real(kind=sp), allocatable, dimension(:) :: f1, f2, phi
  
  !*** Mpi 
  integer(kind=si) :: irecmin, irecmax

end module global_parameters_mod
