!!
!! *PROGRAM* uff_merger
!!
!! *USAGE*: uff_merger MASTER.uff SLAVE.uff z1 z2 buff1 buff2 output.uff
!!
!! @version ecsim1.2.1
!!
!! *SRC_FILE*
!!
!! tools/product_tools/src/uff_merger.f90
!!
!! *LAST CHANGES*
!!
!! -Nov 30, 2011, Changed length of string1 and string2 to 256 to match proper length of scattering libs string and changed gass name writing
!! -Nov 26, 2011, Fixed errors in merged density/pressure profile calculations
!! -Nov 24, 2011, Replaced read_global_size_dist by interpolating version so that mergeing of uffs with different altutides grids would be done properly.  
!! -Nov 22, 2010, Added bounds check in merging of global scattering size dists section. Added missing deallocation statements (IMPORTANT FIX)
!! -Mar 10, 2008, D.D. Added exit(0) on exit
!! -Dec 18, 2007, R.V. Added 'header_' to UFF filename
!! -Nov 16, 2007, D.D. Added error checking to see if master==slave
!! -Nov 12, 2007, D.D. Converted from the old utility
!!
!! *DESCRIPTION* 
!!
!! This program merges two UFF files in the vertical coordinate.
!! This program reads a master UFF file and a slave UFF file. The information in the 
!! slave UFF file is then put on the same vertical grid as the master.
!! The master and slave files must have the same X and Y dimensions.
!!
!! Between z1 and z2 the cloud/aerosol, temperature and gas volume mixing
!! ratios of the slave file are used. Between z1 and 
!! z1-buff1 and  z2 and  z2+buff2
!! the temperature profile is obtained by linearly interpolating between
!! T_master(z1-buf) and T_slave(z1) and T_slave(z2) and
!! T_slave(z2+buf2) respectively. A similar procedure is used for the 
!! gas volume mixing ratios. The pressure profile is then constructed from the surface pressure 
!! and the resulting temperature and water vapor profile assuming
!! hydrostatic equilibrium conditions.
!!
!! The following arguments are expected:
!!
!! -MASTER          : Master UFF file.
!! -SLAVE           : Slave UFF file.
!! -z1              : Starting altitude to use slave info.
!! -z2              : Ending altitude to use slave info.
!! -buff1           : Buffer distance over which t,p,rho etc fields are 
!!                  : merged (z1-buff1 to z1). 
!! -buff2           : Buffer distance over which t,p,rho etc fields are 
!!                  : merged (z2 to z2+buff2). 
!! -OUTFILE          : name of output UFF file.
!!
!!
!!---------
!! OUTPUT |
!!---------
!! A single file OUTFILE.uff is generated.
!!
program uff_merger
  !
  Use write_messages
  Use scene_creator_types
  Use Data_types
  Use physical_parameters
  Use read_uff
  !
  Implicit None
  !
  !----------------------
  ! Character variables
  !----------------------
  !
  Character(len=256) :: master_uff
  Character(len=256) :: slave_uff
  Character(len=256) :: outfile_uff
  Character(len=256) :: work_str
  !
  Character(len=256)                   :: ecsim_home, scatt_lib,error_str
  !------------------------------------
  ! Altitudes to use data in secondary
  ! uff file from
  !------------------------------------
  !
  real              :: splice_z1,splice_z2
  !
  !-----------------------
  ! Global size-dist info
  !-----------------------
  !
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists_master
  Real,Dimension(:),Allocatable                      :: z_global_master
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists_slave
  Real,Dimension(:),Allocatable                      :: z_global_slave
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists_out
  Real,Dimension(:),Allocatable                      :: z_global_out
  !
  !-------------------
  ! Atmospheric Data
  !------------------
  !
  Type(atmos_point),Dimension(:),Allocatable :: data_point_master 
  Type(atmos_point),Dimension(:),Allocatable :: data_point_slave  
  Type(atmos_point),Dimension(:),Allocatable :: data_point_out    
  !
  !
  !---------------------------
  ! Gass+scatter information
  !---------------------------
  !
  Integer                                   :: nscatt_types_master,n_gasses_master
  Character(len=256), Dimension (:), Pointer :: scatt_list_names_master
  Character(len=10), Dimension (:), Pointer :: gasses_master
  !
  Integer                                   :: nscatt_types_slave,n_gasses_slave
  Character(len=256), Dimension (:), Pointer :: scatt_list_names_slave
  Character(len=10), Dimension (:), Pointer :: gasses_slave
  !
  Integer                                   :: nscatt_types_out,n_gasses_out
  Character(len=256), Dimension (:), Pointer :: scatt_list_names_out
  Character(len=10), Dimension (:), Pointer :: gasses_out
  !
  Type(scatt_prop_master),Dimension(:),Allocatable       :: scatt_master_info_master  ! structure containg scattering list info
  Type(scatt_prop_master),Dimension(:),Allocatable       :: scatt_master_info_slave  ! structure containg scattering list info

  !--------------------
  ! UFF header info
  !--------------------
  !
  integer   :: scene_nx_master,scene_ny_master,scene_nz_master
  integer   :: scene_nx_slave,scene_ny_slave,scene_nz_slave,scene_nz_slave_m
  integer   :: scene_nx_out,scene_ny_out,scene_nz_out
  Real      :: x1_master,x2_master,y1_master,y2_master,z1_master,z2_master,&
       &lat1_master,lat2_master,long1_master,long2_master
  Real      :: x1_slave,x2_slave,y1_slave,y2_slave,z1_slave,z2_slave,&
       &lat1_slave,lat2_slave,long1_slave,long2_slave
  Real      :: x1_out,x2_out,y1_out,y2_out,z1_out,z2_out,&
       &lat1_out,lat2_out,long1_out,long2_out
  !
  Integer   :: ascii_or_bi_master,ascii_or_bi_slave
  Integer   :: ascii_or_bi_out=2
  !
  !
  !--------------------
  ! Merging variables
  !--------------------
  !
  Real      :: buff1,buff2 !,deltaz,dis,po,h,Tavg,g
  Real      :: g,c1,slopet,R,X_wv,deltaz,h,dis
  !
  !------------------------
  ! Misc Working variables
  !------------------------
  !
  Integer                        :: status,i,j,k,l,n,ii,im
  Integer                        :: ix,iy,iz,isc,insr,io,jo
  real                           :: zero=0.0
  Integer                        :: one_ar(2)
  real                           :: ztemp
  !
  Integer                        :: iz1,iz1mb,iz2,iz2pb,izstart
  !
  real                           :: res_master,res_slave
!!$  Integer                        :: nz_out,nz_master,nz_slave
!!$  Integer                        :: nx_out,nx_master,nx_slave
!!$  Integer                        :: ny_out,ny_master,ny_slave
!!$  Integer                        :: ix1_master,ix2_master,iy1_master,iy2_master
!!$  Integer                        :: ix1_slave,ix2_slave,iy1_slave,iy2_slave
  !
  ! UFF header update
  Integer :: loc
  Character(len=256) :: char1, char2
  !
  !---------------------------------------------------
  ! Interfaces for external functions and subroutines
  !---------------------------------------------------
  !
  Interface
     !
     Integer Function iargc()
     End Function iargc
     !
     Subroutine getarg(n,arg)
       Integer,Intent(in)                  :: n
       Character(len=*),Intent(out)        :: arg
     End Subroutine getarg
     !   
  End Interface
  !
  !-----------------------------------------
  ! GET  environment variable ECSIM_HOME 
  !-----------------------------------------
  !
  call getenv('ECSIM_HOME', ecsim_home)
  call getenv('SCATT_LIB', scatt_lib)
  !
  uff_index=0
  !
  !-----------------------------------------
  ! Get the Arguments from the command line
  !-----------------------------------------
  !
  Call get_data(status)
  !
  If (status.Ne.0) Then 
     call write_error('Command line arguments not present or incomplete')
     call exit(1)
  Endif
  !
  !-----------------------------------
  ! Is the master the same as the
  ! slave ?
  !----------------------------------
  !
  if (trim(adjustl(master_uff))==trim(adjustl(slave_uff))) then
     call write_error('MASTER == SLAVE: Not premitted !')
     call exit(1)
  endif
  !
  !----------------------------------
  ! Add ECSIM_HOME to the input and 
  ! output files
  !----------------------------------
  !
  loc = index(master_uff,'/',.true.)
  char1 = master_uff(1:loc)
  char2 = master_uff(loc+1:)
  master_uff = trim(adjustl(char1))//'header_'//trim(adjustl(char2))
  master_uff = trim(adjustl(ecsim_home))//master_uff

  loc = index(slave_uff,'/',.true.)
  char1 = slave_uff(1:loc)
  char2 = slave_uff(loc+1:)
  slave_uff = trim(adjustl(char1))//'header_'//trim(adjustl(char2))
  slave_uff = trim(adjustl(ecsim_home))//slave_uff

  outfile_uff=trim(adjustl(ecsim_home))//outfile_uff
  !
  !
  !---------------------------------------------------
  ! Open and read the the UFF file header information
  !---------------------------------------------------
  !
  call write_info('Reading the UFF file: '//trim(adjustl(master_uff)))
  !
  Call read_uff_header_info(master_uff,nscatt_types_master,scatt_list_names_master,&
       & n_gasses_master,gasses_master,scene_nx_master,scene_ny_master,scene_nz_master,&
       & x1_master,y1_master,lat1_master,long1_master,z1_master,&
       & x2_master,y2_master,lat2_master,long2_master,z2_master,ascii_or_bi_master)
  !
  call write_info('Reading the UFF file: '//trim(adjustl(slave_uff)))
  !
  Nullify(scatt_list_names_slave)
  Nullify(gasses_slave)
  !
  Call read_uff_header_info(slave_uff,nscatt_types_slave,scatt_list_names_slave,&
       & n_gasses_slave,gasses_slave,scene_nx_slave,scene_ny_slave,scene_nz_slave,&
       & x1_slave,y1_slave,lat1_slave,long1_slave,z1_slave,&
       & x2_slave,y2_slave,lat2_slave,long2_slave,z2_slave,ascii_or_bi_slave)
  !
  if (splice_z2.gt.z2_slave) then
     splice_z2=z2_slave
  endif
  !
  !-----------------------------------------------------
  ! Compare domains and other things to see if we can
  ! proceed
  !-----------------------------------------------------
  !
  call write_info('=======================================================')
  write(work_str,'(a,1x,i)') 'nscatt_types_master= ',nscatt_types_master
  call write_info(trim(adjustl(work_str)))
  do i=1,n_gasses_master
     call write_info('gasses_master= '//trim(adjustl(gasses_master(i))))
  enddo
  write(work_str,'(a,1x,i)') 'scene_nx_master=',scene_nx_master
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a,1x,i)') 'scene_ny_master=',scene_ny_master
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a,1x,i)') 'scene_nz_master=',scene_nz_master
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a)') 'x1_master,y1_master,lat1_master,long1_master,z1_master='
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(f,1x,f,1x,f,1x,f,1x,f)') x1_master,y1_master,lat1_master,long1_master,z1_master
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a)') 'x2_master,y2_master,lat2_master,long2_master,z2_master='
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(f,1x,f,1x,f,1x,f,1x,f)') x2_master,y2_master,lat2_master,long2_master,z2_master
  call write_info(trim(adjustl(work_str)))
  !
  call write_info('=======================================================')
  write(work_str,'(a,1x,i)') 'nscatt_types_slave= ',nscatt_types_slave
  call write_info(trim(adjustl(work_str)))
  do i=1,n_gasses_slave
     call write_info('gasses_slave= '//trim(adjustl(gasses_slave(i))))
  enddo
  write(work_str,'(a,1x,i)') 'scene_nx_slave=',scene_nx_slave
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a,1x,i)') 'scene_ny_slave=',scene_ny_slave
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a,1x,i)') 'scene_nz_slave=',scene_nz_slave
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a)') 'x1_slave,y1_slave,lat1_slave,long1_slave,z1_slave='
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(f,1x,f,1x,f,1x,f,1x,f)') x1_slave,y1_slave,lat1_slave,long1_slave,z1_slave
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(a)') 'x2_slave,y2_slave,lat2_slave,long2_slave,z2_slave='
  call write_info(trim(adjustl(work_str)))
  write(work_str,'(f,1x,f,1x,f,1x,f,1x,f)') x2_slave,y2_slave,lat2_slave,long2_slave,z2_slave
  call write_info(trim(adjustl(work_str)))
  !
  if (nscatt_types_master.lt.nscatt_types_slave) then
     call write_error('nscatt_types_master.lt.nscatt_types_slave')
     call exit(1)
  endif
  !
  call write_info('Checking to see if the scattering types match')
  !
  do i=1,nscatt_types_slave
     !
     call write_info(trim(adjustl(scatt_list_names_master(i))))
     call write_info(trim(adjustl(scatt_list_names_slave(i))))
     !
     if (trim(adjustl(scatt_list_names_master(i))).ne.trim((adjustl(scatt_list_names_slave(i))))) then
        call write_error('scatt_list_names_master(i).ne.scatt_list_names_slave(i)')
        call exit(1)
     endif
  enddo
  !
  if (n_gasses_master.ne.n_gasses_master) then
     call write_error('n_gasses_master.ne.n_gasses_master')
     call exit(1)
  endif
  !
  do i=1,n_gasses_master
     if (gasses_master(i).ne.gasses_master(i)) then
        call write_error('gasses_master.ne.gasses_master')
        call exit(1)
     endif
  enddo
  !
!!$  if (x1_master.gt.x1_slave) then
!!$     call write_error('x1_master.gt.x1_slave')
!!$     call exit(1)
!!$  endif
!!$  !
!!$  if (y1_master.gt.y1_slave) then
!!$     call write_error('y1_master.gt.y1_slave')
!!$     call exit(1)
!!$  endif
!!$  !
!!$  if (x2_master.lt.x2_slave) then
!!$     call write_error('x2_master.lt.x2_slave')
!!$     call exit(1)
!!$  endif
!!$  !
!!$  if (y2_master.lt.y2_slave) then
!!$     call write_error('y2_master.lt.y2_slave')
!!$     call exit(1)
!!$  endif
!!$  !
!!$  if (scene_nx_master.gt.scene_nx_slave) then
!!$     call write_error('scene_nx_master.gt.scene_nx_slave')
!!$     call exit(1)
!!$  endif
  !
  if (scene_ny_master.gt.scene_ny_slave) then
     call write_error('scene_ny_master.gt.scene_ny_slave')
     call exit(1)
  endif
  !
  res_master=(x2_master-x1_master)/scene_nx_master
  res_slave=(x2_slave-x1_slave)/scene_nx_slave
  !
  if (res_master.ne.res_slave) then
     call write_warning('res_master.ne.res_slave (x)')
  endif
  !
  res_master=(y2_master-y1_master)/scene_ny_master
  res_slave=(y2_slave-y1_slave)/scene_ny_slave
  !
  if (res_master.ne.res_slave) then
     call write_warning('res_master.ne.res_slave (y)')
  endif
  !
  call write_info('======================================================')
  call write_info('       Things match enough to continue                ')
  call write_info('======================================================')
  !
  nscatt_types_out=nscatt_types_master
  scatt_list_names_out=scatt_list_names_master
  n_gasses_out=n_gasses_master
  gasses_out=gasses_master
  scene_nx_out=scene_nx_master
  scene_ny_out=scene_ny_master
  scene_nz_out=scene_nz_master
  !
  x1_out=x1_master
  x2_out=x2_master
  y1_out=y1_master
  y2_out=y2_master
  !
  long1_out=long1_master
  long2_out=long2_master
  lat1_out=lat1_master
  lat2_out=lat2_master
  !
  z1_out=z1_master
  z2_out=z2_master
  !
  !-----------------------------------
  ! Read the scattering master lists
  !-----------------------------------
  !
  Allocate(scatt_master_info_master(1:nscatt_types_master))
  Allocate(scatt_master_info_slave(1:nscatt_types_slave))
  !
  Do i=1,nscatt_types_master
     call write_info('Reading '//trim(scatt_list_names_master(i)))
     Call read_scatt_list_xml(5,scatt_list_names_master(i),scatt_master_info_master(i),status)
     if (status.ne.0) then  
        error_str='Error in reading scattering info'
        call exit(1)
     endif
  Enddo
  !
  Do i=1,nscatt_types_slave
     call write_info('Reading '//trim(scatt_list_names_slave(i)))
     Call read_scatt_list_xml(5,scatt_list_names_slave(i),scatt_master_info_slave(i),status)
     if (status.ne.0) then  
        error_str='Error in reading scattering info'
        call exit(1)
     endif
  Enddo
  !
  !
  !---------------------------------------------
  ! Read the global scattering property section 
  ! of the master UFF
  !---------------------------------------------
  !
  Allocate(global_size_dists_master(1:nscatt_types_master,1:scene_nz_master))
  Allocate(z_global_master(1:scene_nz_master))
  !
  call write_info('Reading global size dists from master uff file '//trim(master_uff))
  !
  Call read_global_size_dists(master_uff,ascii_or_bi_master,scene_nz_master,z_global_master,&
       & nscatt_types_master,scatt_master_info_master,global_size_dists_master)
  !
  !
  !
  ! D.D. Nov 24, 2011
  !
  scene_nz_slave_m=0
  do iz=1,scene_nz_out
     if (z_global_master(iz).le.z2_slave) then
        scene_nz_slave_m=scene_nz_slave_m+1
     else
        !        scene_nz_slave_m=scene_nz_slave_m-1
        exit
     endif
  enddo
  !
  !
  Allocate(global_size_dists_slave(1:nscatt_types_out,1:scene_nz_slave_m))
  Allocate(z_global_slave(1:scene_nz_slave_m))
  !
  call write_info('Reading global size dists from slave uff file '//trim(slave_uff))
  !
  !
  Call read_global_size_dists_intpol(slave_uff,ascii_or_bi_slave,scene_nz_slave,z_global_slave,&
       & scene_nz_slave_m,z_global_master(1:scene_nz_slave_m),nscatt_types_slave,scatt_master_info_slave,global_size_dists_slave)
  !
  ! D.D. Nov 24, 2011
  !
  !Subroutine read_global_size_dists_intpol(infilename,ascii_or_bi,nz,z,&
  !        nz_w,z_w,nscatt_types,scatt_master_info,global_size_dists)
  !
  Allocate(global_size_dists_out(1:nscatt_types_out,1:scene_nz_out))
  Allocate(z_global_out(1:scene_nz_out))
  z_global_out=z_global_master
  !
  do k=1,scene_nz_out
     !
     !write(*,*) k,scene_nz_out,z_global_out(k),splice_z1,splice_z2
     !
     if ((z_global_out(k).lt.splice_z1).or.(z_global_out(k).gt.splice_z2)) then
        do isc=1,nscatt_types_master
           Call Nullify_size_dist(global_size_dists_out(isc,k))
           if (Associated(global_size_dists_master(isc,k)%N_bin)) then
              Call allocate_size_dist(global_size_dists_out(isc,k),&
                   & scatt_master_info_master(isc)%num_sizes_want)
              global_size_dists_out(isc,k)%N_bin(:)=global_size_dists_master(isc,k)%N_bin(:)
              global_size_dists_out(isc,k)%n_sizes=scatt_master_info_master(isc)%num_sizes_want
           endif
        enddo
     else if (z_global_out(k).ge.splice_z1) then
        do isc=1,nscatt_types_master
           Call Nullify_size_dist(global_size_dists_out(isc,k))
           !
           ! D.D. Nov 22, 2011
           !
           if (k.le.scene_nz_slave_m) then
              if (Associated(global_size_dists_slave(isc,k)%N_bin)) then
                 Call allocate_size_dist(global_size_dists_out(isc,k),&
                      & scatt_master_info_master(isc)%num_sizes_want)
                 global_size_dists_out(isc,k)%N_bin(:)=global_size_dists_slave(isc,k)%N_bin(:)
                 global_size_dists_out(isc,k)%n_sizes=scatt_master_info_master(isc)%num_sizes_want
              endif
           endif
        enddo
     endif
  enddo
  !
  !-------------------------
  ! Allocate+nullify things
  !-------------------------
  !
  allocate(data_point_master(1:scene_nz_master))
  do iz=1,scene_nz_master
     Call Nullify_atmos_point(data_point_master(iz))
  enddo
  allocate(data_point_out(1:scene_nz_out))
  do iz=1,scene_nz_out
     Call Nullify_atmos_point(data_point_out(iz))
  enddo
  !
  !
  !--------------
  ! Main loop !!
  !--------------
  !
  iz1mb=1
  !
  do i=1,scene_nz_out
     !
     if (z_global_out(i).lt.splice_z2) then
        im=i
     endif
     !
     if (z_global_out(i).lt.splice_z1-buff1) then
        iz1mb=i
     endif
     !
     if (z_global_out(i).lt.splice_z1) then
        iz1=i
     endif
     !
     if (z_global_out(i).lt.splice_z2+buff2) then
        iz2pb=i
     endif
     !
     if (z_global_out(i).lt.splice_z2) then
        iz2=i
     endif
     !
  enddo
  !
  ! D.D. Nov 24, 2011
  !
  im=scene_nz_slave_m
  !
  allocate(data_point_slave(1:im))
  do iz=1,im
     Call Nullify_atmos_point(data_point_slave(iz))
  enddo
  !
  do ix=1,scene_nx_out
     !
     call write_prog(ix,scene_nx_out)
     !
     do iy=1,scene_ny_out
        !
        !
        !-------------------------
        ! Read the master column
        !-------------------------
        !
        Call Read_uff_column(master_uff,ascii_or_bi_master,3,ix,iy,scene_nx_master,scene_ny_master,&
             & scene_nz_master,data_point_master,&
             & n_gasses_master,nscatt_types_master,scatt_master_info_master)
        !
        !if ((data_point_master(1)%x.ge.x1_slave).and.(data_point_master(1)%x.le.x2_slave).and.&
        !     &(data_point_master(1)%y.ge.y1_slave).and.(data_point_master(1)%y.le.y2_slave)) then
        !
        !----------------------------
        ! Read the slave info at the
        ! master resolution
        !----------------------------
        !
        Call Read_uff_column_intpol(slave_uff,ascii_or_bi_slave,4,ix,iy,&
             & scene_nx_slave,scene_ny_slave,&
             & scene_nz_slave,data_point_slave,&
             & im,z_global_out(1:im),n_gasses_out,nscatt_types_slave,scatt_master_info_slave)        
        !
        !
        do iz=1,im
           !
           if ((data_point_slave(iz)%z.ge.splice_z1).and.&
                &(data_point_slave(iz)%z.le.splice_z2)) then
              call copy_atmos_point_diff_nscatt(data_point_slave(iz),data_point_out(iz),&
                   & nscatt_types_out)
           else
              call copy_atmos_point(data_point_master(iz),data_point_out(iz))
           endif
           !
           data_point_out(iz)%x=data_point_master(1)%x
           data_point_out(iz)%y=data_point_master(1)%y
           !
        enddo
        !
        do iz=im+1,scene_nz_out
           call copy_atmos_point(data_point_master(iz),data_point_out(iz))
        enddo
        !
        !else
        !   !
        !   do iz=1,scene_nz_out
        !      call copy_atmos_point(data_point_master(iz),data_point_out(iz))
        !   enddo
        !   !
        !endif
        !
        !--------------------------------------------------
        ! In the buffer regions use linear interpolation
        ! to determine t,RH, and gas volumne mixing ratios
        !---------------------------------------------------
        !
        !  z1-buff to z1
        !
        !write(*,*) data_point_out(iz1)%z
        !write(*,*) data_point_out(iz1mb)%z
        !write(*,*) data_point_out(iz)%z-data_point_out(iz1mb)%z
        !
        if ((iz1mb.gt.1).and.(iz1mb.lt.iz1)) then
           !
           deltaz=data_point_out(iz1)%z-data_point_out(iz1mb)%z
           !
           do iz=iz1mb,iz1
              !
              dis=data_point_out(iz)%z-data_point_out(iz1mb)%z
              !
              data_point_out(iz)%t=(data_point_out(iz1)%t-data_point_out(iz1mb)%t)/deltaz*dis+&
                   & data_point_out(iz1mb)%t
              !
              data_point_out(iz)%RH=(data_point_out(iz1)%RH-data_point_out(iz1mb)%RH)/deltaz*dis+&
                   & data_point_out(iz1mb)%RH
              !
              do i=1,data_point_master(iz)%n_gasses
                 data_point_out(iz)%X_vol(i)=(data_point_out(iz1)%X_vol(i)-&
                      & data_point_out(iz1mb)%X_vol(i))/deltaz*dis+&
                      & data_point_out(iz1mb)%X_vol(i)
              enddo
              !
           enddo
        endif
        !
        ! z2 to z2+buf
        !
        deltaz=data_point_out(iz2pb)%z-data_point_out(iz2)%z
        do iz=iz2,iz2pb
           !
           ! D.D. Nov 28, 2011
           !
           dis=data_point_out(iz)%z-data_point_out(iz2)%z
           !
           data_point_out(iz)%t=(data_point_out(iz2pb)%t-data_point_out(iz2)%t)/deltaz*dis+&
                & data_point_out(iz2)%t
           !
           data_point_out(iz)%RH=(data_point_out(iz2pb)%RH-data_point_out(iz2pb)%RH)/deltaz*dis+&
                & data_point_out(iz2)%RH
           !
           do i=1,data_point_master(iz)%n_gasses
              data_point_out(iz)%X_vol(i)=(data_point_out(iz2pb)%X_vol(i)-&
                   & data_point_out(iz2)%X_vol(i))/deltaz*dis+&
                   & data_point_out(iz2)%X_vol(i)
           enddo
           !
        enddo
        !
        !-----------------------------------------------------------------------
        ! Now calculate the entrie pressure and density starting from the 
        ! surface pressure (or top of slave) and the merged temperature profile
        !------------------------------------------------------------------------
        !
        if (iz1mb.gt.1) then
           izstart=2
        else
           izstart=iz2+1
        endif
        !
        !
        do iz=izstart,scene_nz_out
           !
           h=data_point_out(iz)%z*1000.0
           !
           ! gravity profile (see CIRA 1965, page 5)
           g= 9.793244 - h*3.086597E-6 &
                &  + (h**2)*7.259E-13
           !
           X_wv=(data_point_out(iz-1)%X_vol(iwater)+data_point_out(iz)%X_vol(iwater))/2.0
           !
           R=Rg*(1.0+0.61*X_wv)      ! J/mol/K
           slopet=(data_point_out(iz)%t-data_point_out(iz-1)%t)/(data_point_out(iz)%z-data_point_out(iz-1)%z) ! in K/km
           !                       
           c1=-(mw_air*g)/(R*slopet)   ! mw_air in g/mol
           c1=c1-1
           !
           data_point_out(iz)%rho=data_point_out(iz-1)%rho*(data_point_out(iz)%t/data_point_out(iz-1)%t)**c1 ! in #/cm^3
           !
           data_point_out(iz)%p=(data_point_out(iz)%rho*1.0e+6/av_num)*Rg*data_point_out(iz)%T/100.0    ! im mb
           !
        enddo
        !
        !----------------------------------
        ! Write info out to the UFF file !
        !----------------------------------
        !
        io=ix
        jo=iy
        !
        call write_bi_di(outfile_uff)
        !
        !
        do iz=1,scene_nz_master
           Call deallocate_atmos_point(data_point_master(iz))
        enddo
        !
        do iz=1,scene_nz_slave
           Call deallocate_atmos_point(data_point_slave(iz))
        enddo
        !
        do iz=1,scene_nz_out
           Call deallocate_atmos_point(data_point_out(iz))
        enddo
        !
        !
        !if (ascii_or_bi_out.eq.0) then 
        !   call write_bi(outfile_uff,10)
        !else
        !   call write_ascii(outfile_uff,10)
        !endif
        !
     enddo
  enddo
  !
  Close(30)
  Do isc=1,nscatt_types_out
     Close(30+isc+10)
  Enddo
  !           
  Call exit(0)
  !
Contains
  !
  !
  Subroutine write_bi_di(filename)
    !
    Character(len=*),intent(in)          :: filename
    !
    Character(len=120),dimension(:),allocatable,save  :: size_dist_uff_files
    Character(len=120),save                           :: header_uff_file
    Character(len=120),save                           :: main_uff_file
    !
    integer                                           :: isc,irm,worki,id,id2,n
    integer,save                                      :: icheck
    real                                              :: workr
    Character(len=4)                                  :: str
    integer,dimension(:),allocatable,save             :: size_record_index
    Integer,dimension(:),allocatable,save             :: index_write 
    Character(len=256)                                :: string1
    Character(len=256)                                :: string2
    Logical                                           :: init
    !
    !------------------------------
    ! Write info out the UFF file !
    !------------------------------
    !
    If ((io.Eq.1).And.(jo.Eq.1)) Then
       !
       !
       ! Build the names
       !
       main_uff_file=filename
       !
       call get_path(main_uff_file,string1)
       !write(*,*) string1
       call get_filename(main_uff_file,string2)
       !write(*,*) string2
       !
       header_uff_file=Trim(string1)//'header_'//Trim(string2)
       !
       !write(*,*) header_uff_file
       !pause
       !
       allocate(size_dist_uff_files(1:nscatt_types_out))
       !
       do isc=1,nscatt_types_out
          !
          write(str,'(i4)') isc
          size_dist_uff_files(isc)=Trim(string1)//'size_dists_'//trim(trim(adjustl(str)))//'_'//Trim(string2)
          !
          !
       enddo
       !
       !----------------------------
       ! Open and write the header
       ! file
       !----------------------------
       !
       call write_info('Writing header file '//header_uff_file)
       open(unit=3,file=header_uff_file,status='unknown',iostat=status,form='formatted')
       !
       !-----------------------
       ! Register file names
       !-----------------------
       !
       uff_index=uff_index+1
       nx_uff(uff_index)=scene_nx_out
       ny_uff(uff_index)=scene_ny_out
       nz_uff(uff_index)=scene_nz_out
       !
       bi_dir_uff_file_header(uff_index)=header_uff_file
       bi_dir_uff_file_main(uff_index)=main_uff_file
       do isc=1,nscatt_types_out
          bi_dir_uff_size_dists_files(uff_index,isc)=size_dist_uff_files(isc)
       enddo
       !
       !----------------------------
       ! Write out the names of the 
       ! Scattering master lists 
       ! this scene requires
       !----------------------------
       !
       Write(3,*) 'ddddddd----------------------------------'
       Write(3,*) nscatt_types_out
       Write(3,*) '-----------------------------------------'
       Do i=1,nscatt_types_out
          call get_filename(scatt_list_names_master(i),string2)
          Write(3,'(a)') string2   ! Master list file for each scatterer type
       Enddo
       Write(3,*) '------------------------------------------'
       Write(3,*) n_gasses_out
       Write(3,*) '------------------------------------------'
       Do i=1,n_gasses_out
          Write(3,*) trim(adjustl(gasses_master(i)))
       Enddo
       Write(3,*) '-------------------------------------------------'
       Write(3,*) '           nx,ny,nz (scene dimensions)'
       Write(3,*) scene_nx_out,scene_ny_out,scene_nz_out
       Write(3,*) '--------------x1,y1,lat1,long1,z1----------------'
       Write(3,'(5f9.4)') x1_out,y1_out,lat1_out,long1_out,z1_out
       Write(3,*) '--------------x2,y2,lat2,long2,z2----------------'
       Write(3,'(5f9.4)') x2_out,y2_out,lat2_out,long2_out,z2_out
       Write(3,*) '-------------------------------------------------'
       Write(3,*) 'Names of main and size dist direct access files  '
       Write(3,*) '-------------------------------------------------'
       !
       call get_filename(main_uff_file,string2)
       Write(3,'(a)') string2
       !
       Do isc=1,nscatt_types_out
          call get_filename(size_dist_uff_files(isc),string2)
          Write(3,'(a)') string2
       Enddo
       !
       Write(3,*) '-------------------------------------------------'
       Write(3,*) 'Backgrnd scattering info (function of z only)'
       Write(3,*) '-------------------------------------------------'
       !
       !--------------------
       ! Write the global
       ! size dist section
       !--------------------
       !
       Do iz=1,scene_nz_out
          write(3,*) data_point_out(iz)%z 
          Do isc=1,nscatt_types_out
             if (associated(global_size_dists_out(isc,iz)%N_bin)) then
                n=global_size_dists_out(isc,iz)%n_sizes
                Write(3,'(i6,1x,i6)') isc,n
                do ii=1,n
                   write(3,200,advance='NO') global_size_dists_out(isc,iz)%N_bin(ii)/1.0e+6
                enddo
                write(3,*) 
             else
                n=0
                Write(3,'(i6,1x,i6)') isc,n
             Endif
          enddo
       enddo
       !
200    Format(1x,e12.6)    
       !
       Close(3)
       !
       !------------------------
       ! Open the main uff file
       !------------------------
       !
       id=30
       worki=(9+n_gasses_out)*8+(nscatt_types_out+1)*4
       open(unit=id,file=main_uff_file,status='replace',iostat=status,form='unformatted',&
            & access='direct',recl=worki)
       !
       !----------------------------
       ! Open the size distribution
       ! files
       !----------------------------
       !
       id2=1
       !
       do isc=1,nscatt_types_out
          worki=(scatt_master_info_master(isc)%num_sizes_want)*8
          open(unit=id+id2*10+isc,file=size_dist_uff_files(isc),status='replace',iostat=status,&
               & form='unformatted',&
               & access='direct',recl=worki)
       enddo
       !
       init=.true.
       !
       call write_uff_direct(id,id2,io,jo,0,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
            & data_point_out(1),init)
       !
       init=.false.
       do iz=1,scene_nz_out
          call write_uff_direct(id,id2,io,jo,iz,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
               & data_point_out(iz),init)    
       enddo
       !
    Else
       !
       id=30
       id2=1
       init=.false.
       call write_uff_direct(id,id2,io,jo,0,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
            & data_point_out(1),init)   
       !
       do iz=1,scene_nz_out
          !
          call write_uff_direct(id,id2,io,jo,iz,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
               & data_point_out(iz),init)    
       enddo
       !
    Endif
    !
    return
    !
  End Subroutine write_bi_di
  !
  Subroutine write_ascii(filename,id)
    !
    Integer,intent(in)                   :: id
    Character(len=*),intent(in)          :: filename
    !
    !----------------------------------
    ! Write info out the the UFF file !
    !----------------------------------
    !
    If ((ix.Eq.1).And.(iy.Eq.1)) Then
       !
       !
       Open(unit=id,file=filename,status='unknown',iostat=status,form='formatted')
       !
       !----------------------------
       ! Write out the names of the 
       ! Scattering master lists 
       ! this scene requires
       !----------------------------
       !
       zero=0.0
       !
       Write(id,*) '-----------------------------------------'
       Write(id,*) nscatt_types_out,'  ,Number of scatter types'  
       Write(id,*) '-----------------------------------------'
       Do i=1,nscatt_types_out
          Write(id,'(a80)') scatt_list_names_master(i)   ! Master list file for each scatterer type
       Enddo
       Write(id,*) '------------------------------------------'
       Write(id,*) n_gasses_out,'  ,Number of atmospheric gasses'
       Write(id,*) '------------------------------------------'
       Do i=1,n_gasses_out
          Write(id,*) gasses_master(i) ! name of gass
       Enddo
       Write(id,*) '-------------------------------------------------'
       Write(id,*) '           nx,ny,nz (scene dimensions)'
       Write(id,*) scene_nx_out,scene_ny_out,scene_nz_out
       Write(id,*) '--------------x1,y1,lat1,long1,z1----------------'
       Write(id,'(5f9.4)') x1_out,y1_out,zero,zero,data_point_out(1)%z 
       Write(id,*) '--------------x2,y2,lat2,long2,z2----------------'
       Write(id,'(5f9.4)') x2_out,y2_out,zero,zero,data_point_out(scene_nz_out)%z
       Write(id,*) '-------------------------------------------------'
       Write(id,*) 'Backgrnd scattering info (function of z only)'
       Write(id,*) '-------------------------------------------------'
       !
       Do iz=1,scene_nz_out
          write(id,*) data_point_out(iz)%z 
          Do isc=1,nscatt_types_out
             !
             If (Associated(global_size_dists_out(isc,iz)%N_bin)) Then
                !
                j=size(global_size_dists_out(isc,iz)%N_bin)
                !
                Write(id,'(i6,1x,i6)') isc,j
                do ii=1,j
                   write(id,200,advance='NO') global_size_dists_out(isc,iz)%N_bin(ii)*1.0e-6 ! put back in /cm^3
                enddo
                write(id,*)
             Else
                j=0
                Write(id,'(i6,1x,i6)') isc,j
             Endif
          Enddo
       Enddo
    Endif
    !
    Do iz=1,scene_nz_out
       If (iz==1) Then 
          Write(id,*) '-------------------------------------------------'
          Write(id,'(a20,1x,i6,1x,i6)') 'Info for column     ',ix,iy
          Write(id,*) '-------------------------------------------------'
          !
          Write(id,'(a20,1x,i6,1x,i6)') 'Surface type and subtype',data_point_out(1)%Surface,int(zero)
          Write(id,*) '-------------------------------------------------'
          !
       endif
       !
       Write(id,100) data_point_out(iz)%x,data_point_out(iz)%y,&
            &data_point_out(iz)%lat,data_point_out(iz)%long,&
            &data_point_out(iz)%z,&
            &data_point_out(iz)%t,&
            &data_point_out(iz)%p,&
            &data_point_out(iz)%rho,&
            &data_point_out(iz)%n_gasses
       Do i=1,n_gasses_out
          Write(id,200,advance='NO') data_point_out(iz)%X_vol(i)
       Enddo
       i=n_gasses_out
       Write(id,200) data_point_out(iz)%X_vol(i)
       Write(id,300) data_point_out(iz)%x_vel,&
            & data_point_out(iz)%y_vel,&
            & data_point_out(iz)%z_vel,&
            & data_point_out(iz)%turb_dv
       Write(id,'(i6)') nscatt_types_out
       Do i=1,nscatt_types_out
          If (Associated(data_point_out(iz)%size_dists(i)%N_bin)) Then
             j=data_point_out(iz)%size_dists(i)%n_sizes
             Write(id,'(i6,1x,i6)') i,j
             do ii=1,j
                write(id,200,advance='NO') data_point_out(iz)%size_dists(i)%N_bin(ii)*1.0e-6
             enddo
             write(id,*)
          Else
             j=0
             Write(id,'(i6,1x,i6)') i,j
          Endif
       Enddo
    Enddo
100 Format(f7.3,1x,f7.3,1x,f8.4,1x,f8.4,1x,f7.3,1x,f7.3,1x,e10.4,1x,e12.6,1x,i3)
200 Format(1x,e12.6)    
300 Format(f7.3,1x,f7.3,1x,f7.3,1x,f7.3)
    !
  End Subroutine write_ascii
  !
  !
  Subroutine write_bi(filename,id)
    !
    Integer,intent(in)                   :: id
    Character(len=*),intent(in)          :: filename
    !
    !----------------------------------
    ! Write info out the the UFF file !
    !----------------------------------
    !
    If ((ix.Eq.1).And.(iy.Eq.1)) Then
       !
       !
       Open(unit=id,file=filename,status='unknown',iostat=status,form='unformatted',access='sequential')
       !
       !----------------------------
       ! Write out the names of the 
       ! Scattering master lists 
       ! this scene requires
       !----------------------------
       !
       Write(id) 'bbbbbbb----------------------------------'
       Write(id) nscatt_types_out
       Write(id) '-----------------------------------------'
       Do i=1,nscatt_types_out
          Write(id) scatt_list_names_master(i)   ! Master list file for each scatterer type
       Enddo
       Write(id) '------------------------------------------'
       Write(id) n_gasses_out
       Write(id) '------------------------------------------'
       Do i=1,n_gasses_out
          Write(id) gasses_master(i)
       Enddo
       Write(id) '-------------------------------------------------'
       Write(id) '           nx,ny,nz (scene dimensions)'
       Write(id) scene_nx_out,scene_ny_out,scene_nz_out 
       Write(id) '--------------x1,y1,lat1,long1,z1----------------'
       Write(id) x1_out,y1_out,zero,zero,data_point_out(1)%z
       Write(id) '--------------x2,y2,lat2,long2,z2----------------'
       Write(id) x2_out,y2_out,zero,zero,data_point_out(scene_nz_out)%z
       Write(id) '-------------------------------------------------'
       Write(id) 'Backgrnd scattering info (function of z only)'
       Write(id) '-------------------------------------------------'
       !
       Do iz=1,scene_nz_out
          write(id) data_point_out(iz)%z 
          Do isc=1,nscatt_types_out
             !
             If (Associated(global_size_dists_out(isc,iz)%N_bin)) Then
                !
                j=size(global_size_dists_out(isc,iz)%N_bin)
                !
                Write(id) isc,j
                write(id) global_size_dists_out(isc,iz)%N_bin(1:j)*1.0e-6
                !
             Else
                j=0
                Write(id) isc,j
             Endif
          Enddo
       Enddo
    Endif
    !
    Do iz=1,scene_nz_out
       If (iz==1) Then 
          Write(id) 'Info for column----------------------------------'
          Write(id) ix,iy
          Write(id) '-------------------------------------------------'
          !
          write(id) data_point_out(1)%Surface,int(zero)
          !
          Write(id) '-------------------------------------------------'
          !
       endif
       !
       Write(id) data_point_out(iz)%x,data_point_out(iz)%y,&
            &data_point_out(iz)%lat,data_point_out(iz)%long,&
            &data_point_out(iz)%z,&
            &data_point_out(iz)%t,&
            &data_point_out(iz)%p,&
            &data_point_out(iz)%rho,&
            &data_point_out(iz)%n_gasses
       Write(id) data_point_out(iz)%X_vol
       i=n_gasses_out
       Write(id) data_point_out(iz)%x_vel,&
            & data_point_out(iz)%y_vel,&
            & data_point_out(iz)%z_vel,&
            & data_point_out(iz)%turb_dv
       Write(id) nscatt_types_out 
       Do i=1,nscatt_types_master
          If (Associated(data_point_out(iz)%size_dists(i)%N_bin)) Then
             j=data_point_out(iz)%size_dists(i)%n_sizes
             Write(id) i,j
             if (i==3) then
                data_point_out(iz)%size_dists(i)%N_bin(20:j)=0.0
             endif
             write(id) data_point_out(iz)%size_dists(i)%N_bin(1:j)*1.0e-6
          Else
             j=0
             Write(id) i,j
          Endif
       Enddo
    Enddo
    !
  End Subroutine write_bi
  !
  !
  Subroutine get_data(status)
    !
    !-------------------------------------
    ! READS INPUT FROM THE COMMAND LINE
    !-------------------------------------
    ! STATUS=0 IF SUCESSFUL, 1 OTHERWISE
    !-------------------------------------
    !
    !-----------------
    ! Passed variables
    !-----------------
    !
    Integer    :: status
    !
    !-----------------
    ! Local variables
    !-----------------
    !
    Integer                          :: i,nargs
    !
    Character(len=130)                :: arg_str
    Character(len=100)               :: error_str
    Character(len=100)               :: error_str2
    !
    error_str='Usage: uff_merger header_master.uff header_slave.uff z1 z2 buff1 buff2 outfile.uff'
    !
    status=0
    !
    nargs=iargc()  
    If (nargs.Ne.7) Then
       status=1
       error_str2='Wrong number of arguments'
    Endif
    !
!!$    If (status == 0) Then
!!$       Call getarg(1,arg_str)
!!$       error_str2='Error in data filename specification'
!!$       Read(arg_str,'(a130)',iostat=status,err=100) infilename
!!$    Endif
    !
     If (status == 0) Then
        Call getarg(1,arg_str)
        error_str2='Attempting to read master UFF file'
        Read(arg_str,'(a)',iostat=status,err=100) master_uff
     Endif
    !
     If (status == 0) Then
        Call getarg(2,arg_str)
        error_str2='Attempting to read slave UFF file'
        Read(arg_str,'(a)',iostat=status,err=100) slave_uff
     Endif
    !
     If (status == 0) Then
        Call getarg(3,arg_str)
        error_str2='Attempting to read z1'
        Read(arg_str,*,iostat=status,err=100) splice_z1
     Endif
    !
    !
     If (status == 0) Then
        Call getarg(4,arg_str)
        error_str2='Attempting to read z2'
        Read(arg_str,*,iostat=status,err=100) splice_z2
     Endif
    !
     If (status == 0) Then
        Call getarg(5,arg_str)
        error_str2='Attempting to read buff1'
        Read(arg_str,*,iostat=status,err=100) buff1
     Endif
     !
     If (status == 0) Then
        Call getarg(6,arg_str)
        error_str2='Attempting to read buff2'
        Read(arg_str,*,iostat=status,err=100) buff2
     Endif
     !
     If (status == 0) Then
        Call getarg(7,arg_str)
        error_str2='Attempting to output UFF file'
        Read(arg_str,'(a)',iostat=status,err=100) outfile_uff
     Endif
     !
     !
100 If (status.Ne.0) Then 
       call Write_error(error_str)
       call Write_error(error_str2)
       call exit(1)
    Endif
    !
  End Subroutine get_data
  !
  !
End program uff_merger
