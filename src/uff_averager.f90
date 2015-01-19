!!
!! *PROGRAM* uff_averager
!!
!! *USAGE*: uff_averager input_header.uff output.uff hor_res[km] v_res[km] 
!!
!! @version ecsim0.3
!!
!! *SRC_FILE*
!!
!! tools/product_tools/src/uff_averager.f90
!!
!! *LAST CHANGES*
!!
!! - Nov 30, 2011, D.D. slight changes in the string length for gasses and scatt libs being written out
!! - Dec 18, 2007, RV: Added 'header_' to UFF filename
!! - Nov 09, 2007, D.D. Converted from the old utility
!! 
!! *DESCRIPTION*
!!
!! This program reads an UFF file and using simple averaging it outputs
!! a UFF file of lesser resolution. Size distributions are averaged 
!! on a bin-by-bin basis over the desired domain.
!!
!! The following arguments are expected
!!
!! input_data_file          : Orginal uff file (header file)
!!
!! output_data_file         : Output uff file.
!!
!! hor_res                  : Desired output horizontal resolution in km.
!!
!! v_res                    : Desired output vertical resolution in km.
!!
!!
program uff_averager
  !
  !
  Use write_messages
  Use scene_creator_types
  Use Data_types
  Use read_uff
  !
  Implicit None
  !
  !----------------------
  ! Character variables
  !----------------------
  !
  Character(len=256) :: uff_in
  Character(len=256) :: uff_out
  !
  !-------------
  ! Resolution
  !-------------
  !
  real              :: h_res,v_res
  Integer           :: n_hor,n_ver
  Integer           :: nh_to_avg,max_n_h
  !
  !
  !-----------------------
  ! Global size-dist info
  !-----------------------
  !
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists_master
  Real,Dimension(:),Allocatable                      :: z_global_master
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists_out
  Real,Dimension(:),Allocatable                      :: z_global_out
  !
  !-------------------
  ! Atmospheric Data
  !------------------
  !
  Type(atmos_point),Dimension(:),Allocatable :: data_point_master 
  Type(atmos_point),Dimension(:,:,:),Allocatable :: data_point_out
  !
  !
  !---------------------------
  ! Gass+scatter information
  !---------------------------
  !
  Integer                                   :: nscatt_types_master,n_gasses_master
  Integer                                   :: nscatt_types_out,n_gasses_out
  Character(len=256), Dimension (:), Pointer :: scatt_list_names_master
  Character(len=10), Dimension (:), Pointer :: gasses_master
  Type(scatt_prop_master),Dimension(:),Allocatable       :: scatt_master_info_master
  !
  !--------------------
  ! UFF header info
  !--------------------
  !
  integer   :: scene_nx_master,scene_ny_master,scene_nz_master
  integer   :: scene_nx_out,scene_ny_out,scene_nz_out
  Real      :: x1_master,x2_master,y1_master,y2_master,z1_master,z2_master,&
       &lat1_master,lat2_master,long1_master,long2_master
  Real      :: x1_out,x2_out,y1_out,y2_out,z1_out,z2_out,&
       &lat1_out,lat2_out,long1_out,long2_out
  !
  Integer   :: ascii_or_bi_master,ascii_or_bi_out
  !
  !------------------------
  ! Misc Working variables
  !------------------------
  !
  Integer,dimension(:,:,:),allocatable         :: n_hits                        
  Integer,dimension(:,:,:,:),allocatable       :: n_hits_st                        
  Integer,dimension(:,:),allocatable           :: n_hits_global             
  Integer                        :: status,i,j,k,l,n,ii,im,isc,io,jo,ko,ig,k1,kl
  Integer                        :: maxio,maxjo
  Integer                        :: ix1,ix2,iy1,iy2
  Integer                        :: iix1,iix2,iiy1,iiy2
  Real                           :: z1,z2
  Integer                        :: ix,iy,iz
  real                           :: zero=0.0
  Real                           :: h_res_master,v_res_master
  Real,dimension(:),allocatable  :: temp_dist,work_gasses
  Real                           :: work
  !
  Character(len=256)                   :: ecsim_home, scatt_lib,error_str
  !
  !
  ! UFF header update
  Integer :: loc
  Character(len=256) :: char1, char2
  !
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
       Character(len=*),Intent(out)           :: arg
     End Subroutine getarg
     !   
  End Interface
  !
  uff_index=0
  !
  !-----------------------------------------
  ! GET  environment variable ECSIM_HOME 
  !-----------------------------------------
  !
  call getenv('ECSIM_HOME', ecsim_home)
  call getenv('SCATT_LIB', scatt_lib)
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
  !----------------------------------
  ! Add ECSIM_HOME to the input and 
  ! output files
  !----------------------------------
  !
  loc = index(uff_in,'/',.true.)
  char1 = uff_in(1:loc)
  char2 = uff_in(loc+1:)
  uff_in = trim(adjustl(char1))//'header_'//trim(adjustl(char2))
  uff_in = trim(adjustl(ecsim_home))//uff_in
  uff_out= trim(adjustl(ecsim_home))//uff_out
  !
  !
  !---------------------------------------------------
  ! Open and read the the UFF file header information
  !---------------------------------------------------
  !
  call write_info('Reading the UFF file:'//trim(adjustl(uff_in)))
  !
  Call read_uff_header_info(uff_in,nscatt_types_master,scatt_list_names_master,&
       & n_gasses_master,gasses_master,scene_nx_master,scene_ny_master,scene_nz_master,&
       & x1_master,y1_master,lat1_master,long1_master,z1_master,&
       & x2_master,y2_master,lat2_master,long2_master,z2_master,ascii_or_bi_master)
  !
  nscatt_types_out=nscatt_types_master
  n_gasses_out=n_gasses_master
  !
  !-----------------------------------
  ! Read the scattering master lists
  !-----------------------------------
  !
  Allocate(scatt_master_info_master(1:nscatt_types_master))
  Do i=1,nscatt_types_master
     write(error_str,*) i,':',trim(adjustl(scatt_list_names_master(i)))
     call write_info(error_str)
     Call read_scatt_list_xml(5,scatt_list_names_master(i),scatt_master_info_master(i),status)
     if (status.ne.0) then  
        error_str='Error in reading scattering info'
        call write_info(error_str)
        call exit(1)
     endif
  Enddo
  !
  !
  !---------------------------------------------
  ! Read the global scattering property section 
  ! of the UFF
  !---------------------------------------------
  !
  Allocate(global_size_dists_master(1:nscatt_types_master,1:scene_nz_master))
  Allocate(z_global_master(1:scene_nz_master))
  !
  call write_info('Reading global size dists from master uff file '//trim(adjustl(uff_in)))
  !
  Call read_global_size_dists(uff_in,ascii_or_bi_master,scene_nz_master,z_global_master,&
       & nscatt_types_master,scatt_master_info_master,global_size_dists_master)
  !
  !
  !---------------------------------------
  ! See how many high res points need to 
  ! be averaged in the horizontal
  !---------------------------------------
  !
  h_res_master=(y2_master-y1_master)/scene_ny_master
  nh_to_avg=int(h_res/h_res_master+0.001)
  !
  if (nh_to_avg.lt.1) then
     call write_warning('Warning h_res too small')
     call write_warning('h_res ==> h_res_master')
     nh_to_avg=1
  endif
  !
  !
  !----------------------------
  ! What will the size of the
  ! output vertical array be ?
  !----------------------------
  !
  !
  io=1
  j=1
  !
  do i=2,scene_nz_master
     if ((z_global_master(i)-z_global_master(io)).ge.v_res*0.999) then
           j=j+1
           io=i
     endif
  enddo
  !
  scene_nz_out=j
  !
  !------------------
  ! Allocate storage
  !------------------
  !
  Allocate(z_global_out(1:scene_nz_out))
  Allocate(global_size_dists_out(nscatt_types_master,scene_nz_out))
  do i=1,scene_nz_out
     do j=1,nscatt_types_master
        Call Nullify_size_dist(global_size_dists_out(j,i))
     enddo
  enddo
  !
  !----------------------------
  ! Build the output vertical
  ! array
  !---------------------------
  !
  j=1
  z1=0.0
  z2=0.0
  n=0
  !
  z_global_out=0.0
  !
  z_global_out(1)=z_global_master(1)    ! the surface point
  io=1
  j=2
  !
  do i=2,scene_nz_master
     if ((z_global_master(i)-z_global_master(io)).lt.v_res*1.0001) then
        z_global_out(j)=z_global_out(j)+z_global_master(i)
        n=n+1
     else
        if (n.ge.1) then
           z_global_out(j)=z_global_out(j)/n
           j=j+1
           io=i
           n=0
        else
           z_global_out(j)=z_global_master(i)
           io=i
           j=j+1
           n=0
        endif
     endif
  enddo
  !
  max_n_h=2*scene_nx_master
  !
  scene_nx_out=int(min(scene_nx_master,max_n_h)/nh_to_avg)
  scene_ny_out=int(min(scene_ny_master,max_n_h)/nh_to_avg)
  !
  allocate(work_gasses(1:n_gasses_master))
  !
  x1_out=x1_master
  x2_out=x1_out+scene_nx_out*h_res
  !
  y1_out=y1_master
  y2_out=y1_out+scene_ny_out*h_res
  !
  z1_out=z1_master
  z2_out=z2_master
  !
  Allocate(data_point_out(1:scene_nx_out,1:scene_ny_out,1:scene_nz_out))
  Allocate(n_hits(1:scene_nx_out,1:scene_ny_out,1:scene_nz_out))
  Allocate(n_hits_st(1:scene_nx_out,1:scene_ny_out,1:scene_nz_out,1:nscatt_types_master))
  Allocate(n_hits_global(1:scene_nz_out,1:nscatt_types_master))
  !
  n_hits=0
  n_hits_st=0
  n_hits_global=0
  !
  do ix=1,scene_nx_out
     do iy=1,scene_ny_out
        do iz=1,scene_nz_out
           call Nullify_atmos_point(data_point_out(ix,iy,iz))
           Call allocate_atmos_point(data_point_out(ix,iy,iz),n_gasses_master,nscatt_types_master)
           !
           data_point_out(ix,iy,iz)%x=0.0
           data_point_out(ix,iy,iz)%y=0.0
           data_point_out(ix,iy,iz)%z=0.0
           data_point_out(ix,iy,iz)%lat=0.0
           data_point_out(ix,iy,iz)%long=0.0
           data_point_out(ix,iy,iz)%t=0.0
           data_point_out(ix,iy,iz)%p=0.0
           data_point_out(ix,iy,iz)%rho=0.0
           data_point_out(ix,iy,iz)%rh=0.0
           data_point_out(ix,iy,iz)%X_vol=0.0
           data_point_out(ix,iy,iz)%x_vel=0.0
           data_point_out(ix,iy,iz)%y_vel=0.0
           data_point_out(ix,iy,iz)%z_vel=0.0
           data_point_out(ix,iy,iz)%turb_dv=0.0
           !
        enddo
     enddo
  enddo
  !
  !-----------
  ! Main loop 
  !-----------
  !
  do i=1,min(scene_nx_master,max_n_h)
  !do i=1,scene_nx_master
     io=1+int(i/nh_to_avg-0.01)
     !
     call write_prog(i,min(scene_nx_master,max_n_h))
     !
     maxio=io
     !
     do j=1,min(scene_ny_master,max_n_h)
     !do j=1,scene_ny_master
        !
        jo=1+int(j/nh_to_avg-0.01)
        maxjo=jo
        !
        if (allocated(data_point_master)) then
           deallocate(data_point_master)
        endif
        !
        Allocate(data_point_master(1:scene_nz_master))
        do iz=1,scene_nz_master
           Call Nullify_atmos_point(data_point_master(iz))
        enddo
        !
        Call Read_uff_column(uff_in,ascii_or_bi_master,3,i,j,scene_nx_master,scene_ny_master,&
             & scene_nz_master,data_point_master,&
             & n_gasses_master,nscatt_types_master,scatt_master_info_master)
        !
        !
        !----------------
        ! Average things
        !----------------
        !
        kl=1
        ko=1
        do k=1,scene_nz_master
           !
           if (k.gt.1) then
              if ((z_global_master(k)-z_global_master(kl)).ge.v_res*0.9999) then 
                 kl=k
                 ko=ko+1
              endif
           else
              ko=1
           endif
           !
           n_hits(io,jo,ko)=n_hits(io,jo,ko)+1
           !
           !
           data_point_out(io,jo,ko)%x=data_point_out(io,jo,ko)%x+data_point_master(1)%x
           data_point_out(io,jo,ko)%y=data_point_out(io,jo,ko)%y+data_point_master(1)%y
           !
           !
           if (ko==1) then
              data_point_out(io,jo,ko)%z=-v_res/2.0
           else
              data_point_out(io,jo,ko)%z=data_point_out(io,jo,ko)%z+data_point_master(k)%z
           endif
           !
           data_point_out(io,jo,ko)%lat=data_point_out(io,jo,ko)%lat+data_point_master(k)%lat
           data_point_out(io,jo,ko)%long=data_point_out(io,jo,ko)%long+data_point_master(k)%long
           data_point_out(io,jo,ko)%t=data_point_out(io,jo,ko)%t+data_point_master(k)%t
           data_point_out(io,jo,ko)%p=data_point_out(io,jo,ko)%p+data_point_master(k)%p
           data_point_out(io,jo,ko)%rho=data_point_out(io,jo,ko)%rho+data_point_master(k)%rho
           !           data_point_out(io,jo,ko)%rh=data_point_out(io,jo,ko)%rh+data_point_master(k)%rh
           data_point_out(io,jo,ko)%X_vol=data_point_out(io,jo,ko)%X_vol+data_point_master(k)%X_vol
           data_point_out(io,jo,ko)%x_vel=data_point_out(io,jo,ko)%x_vel+data_point_master(k)%x_vel
           data_point_out(io,jo,ko)%y_vel=data_point_out(io,jo,ko)%y_vel+data_point_master(k)%y_vel
           data_point_out(io,jo,ko)%z_vel=data_point_out(io,jo,ko)%z_vel+data_point_master(k)%z_vel
           data_point_out(io,jo,ko)%turb_dv=data_point_out(io,jo,ko)%turb_dv+data_point_master(k)%turb_dv
           !
           !----------
           ! Surface
           !----------
           !
           data_point_out(io,jo,ko)%Surface=data_point_master(1)%surface
           !
           !-------------------
           ! The global size 
           ! Dists (only once)
           !--------------------
           !
           if ((i.eq.1).and.(j.eq.1)) then
              do isc=1,nscatt_types_master
                 !
                 if (Associated(global_size_dists_master(isc,k)%N_bin)) then
                    !
                    if (.not.(Associated(global_size_dists_out(isc,ko)%N_bin))) then
                       Call allocate_size_dist(global_size_dists_out(isc,ko),&
                            & scatt_master_info_master(isc)%num_sizes_want)
                       global_size_dists_out(isc,k)%N_bin(:)=0.0  
                    endif
                    !
                    n_hits_global(ko,isc)=n_hits_global(ko,isc)+1
                    !
                    global_size_dists_out(isc,ko)%N_bin(:)=&
                         &global_size_dists_out(isc,ko)%N_bin(:)+&
                         &global_size_dists_master(isc,k)%N_bin(:)
                    !
                 endif
              enddo
           endif
           !
           !
           !----------------
           ! The size dists
           !----------------
           !
           do isc=1,nscatt_types_master
              !
              if (Associated(data_point_master(k)%size_dists(isc)%N_bin)) then
                 !
                 if (.not.(Associated(data_point_out(io,jo,ko)%size_dists(isc)%N_bin))) then
                    Call allocate_size_dist(data_point_out(io,jo,ko)%size_dists(isc),&
                         & scatt_master_info_master(isc)%num_sizes_want)
                    data_point_out(io,jo,ko)%size_dists(isc)%N_bin=0.0
                 endif
                 !
                 n_hits_st(io,jo,ko,isc)=n_hits_st(io,jo,ko,isc)+1
                 !
                 data_point_out(io,jo,ko)%size_dists(isc)%N_bin=&
                      & data_point_out(io,jo,ko)%size_dists(isc)%N_bin+&
                      & data_point_master(k)%size_dists(isc)%N_bin
                 !
              endif
              !
           enddo
           !
        enddo
        !
        do k=1,scene_nz_master
           Call deallocate_atmos_point(data_point_master(k))
           Call nullify_atmos_point(data_point_master(k))
        enddo
        !
     enddo
     !
  enddo
  !
  call close_uff(3)
  !
  do io=1,scene_nx_out
     !
     do jo=1,scene_ny_out
        !
        do ko=1,scene_nz_out
           !
!           data_point_out(io,jo,ko)%x=data_point_out(io,jo,ko)%x/n_hits(io,jo,ko)
!           data_point_out(io,jo,ko)%y=data_point_out(io,jo,ko)%y/n_hits(io,jo,ko)
           !
           data_point_out(io,jo,ko)%x=io*h_res-h_res/2.0
           data_point_out(io,jo,ko)%y=jo*h_res-h_res/2.0
!
           if (ko.gt.1) then
              data_point_out(io,jo,ko)%z=data_point_out(io,jo,ko)%z/n_hits(io,jo,ko)
           endif
           !
           data_point_out(io,jo,ko)%lat=data_point_out(io,jo,ko)%lat/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%long=data_point_out(io,jo,ko)%long/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%t=data_point_out(io,jo,ko)%t/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%p=data_point_out(io,jo,ko)%p/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%rho=data_point_out(io,jo,ko)%rho/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%X_vol=data_point_out(io,jo,ko)%X_vol/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%x_vel=data_point_out(io,jo,ko)%x_vel/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%y_vel=data_point_out(io,jo,ko)%y_vel/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%z_vel=data_point_out(io,jo,ko)%z_vel/n_hits(io,jo,ko)
           data_point_out(io,jo,ko)%turb_dv=data_point_out(io,jo,ko)%turb_dv/n_hits(io,jo,ko)
           !
           if ((io==1).and.(jo==1)) then
              do isc=1,nscatt_types_master
                 if (n_hits_global(ko,isc).gt.0) then
                    global_size_dists_out(isc,ko)%N_bin=&
                         global_size_dists_out(isc,ko)%N_bin/n_hits_global(ko,isc)
                 endif
              enddo
           endif
           !
           do isc=1,nscatt_types_master
              if (Associated(data_point_out(io,jo,ko)%size_dists(isc)%N_bin)) then
                 data_point_out(io,jo,ko)%size_dists(isc)%N_bin=&
                      &data_point_out(io,jo,ko)%size_dists(isc)%N_bin/n_hits_st(io,jo,ko,isc)
              endif
           enddo
           !
        enddo
        !
        ascii_or_bi_out=3
        !
        call write_bi_di(uff_out)
        !
        !if (ascii_or_bi_out.eq.0) then 
        !   call write_bi(uff_out,10)
        !else
        !call write_ascii(uff_out,10)
        !endif
        !
     enddo
  enddo
  !
  call write_info('Uff_averager Finished')
  !
  Close(4)
  Do isc=1,nscatt_types_out
     Close(4+isc+10)
  Enddo
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
    Character(len=256)                                :: string1,string2
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
          !write(*,*) isc,size_dist_uff_files(isc)
          !pause
          !
       enddo
       !
       !----------------------------
       ! Open and write the header
       ! file
       !----------------------------
       !
       open(unit=3,file=header_uff_file,status='unknown',iostat=status,form='formatted')
       !
       !-----------------------
       ! Register file names
       !-----------------------
       !
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
          Write(3,*) gasses_master(i)
       Enddo
       Write(3,*) '-------------------------------------------------'
       Write(3,*) '           nx,ny,nz (scene dimensions)'
       Write(3,*) scene_nx_out,scene_ny_out,scene_nz_out
       Write(3,*) '--------------x1,y1,lat1,long1,z1----------------'
       Write(3,'(5f9.4)') x1_out,y1_out,lat1_out,long1_out,data_point_out(io,jo,1)%z
       Write(3,*) '--------------x2,y2,lat2,long2,z2----------------'
       Write(3,'(5f9.4)') x2_out,y2_out,lat2_out,long2_out,data_point_out(io,jo,scene_nz_out)%z
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
          write(3,*) data_point_out(io,jo,iz)%z 
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
       id=4
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
!          write(*,*) isc,worki,size_dist_uff_files(isc)
          open(unit=id+id2*10+isc,file=size_dist_uff_files(isc),status='replace',iostat=status,&
               & form='unformatted',&
               & access='direct',recl=worki)
!	          write(*,*) id+id2*10+isc,size_dist_uff_files(isc)
       enddo
       !
       init=.true.
       !
       call write_uff_direct(id,id2,io,jo,0,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
            & data_point_out(io,jo,1),init)
       !
       init=.false.
       do iz=1,scene_nz_out
          call write_uff_direct(id,id2,io,jo,iz,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
               & data_point_out(io,jo,iz),init)    
       enddo
       !
    Else
       !
       id=4
       id2=1
       init=.false.
       call write_uff_direct(id,id2,io,jo,0,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
            & data_point_out(io,jo,1),init)   
       !
       do iz=1,scene_nz_out
          call write_uff_direct(id,id2,io,jo,iz,n_gasses_out,nscatt_types_out,scatt_master_info_master,&
               & data_point_out(io,jo,iz),init)    
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
    If ((io.Eq.1).And.(jo.Eq.1)) Then
       !
       !
       Open(unit=id,file=filename,status='unknown',iostat=status,form='formatted',access='sequential')
       !
       !----------------------------
       ! Write out the names of the 
       ! Scattering master lists 
       ! this scene requires
       !----------------------------
       !
       zero=0.0
       !
       Write(id,'(a)') '-----------------------------------------'
       Write(id,*) nscatt_types_out,'  ,Number of scatter types'  
       Write(id,'(a)') '-----------------------------------------'
       Do i=1,nscatt_types_out
          Write(id,'(a80)') scatt_list_names_master(i)   ! Master list file for each scatterer type
       Enddo
       Write(id,'(a)') '------------------------------------------'
       Write(id,*) n_gasses_out,'  ,Number of atmospheric gasses'
       Write(id,*) '------------------------------------------'
       Do i=1,n_gasses_out
          Write(id,*) gasses_master(i) ! name of gass
       Enddo
       Write(id,'(a)') '-------------------------------------------------'
       Write(id,'(a)') '           nx,ny,nz (scene dimensions)'
       Write(id,*) scene_nx_out,scene_ny_out,scene_nz_out-1
       Write(id,'(a)') '--------------x1,y1,lat1,long1,z1----------------'
       Write(id,'(5f9.4)') x1_out,y1_out,zero,zero,data_point_out(io,jo,1)%z 
       Write(id,'(a)') '--------------x2,y2,lat2,long2,z2----------------'
       Write(id,'(5f9.4)') x2_out,y2_out,zero,zero,data_point_out(io,jo,scene_nz_out-1)%z
       Write(id,'(a)') '-------------------------------------------------'
       Write(id,'(a)') 'Backgrnd scattering info (function of z only)'
       Write(id,'(a)') '-------------------------------------------------'
       !
       Do iz=1,scene_nz_out-1
          write(id,*) data_point_out(io,jo,iz)%z 
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
    Do iz=1,scene_nz_out-1
       If (iz==1) Then 
          Write(id,*) '-------------------------------------------------'
          Write(id,'(a20,1x,i6,1x,i6)') 'Info for column     ',io,jo
          Write(id,*) '-------------------------------------------------'
          Write(id,'(a20,1x,i6,1x,i6)') 'Surface type and subtype',data_point_out(io,jo,1)%Surface,int(zero)
          Write(id,*) '-------------------------------------------------'
          !
       endif
       !
       Write(id,100) data_point_out(io,jo,iz)%x,data_point_out(io,jo,iz)%y,&
            &data_point_out(io,jo,iz)%lat,data_point_out(io,jo,iz)%long,&
            &data_point_out(io,jo,iz)%z,&
            &data_point_out(io,jo,iz)%t,&
            &data_point_out(io,jo,iz)%p,&
            &data_point_out(io,jo,iz)%rho,&
            &data_point_out(io,jo,iz)%n_gasses
       Do i=1,n_gasses_out
          Write(id,200,advance='NO') data_point_out(io,jo,iz)%X_vol(i)
       Enddo
       i=n_gasses_out
       Write(id,200) data_point_out(io,jo,iz)%X_vol(i)
       Write(id,300) data_point_out(io,jo,iz)%x_vel,&
            & data_point_out(io,jo,iz)%y_vel,&
            & data_point_out(io,jo,iz)%z_vel,&
            & data_point_out(io,jo,iz)%turb_dv
       Write(id,'(i6)') nscatt_types_out
       Do i=1,nscatt_types_out
          If (Associated(data_point_out(io,jo,iz)%size_dists(i)%N_bin)) Then
             j=data_point_out(io,jo,iz)%size_dists(i)%n_sizes
             Write(id,'(i6,1x,i6)') i,j
             do ii=1,j
                write(id,200,advance='NO') data_point_out(io,jo,iz)%size_dists(i)%N_bin(ii)*1.0e-6
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
    If ((io.Eq.1).And.(jo.Eq.1)) Then
       !
       !
       Open(unit=id,file=filename,status='unknown',iostat=status,form='unformatted',access='sequential')
       !
       zero=0.0
       !
       !----------------------------
       ! Write out the names of the 
       ! Scattering master lists 
       ! this scene requires
       !----------------------------
       !
       Write(id) 'bbbbbbbb---------------------------------'
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
       Write(id) scene_nx_out,scene_ny_out,scene_nz_out-1
       Write(id) '--------------x1,y1,lat1,long1,z1----------------'
       Write(id) x1_out,y1_out,zero,zero,data_point_out(io,jo,1)%z
       Write(id) '--------------x2,y2,lat2,long2,z2----------------'
       Write(id) x2_out,y2_out,zero,zero,data_point_out(io,jo,scene_nz_out-1)%z
       Write(id) '-------------------------------------------------'
       Write(id) 'Backgrnd scattering info (function of z only)'
       Write(id) '-------------------------------------------------'
       !
       Do iz=1,scene_nz_out-1
          write(id) data_point_out(io,jo,iz)%z 
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
    Do iz=1,scene_nz_out-1
       If (iz==1) Then 
          Write(id) 'Info for column----------------------------------'
          Write(id) io,jo
          Write(id) '-------------------------------------------------'
          !
          write(id) data_point_out(io,jo,1)%Surface,int(zero)
          !
          Write(id) '-------------------------------------------------'
          !
       endif
       !
       Write(id) data_point_out(io,jo,iz)%x,data_point_out(io,jo,iz)%y,&
            &data_point_out(io,jo,iz)%lat,data_point_out(io,jo,iz)%long,&
            &data_point_out(io,jo,iz)%z,&
            &data_point_out(io,jo,iz)%t,&
            &data_point_out(io,jo,iz)%p,&
            &data_point_out(io,jo,iz)%rho,&
            &data_point_out(io,jo,iz)%n_gasses
       Write(id) data_point_out(io,jo,iz)%X_vol
       i=n_gasses_out
       Write(id) data_point_out(io,jo,iz)%x_vel,&
            & data_point_out(io,jo,iz)%y_vel,&
            & data_point_out(io,jo,iz)%z_vel,&
            & data_point_out(io,jo,iz)%turb_dv
       Write(id) nscatt_types_out 
       Do i=1,nscatt_types_master
          If (Associated(data_point_out(io,jo,iz)%size_dists(i)%N_bin)) Then
             j=data_point_out(io,jo,iz)%size_dists(i)%n_sizes
             Write(id) i,j
             write(id) data_point_out(io,jo,iz)%size_dists(i)%N_bin(1:j)*1.0e-6
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
    error_str='Usage: uff_averager header_input.uff output.uff hor_res[km] v_res[km] '
    !
    status=0
    !
    nargs=iargc()  
    If (nargs.Ne.4) Then
       status=1
       error_str2='Wrong number of arguments'
    Endif
    !
    If (status == 0) Then
       Call getarg(1,arg_str)
       error_str2='Error in infilename specification'
       Read(arg_str,'(a130)',iostat=status,err=100) uff_in
    Endif
    !
    If (status == 0) Then
       Call getarg(2,arg_str)
       error_str2='Error in infilename specification'
       Read(arg_str,'(a130)',iostat=status,err=100) uff_out
    Endif
    !
    If (status == 0) Then
       Call getarg(3,arg_str)
       error_str2='Error in hor_res specification'
       Read(arg_str,'(f)',iostat=status,err=100) h_res
    Endif
    !
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str2='Error in hor_res specification'
       Read(arg_str,'(f)',iostat=status,err=100) v_res
    Endif
    !
!    If (status == 0) Then
!       Call getarg(5,arg_str)
!       error_str2='Error in maximum number of horizontal points to use as input'
!       Read(arg_str,'(i)',iostat=status,err=100) max_n_h
!    Endif
!    !
!    If (status == 0) Then
!       Call getarg(6,arg_str)
!       error_str2='1==> ascii output, 0==> binary output'
!       Read(arg_str,'(i)',iostat=status,err=100) ascii_or_bi_out
!    Endif
    !
100 If (status.Ne.0) Then 
       call Write_info(error_str2)
       call Write_error(error_str)
    Endif
    !
  End Subroutine get_data
  !
End program uff_averager



