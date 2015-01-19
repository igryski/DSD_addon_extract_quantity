!!
!! *PROGRAM*  extract_env
!!
!! *USAGE*: extract_env input.uff env_data.nc x_start(km) x_finish(km)'
!!            ' y_start(km) y_finish(km) hor_res(km) z1(km) z2(km) vert_res(km)'
!!
!! @version ecsim 1.1.2
!!
!! *SRC_FILE*
!!
!! tools/product_tools/src/extract_env.f90
!!
!! *LAST CHANGES*
!! 
!! -Nov 30, 2011: D.D. Improved checking against gass types
!! -Nov 22, 2011: D.D. Added missing atmos_point deallocation statments..fixed running out of memory on large scenes.
!!- Nov 12, 2010, D.D.: Fixed ny_ins that should have been ny_ins. Lead to array out of bounds when requested field was not square
!!- Jan 11, 2008, D.D.: Added some cosmetic write_info statements
!!- Dec 18, 2007, R.V : Added 'header_' to UFF filename
!! -Nov 16, 2007, D.D . Added this header
!!
!! *DESCRIPTION* 
!!
!! This program reads a 3-D domain from a UFF file and 
!! extracts various averaged quantities and outputs them to the specified output file.
!!
!! -input.uff         : UFF header file file to read from.
!! -env_data.nc       : Ncdf file to write data to.
!! -x_start           : Starting X position in km.
!! -x_finish          : Final X position in km.
!! -y_start           : Starting Y position in km.
!! -y_finish          : Final Y position in km.
!! -hor_res           : Horizontial output resolution in km.
!! -z1                : Lower altitude bound in km.
!! -z2                : Upper altitude bound in km.
!! -vert_res          : Vertical output resolution in km.
!!
!! If the incorrect number of commandline arguments are given then
!! an error message is printed along with a help message.
!!
!  Copyright KNMI for ESA
!
!  This source code is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
Program extract_env
  !
  Use scene_creator_types
  Use data_types
  Use physical_parameters
  Use read_uff
  Use Ray_SCATT
  Use write_messages
  !
  Implicit None
  !
  !
  !------------------
  ! Input variables
  !------------------
  !
  Character(len=256)             :: outfilename
  Character(len=256)             :: infilename
  Real                           :: cnst,waves1
  Real                           :: z1_ins,z2_ins,dz_ins,buffer,dh_ins,hor_res
  !
  !----------------
  ! Main vaiables
  !----------------
  !
  Real                                      :: x_start,y_start,x_finish,y_finish
  Character(len=10), Dimension (:), Pointer :: gasses
  Integer                                   :: nscatt_types,n_gasses,scene_nx,scene_ny,scene_nz
  !
  Character(len=256), Dimension (:), Pointer         :: scatt_list_names
  Type(scatt_prop_master),Dimension(:),Allocatable   :: scatt_master_info ! structure containg scattering list info
  Type(scatterer_info),Dimension(:),Allocatable      :: rad_scatt_info    ! where the data is really stored
  Type(scatterer_info_pol),Dimension(:),Allocatable  :: lid_scatt_info    ! where the data is really stored
  !
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists
  Real,Dimension(:),Allocatable                      :: z_global
  !
  Type(atmos_point),Dimension(:),Allocatable :: data_column                ! structure containing column data
  !
  Real,Dimension(:),Allocatable              :: x,y
  Real,Dimension(:,:,:),Allocatable          :: P_3grid
  Real,Dimension(:,:,:),Allocatable          :: T_3grid
  Real,Dimension(:,:,:),Allocatable          :: Wv_3grid
  Real,Dimension(:,:,:),Allocatable          :: O3_3grid
  Real,Dimension(:,:,:),Allocatable          :: u_3grid
  Real,Dimension(:,:,:),Allocatable          :: v_3grid
  Real,Dimension(:,:,:),Allocatable          :: w_3grid
  !
  Real,Dimension(:),Allocatable              :: z_ins ! instrument resolution altitude vector (km)
  Integer                                    :: nz_ins,nx_ins,ny_ins
  !
  !-------------------------
  ! Misc working variables 
  !-------------------------
  !
  Integer                              :: status,i,ix,iy,iz,isc,j,k,iozone
  Integer                              :: ix1,ix2,iy1,iy2,iix,iiy
  Type(scatterer_info)                 :: scatt_info_temp
  Integer                              :: ix_out_min,ix_out_max
  Integer                              :: iy_out_min,iy_out_max
  real                                 :: x1,y1,lat1,long1,z1
  real                                 :: x2,y2,lat2,long2,z2
  Integer                              :: same_res
  Integer                              :: ascii_or_bi
  Real                                 :: freq
  Character(len=100)                   :: title, nc_title,units,plot_title
  character(len=200)                   :: ecsim_home, scatt_lib,error_str
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
       Character(len=*),Intent(out)           :: arg
     End Subroutine getarg
     !   
  End Interface
  !
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
  If (status.eq.1) Then 
     error_str='Command line arguments not present or incomplete'
     goto 200
  Endif
  !
  !----------------------------------
  ! Add ECSIM_HOME to the input and 
  ! output files
  !----------------------------------
  !
  loc = index(infilename,'/',.true.)
  char1 = infilename(1:loc)
  char2 = infilename(loc+1:)
  infilename  = trim(adjustl(char1))//'header_'//trim(adjustl(char2))
  infilename  = trim(adjustl(ecsim_home))//infilename
  outfilename = trim(adjustl(ecsim_home))//outfilename
  !
  !---------------------------------------------------
  ! Open and read the the UFF file header information
  !---------------------------------------------------
  !
  Call write_info('Reading uff header file: '//infilename)
  !
  Call read_uff_header_info(infilename,nscatt_types,scatt_list_names,&
       & n_gasses,gasses,scene_nx,scene_ny,scene_nz,&
       & x1,y1,lat1,long1,z1,&
       & x2,y2,lat2,long2,z2,ascii_or_bi)
  !
  !-----------------------------------------------------
  ! Build the instrument resolution verticle coordinate
  !-----------------------------------------------------
  !
  same_res=0
  !
  if (dz_ins.le.0) then 
     call write_error( 'Must specifiy a positive dz_ins !')
     status=1
  else if (z1_ins.lt.z1) then
     z1_ins=z1
     same_res=0
  else if (z2_ins.gt.z2) then
     z2_ins=z2
     same_res=0
  endif
  !
  nz_ins=int((z2_ins-z1_ins)/dz_ins+0.5)
  Allocate(z_ins(1:nz_ins))
  do iz=1,nz_ins
     z_ins(iz)=(iz-1)*dz_ins+z1
  enddo
  !
  !  scene_nz=scene_nz-1
  !
  !--------------------------------------------------------
  ! Build the instrument resolution horizontial coordinate
  !---------------------------------------------------------
  !
  hor_res=(x2-x1)/(scene_nx)
  nx_ins=int((x2-x1)/hor_res+1)
  ny_ins=int((y2-y1)/hor_res+1)
  !
  Allocate(x(nx_ins))
  Allocate(y(ny_ins))
  !
  do i=1,nx_ins
     X(i)=x1+(x2-x1)/(nx_ins-1)*(i-1)
  enddo
  !
  ! D.D. Nov 12, 2010
  !
  do j=1,ny_ins
     Y(j)=y1+(y2-y1)/(ny_ins-1)*(j-1)
  enddo
  !
  !---------------------------------------------
  ! Now scan through the UFF and build up the 
  ! outputfile
  !---------------------------------------------
  !
  Allocate(data_column(1:nz_ins))
  Allocate(P_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  Allocate(T_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  Allocate(Wv_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  Allocate(O3_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  Allocate(u_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  Allocate(v_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  Allocate(w_3grid(1:nx_ins,1:ny_ins,1:nz_ins))
  !
  Do iz=1,nz_ins
     Call Nullify_atmos_point(data_column(iz))
  Enddo
  !
  !-----------------------------------
  ! Read the scattering master lists
  !-----------------------------------
  !
  Allocate(scatt_master_info(1:nscatt_types))
  !
  Call Write_info('=================================================================')
  Call Write_info('The Following scattering types are referenced in the input file')
  Call Write_info('=================================================================')
  Do i=1,nscatt_types
     write(error_str,*) i,':',trim(adjustl(scatt_list_names(i)))
     call write_info(error_str)
  enddo
  !
  !
  Do i=1,nscatt_types
     Call read_scatt_list_xml(5,scatt_list_names(i),scatt_master_info(i),status)
     if (status.ne.0) then  
        error_str='Error in reading scattering info'
        goto 200
     endif
  Enddo
  !
  ! What index is water
  !
  Do i=1,n_gasses
     If ((trim(adjustl(gasses(i)))=='h2o').Or.(trim(adjustl(gasses(i))).Eq.'H2O')) Then
        iwater=i
     Endif
  Enddo
  !
  ! What index is O3
  !
  Do i=1,n_gasses
     If ((trim(adjustl(gasses(i)))=='o3').Or.(trim(adjustl(gasses(i))).Eq.'O3')) Then
        iozone=i
     Endif
  Enddo

  !---------------------------------
  ! Read the uff_file and populate 
  ! Quantity_grid
  ! Process the continious
  ! part first then take any
  ! wrap around into account
  !---------------------------------
  !
  !
  Do i=1,nx_ins
     iix=int(((x(i)-x1))/hor_res+1)
     if (iix.lt.1) then
        iix=iix+scene_nx
     else if (iix.gt.scene_nx) then
        iix=iix-scene_nx
     endif
     !
     Do j=1,ny_ins
        iiy=int(((y(j)-y1))/hor_res+1)
        if (iiy.lt.1) then
           iiy=iiy+scene_ny
        else if (iiy.gt.scene_ny) then
           iiy=iiy-scene_ny
        endif
        !
        Call Read_uff_column_intpol(infilename,ascii_or_bi,3,iix,iiy,scene_nx,scene_ny,scene_nz,data_column,&
             & nz_ins,z_ins,n_gasses,nscatt_types,scatt_master_info)
        ! 
        do k=1,nz_ins
           P_3grid(i,j,k)=data_column(k)%p
           T_3grid(i,j,k)=data_column(k)%T
           Wv_3grid(i,j,k)=data_column(k)%X_vol(iwater)*mol_wt_H2O*(100.0*data_column(k)%p)/(Rg*data_column(k)%T)
           O3_3grid(i,j,k)=data_column(k)%X_vol(iozone)*mol_wt_O3*(100.0*data_column(k)%p)/(Rg*data_column(k)%T)
           u_3grid(i,j,k)=data_column(k)%x_vel
           v_3grid(i,j,k)=data_column(k)%y_vel
           w_3grid(i,j,k)=data_column(k)%z_vel
        enddo
        !
        ! D.D. Nov 22, 2011
        !
        do iz=1,nz_ins
           Call deallocate_atmos_point(data_column(iz))
        enddo
        !
     Enddo
  Enddo
  !
  if ((ascii_or_bi.ne.2).and.(ascii_or_bi.ne.3)) then
     Call Close_uff(3)
  endif
  !
  call write_results_ncdf(outfilename)
  !
  !
  call Write_info('******************FINISHED******************')
  
200 If (status.ne.0) then
     call write_error(error_str)
     call exit(1)
  endif
  ! 
Contains
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
    Character(len=160)               :: arg_str
    Character(len=190)               :: error_str
    Character(len=100)               :: error_str2
    !
    error_str2='Usage: extract_env input.uff env_data.nc x_start(km) x_finish(km)'//&
         &       ' y_start(km) y_finish(km) hor_res(km) z1(km) z2(km) vert_res(km)'
    !
    status=0
    !
    nargs=iargc()  
    If (nargs.ne.10) Then
       status=2
       error_str='Wrong number of arguments'
    Endif
    !
    If (status == 0) Then
       Call getarg(1,arg_str)
       error_str='Error in infilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) infilename
    Endif
    !
    If (status == 0) Then
       Call getarg(2,arg_str)
       error_str='Error in outfilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) outfilename
    Endif
    !
    If (status == 0) Then
       Call getarg(3,arg_str)
       error_str='Error in x starting value (km)'
       Read(arg_str,*,iostat=status,err=100) x_start
    Endif
    !
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str='Error in x ending value (km)'
       Read(arg_str,*,iostat=status,err=100) x_finish
    Endif
    !
    If (status == 0) Then
       Call getarg(5,arg_str)
       error_str='Error in x starting value (km)'
       Read(arg_str,*,iostat=status,err=100) y_start
    Endif
    !
    If (status == 0) Then
       Call getarg(6,arg_str)
       error_str='Error in x ending value (km)'
       Read(arg_str,*,iostat=status,err=100) y_finish
    Endif
    !
    If (status == 0) Then
       Call getarg(7,arg_str)
       error_str='Error in x ending value (km)'
       Read(arg_str,*,iostat=status,err=100) dh_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(8,arg_str)
       error_str='Error in instrument start altitude(km)'
       Read(arg_str,*,iostat=status,err=100) z1_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(9,arg_str)
       error_str='Error in instrument stop altitude(km)'
       Read(arg_str,*,iostat=status,err=100) z2_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(10,arg_str)
       error_str='Error in desired vert. resolution (km)'
       Read(arg_str,*,iostat=status,err=100) dz_ins
    Endif
    !
100 If (status.ne.0) Then 
       call Write_error(error_str)
       call Write_error(error_str2)
    Endif
    !
  End Subroutine get_data
  !
  Subroutine write_results_ncdf(filename)
    !
    use typeSizes
    use netcdf
    use ncdf_utilities
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename

    Integer                        :: i,j
    !
    Integer :: ncid, status
    CHARACTER(len=200)             :: error_str
    !
    integer :: X_dim,Y_dim,Z_dim
    integer :: XDistId, YDistId, ZDistId
    integer :: Pid,Tid,wvid,o3id,uid,vid,wid
    !

    ! Assume the file does not exist;
    ! If you want to check, you need to use the nf90_open function
    status = nf90_create(filename, 0, ncid)
    if (status /= 0) then
       error_str='error in nf90_create'
       return
    endif
    !    
    ! Defining dimensions
    status = nf90_def_dim(ncid, "nx", nx_ins, X_dim)    
    if (status /= 0) then 
       error_str='error in nf90_def_dim: x_scene'
       return
    endif
    status = nf90_def_dim(ncid, "ny", ny_ins, Y_dim)    
    if (status /= 0) then 
       error_str='error in nf90_def_dim: y_scene'
       return
    endif
    status = nf90_def_dim(ncid, "nz", nZ_ins,  Z_dim)    
    if (status /= 0) then 
       error_str='error in nf90_def_dim: altitude'
       return
    endif
    !
    ! Defining variables
    status = nf90_def_var(ncid, "x", NF90_FLOAT, (/X_dim/), XDistId)
    if (status /= 0) then 
       error_str='error in nf90_def_var1'
       return
    endif
    status = nf90_def_var(ncid, "y", NF90_FLOAT, (/Y_dim/), YDistId)
    if (status /= 0) then 
       error_str='error in nf90_def_var2'
       return
    endif
    status = nf90_def_var(ncid, "z", NF90_FLOAT, (/Z_dim/), ZDistId)
    if (status /= 0) then 
       error_str='error in nf90_def_var2'
       return
    endif
    !
    call my_nf90_def_var(ncid, "Pressure", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), pId,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    call my_nf90_def_var(ncid, "Temperature", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), TId,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    call my_nf90_def_var(ncid, "Water_vapor", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), WvId,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    call my_nf90_def_var(ncid, "O3", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), O3Id,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    call my_nf90_def_var(ncid, "u", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), uId,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    call my_nf90_def_var(ncid, "v", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), vId,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    call my_nf90_def_var(ncid, "w", NF90_FLOAT, (/X_dim,Y_dim,Z_dim/), wId,status)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    !
    ! Add attributes
    !
    call add_attributes(ncid,pid,'pressure','mb','pressure in mb',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    call add_attributes(ncid,tid,'temperature','K','temperature in K',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    call add_attributes(ncid,wvid,'Waver_vapor','gm-3','Water Vapor in gm-3',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    call add_attributes(ncid,O3id,'O3','gm-3','O3 density gm-3',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    call add_attributes(ncid,uid,'u','ms-1','x wind in ms-1',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    call add_attributes(ncid,vid,'v','ms-1','y wind in ms-1',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    call add_attributes(ncid,wid,'w','ms-1','w wind in ms-1',status)
    If (status /= 0) Then
       error_str = 'Error in add_attributes: pressure)'
       return
    Endif
    !
    ! ---------------------------------
    ! End of defining part:
    Call my_nf90_enddef(ncid, status) 
    !
    status=nf90_put_var(ncid, XDistId, x)
    status=nf90_put_var(ncid, YDistId, y)
    status=nf90_put_var(ncid, ZDistId, z_ins)
    !
    Call my_nf90_put_var(ncid, pid,  P_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: pressure'
       return
    Endif
    !
    Call my_nf90_put_var(ncid, tid,  T_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: x_lr'
       return
    Endif
    !
    Call my_nf90_put_var(ncid, wvid,  Wv_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: x_lr'
       return
    Endif
    !
    Call my_nf90_put_var(ncid, O3id,  O3_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: x_lr'
       return
    Endif
    !
    Call my_nf90_put_var(ncid, uid,  u_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: x_lr'
       return
    Endif
    !
    Call my_nf90_put_var(ncid, vid,  v_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: x_lr'
       return
    Endif
    !
    Call my_nf90_put_var(ncid, wid,  w_3grid, status)      
    If (status /= 0) Then
       error_str = 'Error in my_nf90_put_var: x_lr'
       return
    Endif
    !
    status = nf90_close(ncid)
    if (status /= 0) then 
       error_str='Error in nf90_close'
       return
    endif
    !  
   Return
   !
 End Subroutine write_results_ncdf
 !
End Program extract_env
  
