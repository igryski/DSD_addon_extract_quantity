!! *PROGRAM* extract_quantity_hor
!! 
!! *USAGE*: extract_quantity_hor inputfile.uff outfile.nc qindex qindex2.
!! 
!!  @version ecsim 1.2.4
!!
!! *SRC_FILE*
!!
!! tools/product_tools/src/extract_quantity_hor.f90
!!
!! *LAST CHANGES*
!!
!! -Nov 22, 2011: D.D. Added missing atmos_point deallocation statments..fixed running out of memory on large scenes.
!! -Jan 29, 2009: D.D. 'Integrated Unattenuated_Reflectivity_94_GHz'==> 'Integrated_Unattenuated_Reflectivity_94_GHz' missing '_' was added. Did not affect running on my machine but on others it did (different ncdf versions ?)
!! -Oct 17, 2008: D.D. Fixed incorrect labeling of integrated reflectivity
!! -Oct 10, 2008: D.D. Fixed problem with mishandelling of mulitwavlength radar scatt lib info. (SPR 0000168)
!! -Mar 10, 2008: D.D. Modified the header and Copyright info
!! -Feb 17, 2008: D.D: Fixed  Y(:)=Y_grid(ix1,:) not Y(:)=Y_grid(1,:)
!! -Feb 04, 2008, D.D: Cleaned up exit code reporting
!! -Dec 18, 2007, R.V: Added 'header_' to UFF filename
!! -Nov 22, 2007, D.D: Added this header and fixed confusion between along track and x_scene in the output file
!! -Oct 04, 2007, D.D: Extended length of character strings in get_data and changed ESIM_HOME==>ECSIM_HOME
!!
!! *DESCRIPTION*
!!
!! This program reads an entire UFF file and extracts various
!! column integrated quantities and outputs them to the specified output file
!!
!! The following command line arguments are expected:
!! -inputfile.uff   : UFF file header to read from.
!! -outfile.nc   : netcdf file to write data to.
!! -qindex          : Index of quantity to extract (see table below).
!! -qindex2         : Index of scattering type (0=> use all scattering types).
!!
!! Current Choices For qindex are:
!!
!! -10 ===> Mass path    [g/m$^2$].
!! -11 ===> Average Reff [microns].
!! -13 ===> Optical Depth at 353 nm. 
!! -100 ==> Integrated idealized Radar reflectivity at 94 GHz [mm$^6$/m$^3$]*m.
!! 
!!
!! If the incorrect number of commandline arguments are given then
!! an error message is printed along with a help message.
!!
!! *OUTPUT*
!!
!! outfile.nc will contain the output information which may then
!! be plotted using plot_hor. The output is netcdf. 
!!
!! 
!
!  Copyright KNMI for ESA
!
!  This source code is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
Program extract_quantity_hor
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
  Character(len=256)              :: outfilename
  Character(len=256)              :: infilename
  Real                            :: z1_ins,z2_ins,dz_ins,buffer,hor_res
  Integer                         :: qindex,qindex2
  !
  !----------------
  ! Main vaiables
  !----------------
  !
  Real                                      :: x_start,y_start,x_finish,y_finish
  Character(len=256), Dimension (:), Pointer :: scatt_list_names
  Character(len=10), Dimension (:), Pointer :: gasses
  Integer                                   :: nscatt_types,n_gasses,scene_nx,scene_ny,scene_nz
  !
  Type(scatt_prop_master),Dimension(:),Allocatable   :: scatt_master_info ! structure containg scattering list info
  Type(scatterer_info),Dimension(:),Allocatable      :: rad_scatt_info    ! where the data is really stored
  Type(scatterer_info_pol),Dimension(:),Allocatable  :: lid_scatt_info    ! where the data is really stored
  !
  Type(size_dist),Dimension(:,:),Allocatable         :: global_size_dists
  Real,Dimension(:),Allocatable                      :: z_global
  !
  Type(atmos_point),Dimension(:),Allocatable :: data_column                ! structure containing column data
  !
  Real,Dimension(:),Allocatable              :: fall_vel
  Real,Dimension(:),Allocatable              :: Nsize
  Real,Dimension(:),Allocatable              :: Ze_vec
  Real,Dimension(:),Allocatable              :: ext_vec
  Real,Dimension(:),Allocatable              :: x,y
  Real,Dimension(:,:),Allocatable            :: Quantity
  !
  Real,Dimension(:,:),Allocatable            ::  X_grid,Y_grid
  Real,Dimension(:,:,:),Allocatable           :: Quantity_grid 
  !
  Real,Dimension(:),Allocatable              :: z_ins ! instrument resolution altitude vector (km)
  Integer                                    :: nz_ins
  !
  !-------------------------
  ! Misc working variables 
  !-------------------------
  !
  Integer                              :: status,i,ix,iy,iz,isc,j
  Integer                              :: ix1,ix2,iy1,iy2,iix,iiy
  Type(scatterer_info)                 :: scatt_info_temp
  Integer                              :: ix_out_min,ix_out_max
  Integer                              :: iy_out_min,iy_out_max
  real                                 :: x1,y1,lat1,long1,z1
  real                                 :: x2,y2,lat2,long2,z2,fac
  Integer                              :: same_res
  Integer                              :: ascii_or_bi
  Real                                 :: waves(2)
  Real                                 :: freq
  Character(len=80)                    :: title, nc_title,units,plot_title
  !
  character(len=256)                   :: ecsim_home, scatt_lib,error_str
  !
  ! UFF header update
  Integer :: loc
  Character(len=256) :: char1, char2
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
  !-----------------------------------------
  ! GET  environment variable ECSIM_HOME 
  !-----------------------------------------
  !
  call getenv('ECSIM_HOME', ecsim_home)
  call getenv('SCATT_LIB', scatt_lib)
  !
  !
  !-----------------------------------------
  ! Get the Arguments from the command line
  !-----------------------------------------
  !
  Call get_data(status)
  !
  If (status.eq.1) Then 
     call write_error('Command line arguments not present or incomplete')
     call exit(1)
  Endif
  !
  !----------------------------------
  ! Add ECSIM_HOME to the input and 
  ! output files
  !----------------------------------
  !
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
  Call read_uff_header_info(infilename,nscatt_types,scatt_list_names,&
       & n_gasses,gasses,scene_nx,scene_ny,scene_nz,&
       & x1,y1,lat1,long1,z1,&
       & x2,y2,lat2,long2,z2,ascii_or_bi)
  !
  x_start=x1
  x_finish=x2
  y_start=y1
  y_finish=y2
  !
  dz_ins=0.04
  !
  if (z2.gt.15.0) then
     z2=15.0
  endif
  !
  nz_ins=int((z2-z1)/dz_ins+0.5)
  Allocate(z_ins(1:nz_ins))
  do iz=1,nz_ins
     z_ins(iz)=(iz-1)*dz_ins+z1
  enddo
  !     
  !-----------------------------------
  ! Read the scattering master lists
  !-----------------------------------
  !
  Allocate(scatt_master_info(1:nscatt_types))
  Call Write_info('=================================================================')
  Call Write_info('The Following scattering types are referenced in the input file')
  Call Write_info('=================================================================')
  Do i=1,nscatt_types
     write(error_str,*) i,':',trim(adjustl(scatt_list_names(i)))
     call write_info(error_str)
  enddo

  Do i=1,nscatt_types
     Call read_scatt_list_xml(5,scatt_list_names(i),scatt_master_info(i),status)
     if (status.ne.0) then  
        error_str='Error in reading scattering info'
        goto 200
     endif
  Enddo
  !
  !
  if (qindex.ge.100) then
     !
     !---------------------------------------
     ! Read the radar scatering information
     !---------------------------------------
     !
     Allocate(rad_scatt_info(1:nscatt_types))
     !
     freq=94.0
     !
     Do i=1,nscatt_types
        !
        call find_rad_scatt_at_freq_sing(7,scatt_master_info(i),rad_scatt_info(i),freq)
        !
     enddo
     !
  else
     !
     !-------------------------------------
     ! Read the UV/optical/IR scattering
     ! Information
     !-------------------------------------
     !
     Allocate(lid_scatt_info(1:nscatt_types))
     !
     waves=0.353
     Do i=1,nscatt_types
        !
        call find_scatt_at_waves_pol(7,scatt_master_info(i),lid_scatt_info(i),&
             & 2,waves)
        !
     Enddo
     !
  endif
  !
  !---------------------------------------------
  ! Read the global scattering property section 
  ! of the UFF
  !---------------------------------------------
  !
  !
  ! First allocate storage
  !
  Allocate(global_size_dists(1:nscatt_types,1:nz_ins))
  Allocate(z_global(1:scene_nz))
  !
  Call read_global_size_dists_intpol(infilename,ascii_or_bi,scene_nz,z_global,nz_ins,z_ins,&
       &nscatt_types,scatt_master_info,global_size_dists)
  !
  !---------------------------------------------
  ! Now scan through the UFF and build up the 
  ! outputfile
  !---------------------------------------------
  !
  Allocate(data_column(1:nz_ins))
  !
  hor_res=(x2-x1)/(scene_nx)
  !
  !
  ix1=Int((x_start-x1)/hor_res)+1
  ix2=Int((x_finish-x1)/hor_res+0.01)
  iy1=Int((y_start-y1)/hor_res)+1
  iy2=Int((y_finish-y1)/hor_res+0.01)
  !
  !write(*,*) hor_res,ix1,ix2,iy1,iy2
  !
  Allocate(X_grid(ix1:ix2,iy1:iy2))
  Allocate(Y_grid(ix1:ix2,iy1:iy2))
  Allocate(Quantity_grid(ix1:ix2,iy1:iy2,1:nz_ins))
  !
  X_grid=0.0
  Y_grid=0.0
  Quantity_grid=0.0
  !
  Do iz=1,nz_ins
     Call Nullify_atmos_point(data_column(iz))
  Enddo
  !
  ix_out_min=ix1
  ix_out_max=ix2
  iy_out_min=iy1
  iy_out_max=iy2
  !
  !-------------------------
  ! Read and populate the
  ! output structure
  !-------------------------
  !
  Allocate(Quantity(1:scene_nx,1:scene_ny))
  !
  Do ix=1,scene_nx
     Do iy=1,scene_ny
        !
        !
        Call Read_uff_column_intpol(infilename,ascii_or_bi,3,ix,iy,scene_nx,scene_ny,scene_nz,data_column,&
             & nz_ins,z_ins,n_gasses,nscatt_types,scatt_master_info)
!        Call Read_uff_column_intpol(infilename,ascii_or_bi,3,-1,-1,scene_nx,scene_ny,scene_nz,data_column,&
!             & nz_ins,z_ins,n_gasses,nscatt_types,scatt_master_info)
        ! 
        X_grid(ix,iy)=data_column(1)%x
        Y_grid(ix,iy)=data_column(1)%y
        !
        !write(*,*) ix,iy,data_column(1)%x,data_column(1)%y
        !
        fac=0.0
        call populate(ix,iy,qindex,qindex2,fac)
        !
        Quantity(ix,iy)=0.0
        !
        do iz=1,nz_ins-1
           if (qindex.ne.11) then
              Quantity(ix,iy)=Quantity(ix,iy)+Quantity_grid(ix,iy,iz)*(z_ins(iz+1)-z_ins(iz))*1000.0
           else
              if (fac.gt.0.0) then
                 Quantity(ix,iy)=Quantity(ix,iy)+Quantity_grid(ix,iy,iz)/fac
              endif
           endif
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
  Call Close_uff(3)
  !
  !-----------------------------------
  ! Write the Results to the Output file
  !-----------------------------------
  !
  Allocate(x(1:scene_nx))
  Allocate(y(1:scene_ny))
  x=X_grid(:,iy1)
  y=Y_grid(ix1,:)
  !
  !call write_results(outfilename,size(x),size(y),x,y,Quantity,title)

  call write_results_ncdf(outfilename,size(x),size(y),x,y,x_grid,y_grid,Quantity,nc_title,title,units,plot_title)
  !
  !
  call Write_info( '******************FINISHED******************')
  call exit(0)
  
200 call Write_error(error_str)
  ! call Write_error(error_str2)
  call exit(1)

  ! 
Contains
  !
  !
  subroutine populate(ix,iy,qindex,qindex2,fac)
    !
    Implicit none
    !
    Integer,intent(in)             :: ix,iy,qindex,qindex2
    Real,intent(out)               :: fac
    !
    Integer                        :: irh,il,it,iz
    Real                           :: work1,work2
    Real,dimension(:),allocatable  :: Nsize
    !
    !
    fac=0.0
    !
    Do iz=1,nz_ins-1
       !
       ! Populate the output arrays
       !
       !
       if ((qindex.ge.10).and.(qindex.lt.100)) then ! We want some IR/Vis property of the scatterers
          !
          work1=0.0
          work2=0.0
          Do isc=1,nscatt_types
             if ((qindex2.lt.1).or.(qindex2==isc)) then
                If ((Associated(data_column(iz)%size_dists(isc)%N_bin)).Or.&
                     (Associated(global_size_dists(isc,iz)%N_bin))) Then
                   Allocate(Nsize(1:lid_scatt_info(isc)%n_sizes))
                   Nsize(:)=0.0
                   If (Associated(global_size_dists(isc,iz)%N_bin)) Then
                      Nsize=global_size_dists(isc,iz)%N_bin
                   Endif
                   If (Associated(data_column(iz)%size_dists(isc)%N_bin)) Then
                      Nsize=Nsize+data_column(iz)%size_dists(isc)%N_bin
                   Endif
                   !
                   if (qindex==10) then ! We want the total water/ice mass
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      work1=work1+Sum(lid_scatt_info(isc)%Vol(:)*&
                           &Nsize(:))*lid_scatt_info(isc)%density(it,irh)*&
                           &1.e-12  ! now in g/m^3
                      work2=1.0
                      title='Mass Path'
                      nc_title='Mass_Path'
                      units='g m^-2'
                      plot_title='Mass Path [g/m\u2\d]'                      
                      fac=1.0
                      !
                   else if (qindex==11) then ! We want Reff
                      work1=work1+Sum(lid_scatt_info(isc)%Vol*Nsize)*0.75
                      work2=work2+Sum(lid_scatt_info(isc)%Ac*Nsize)
                      !
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      fac=fac+Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      !
                      work1=work1*Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4
                      !
                      plot_title='R\deff\u[microns]'
                      title='reff'
                      nc_title='reff'
                      units='microns'
                      !
                   else if (qindex==12) then  ! Mean maximum dimension
                      work1=work1+Sum(lid_scatt_info(isc)%R*Nsize)*2.0
                      work2=work2+Sum(Nsize)
                      plot_title='<D> [microns]'
                      title='<D>'
                      nc_title='mean_max_dimension'
                      units='microns'
                   else if (qindex==13) then ! Visibale extinction
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      work1=work1+Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      work2=1.0
                      plot_title='Optical Depth 353 nm'
                      title='tau'
                      nc_title='optical_depth_353'
                      units='--'
                      !
                   endif
                   DeAllocate(Nsize)
                endIf
             endif
          enddo
          !
          if (work2==0.0) then
             Quantity_grid(ix,iy,iz)=0.0
          else
             Quantity_grid(ix,iy,iz)=work1/work2
          endif
       else if (qindex.ge.100) then
          if (qindex==100) then   ! we want to calculate the reflectivity
             !
             title='Integrated Unattenuated Reflectivity [94 GHz]'
             nc_title='Integrated_Unattenuated_Reflectivity_94_GHz'
             units='mm^6 m^-3 m'
             plot_title='Integrated Unattenuated Reflectivity [mm\U6\D/m\U3\Dm] [94 GHz]'!
             work1=0.0
             work2=0.0
             !
             Do isc=1,nscatt_types
                if ((qindex2.lt.1).or.(qindex2==isc)) then
                   If ((Associated(data_column(iz)%size_dists(isc)%N_bin)).Or.&
                        (Associated(global_size_dists(isc,iz)%N_bin))) Then
                      Allocate(Nsize(1:rad_scatt_info(isc)%n_sizes))
                      Nsize(:)=0.0
                      If (Associated(global_size_dists(isc,iz)%N_bin)) Then
                         Nsize=global_size_dists(isc,iz)%N_bin
                      Endif
                      If (Associated(data_column(iz)%size_dists(isc)%N_bin)) Then
                         Nsize=Nsize+data_column(iz)%size_dists(isc)%N_bin
                      Endif
                      !
                      ! Find the reflecivity-vs-size (function of temperature) 
                      !
                      Allocate(Ze_vec(1:rad_scatt_info(isc)%n_sizes))
                      Allocate(ext_vec(1:rad_scatt_info(isc)%n_sizes))
                      !
                      Call find_scatt_at_T_rad(data_column(iz)%t,&
                           & rad_scatt_info(isc),Ze_vec,ext_vec)
                      !
                      work1=work1+Sum(Nsize*Ze_vec)
                      work2=1.0
                      !
                      DeAllocate(Nsize)
                      DeAllocate(Ze_vec)
                      Deallocate(ext_vec)
                      !
                   endif
                endif
             enddo
             !
             if (work2==0.0) then
                Quantity_grid(ix,iy,iz)=0.0
             else
                Quantity_grid(ix,iy,iz)=work1/work2
             endif
             !
          endif
          !
       endif
    enddo
    return
  end subroutine populate
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
    Character(len=160)                :: arg_str
    Character(len=190)               :: error_str
    Character(len=100)               :: error_str2
    !
    error_str='Usage: extract_quantity_hor inputfilename outputfilename qindex qindex2'
    !
    status=0
    !
    nargs=iargc()  
    If (nargs.Ne.4) Then
       status=2
       error_str2='Wrong number of arguments'
    Endif
    !
    If (status == 0) Then
       Call getarg(1,arg_str)
       error_str2='Error in infilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) infilename
    Endif
    !
    If (status == 0) Then
       Call getarg(2,arg_str)
       error_str2='Error in outfilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) outfilename
    Endif
    !
    If (status == 0) Then
       Call getarg(3,arg_str)
       error_str2='Error in index of quantity to extract'
       Read(arg_str,*,iostat=status,err=100) qindex
    Endif
    !
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str2='Error in index of quantity to extract'
       Read(arg_str,*,iostat=status,err=100) qindex2
    Endif
    !
100 If (status.ne.0) Then 
      Call Write_error(error_str)
      Call Write_error(error_str2)
    Endif
    !
    If (status.eq.2) then
       Call Write_info('Current Choices For Quantity Index are:')
       Call Write_info('10 ===> mass path')
       call Write_info('11 ===> Average Reff')
       call Write_info('13 ===> Optical Depth at 353 nm ')
       call Write_info('100 ==> Integrated Idealized Radar 94 GHz reflectivity')
       call Write_info('----------------------------------')
       call Write_info('qindex2 ==> scatterer index (set to a specific scatterer index)')
       call Write_info('set to 0 to use all scatterer types present')
       call exit(1)
    endif
    
  End Subroutine get_data
  
end Program extract_quantity_hor
!
Subroutine write_results_ncdf(filename,nx,ny,x,y,x_grid,y_grid,Quant,nc_title,title,units,plot_title)
    !
    use typeSizes
    use netcdf
    Use write_messages
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: nx,ny
    Real,intent(in)                :: x(nx),y(ny)
    Real,intent(in)                :: x_grid(nx,ny)
    Real,intent(in)                :: y_grid(nx,ny)
    Real,intent(in)                :: Quant(nx,ny)
    Character(len=*),intent(in)    :: nc_title,title,units,plot_title
    !
    Integer                        :: i,j,xdir,ydir
    !
    Integer                        :: ncid, status
    !
    integer                          :: XDistId, YDistId,  QuantId,YsceneId,XsceneID
    Character(len=190)               :: error_str
    !
    ! Assume the file does not exist;
    ! If you want to check, you need to use the nf90_open function
    status = nf90_create(filename, 0, ncid)
    if (status /= 0) then 
       error_str='error in nf90_create'
       goto 400
    endif
    !
    ! Defining dimensions
    status = nf90_def_dim(ncid, "nx", nx, Xdir)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: GroundDist'
       goto 400
    endif
    status = nf90_def_dim(ncid, "ny", ny, Ydir)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: Altitude'
       goto 400
    endif
    
    !
    ! Defining variables
    status = nf90_def_var(ncid, "x_scene", NF90_FLOAT, (/Xdir,Ydir/), XsceneId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var1'
       goto 400
    endif
    status = nf90_def_var(ncid, "y_scene", NF90_FLOAT, (/Xdir,Ydir/), YsceneId)
    if (status /= 0) then
       error_str= 'error in nf90_def_var2'
       goto 400
    endif
    !
    status = nf90_def_var(ncid, "along_track", NF90_FLOAT, (/Xdir/), XDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var3'
       goto 400
    endif
    status = nf90_def_var(ncid, "cross_track", NF90_FLOAT, (/Ydir/), YDistId)
    if (status /= 0) then
       error_str= 'error in nf90_def_var4'
       goto 400
    endif
    !
    status = nf90_def_var(ncid, trim(nc_title), NF90_FLOAT, (/Xdir,Ydir/), QuantId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var5'
       goto 400
    endif
    !
    ! Defining attributes
    !
    status=nf90_put_att(ncid,XDistId,"long_name","Along track")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"plot_title","  [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,YDistId,"long_name","Cross track")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef2a'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2b'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"plot_title","Cross track [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,XsceneId,"long_name","X scene")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,XsceneId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,XsceneId,"plot_title","X scene  [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,YsceneId,"long_name","Y scene")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef2a'
       goto 400
    endif
    status=nf90_put_att(ncid,YsceneId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2b'
       goto 400
    endif
    status=nf90_put_att(ncid,YsceneId,"plot_title","Y scene [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,QuantId,"long_name",title)
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5a'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId,"units",units)
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5b'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId,"plot_title",plot_title)
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5c'
       goto 400
    endif
    !
    status = nf90_enddef(ncid)
    if (status /= 0) then 
       error_str= 'Error in nf90_enddef'
       goto 400
    endif
    ! 
    ! Writing variables
    status = nf90_put_var(ncid, XDistId, x(1:nx))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var1'
       goto 400
    endif
    status = nf90_put_var(ncid, YDistId, y(1:ny))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var2'
       goto 400
    endif
    status = nf90_put_var(ncid, XsceneId, x_grid)    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var3'
       goto 400
    endif
    status = nf90_put_var(ncid, YsceneId, y_grid)
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var4'
       goto 400
    endif
    status = nf90_put_var(ncid, QuantId, Quant(1:nx,1:ny))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var5'
       goto 400
    endif
    !
    ! Close the file
    !
    status = nf90_close(ncid)
    if (status /= 0) then 
       error_str= 'Error in nf90_close'
       goto 400
    endif
    
    !  
400 if (status.ne.0) then
       call write_error(error_str)
       call exit(1)
    endif
    
    Return
    
  End Subroutine write_results_ncdf
  !
