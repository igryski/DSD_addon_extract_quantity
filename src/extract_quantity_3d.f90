!!
!! *PROGRAM* extract_quantity_3d
!!
!! *USAGE*: extract_quantity_3d inputfile.uff outfile.nc x1 y1 x2 y2
!!               z1(km) z2(km) vert_res(km)  
!!               qindex qindex2 (wavelength)
!! wavelength is only needed for qindex=13,14 or 15
!!
!! @version ecsim 1.1.1
!!
!! *SRC_FILE*
!!
!! tools/product_tools/src/extract_quantity.f90
!!
!! *LAST CHANGES*
!! 
!!  -Nov 30, 2011: D.D. Improved checking against gass types
!!  -Nov 22, 2011: D.D. Added missing atmos_point deallocation statments..fixed running out of memory on large scenes.
!!  -Oct 09, 2008: Fixed wrong usage message (SPR 0000160)
!!  -Mar 10, 2008: Modified the header and Copyright info
!! - Feb 17, 2008: D.D. Fixed  Y(:)=Y_grid(ix1,:) not Y(:)=Y_grid(1,:)
!! - Feb 04, 2008, D.D. Cleaned up exit code reporting
!! - Dec 18, 2007, RV: Added 'header_' to UFF filename
!! - Oct 4, :D.D. Extended length of character strings in get_data and changed ESIM_HOME==>ECSIM_HOME
!!
!! *DESCRIPTION*
!!
!!
!! This program reads a 3-D domain from a UFF file and 
!! extracts various averaged quantities and outputs them to the specified output file.
!!
!! The following command line arguments are expected:
!!
!! -inputfile.uff : UFF file file to read from.
!! -outfile.nc    : ASCI file to write data to.
!! -x1            : Starting X position in km.
!! -x2            : Final X position in km.
!! -y1            : Starting Y position in km.
!! -y2            : Final Y position in km.
!! -z1            : Lower altitude bound in km.
!! -z2            : Upper altitude bound in km.
!! -vert_res      : Vertical output resolution in km.
!! -qindex        : Index of quantity to extract (see table below).
!! -qindex2       : Index of scattering type (0=> use all scattering types).
!! -wavelength    : Lidar wavelength (microns), only needed for qindex=13,14 or 15
!!
!!
!!
!! Current Choices For qindex are:
!!
!!
!! -1  ===> Temperature [K].
!! -2  ===> Pressure [mb].
!! -3  ===> Density [molecules/cm$^3$].
!! -4  ===> Xvol H2O.
!! -5  ===> Xvol O3.
!! -10 ===> mass [g/m$^3$].
!! -11 ===> Reff [um].
!! -12 ===> Mean maximum dimension [um].
!! -13 ===> Extinction at 353 nm [1/m].
!! -14 ===> Extinction/Backscatter Ratio at 353 nm [1/sr]. 
!! -15 ===> Backscatter at 353 nm [1/m/sr].
!! -16 ===> R'eff [microns].
!! -17 ===> No [#/cm$^3$].
!! -18 ===> R_a [microns].
!! -100 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (94 GHz).
!! -101 ==> Attenuation at 94 GHz [1/m].
!! -102 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (32 GHz).
!! -103 ==> Attenuation at 32 GHz [1/m].
!! -104 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (5 GHz).
!! -105 ==> Attenuation at 5 GHz [1/m].
!!
!!
!! If the incorrect number of commandline arguments are given then
!! an error message is printed along with a help message.
!!
!!
!! *OUTPUT*
!!
!! outfile.nc will contain the output information which may then
!! be plotted using plot_3d.
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
Program extract_quantity_3d
  !!
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
  Integer                        :: slice_x_or_y               !extract slice along 0 ==x,1==y
  Real                           :: cnst,waves1
  Real                           :: z1_ins,z2_ins,dz_ins,buffer,hor_res
  Integer                        :: qindex,qindex2
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
  Real,Dimension(:,:,:),Allocatable              :: Quantity
  !
  Real,Dimension(:,:),Allocatable            :: X_grid,Y_grid
  Real,Dimension(:,:,:),Allocatable          :: Quantity_grid 
  !
  Real,Dimension(:),Allocatable              :: z_ins ! instrument resolution altitude vector (km)
  Integer                                    :: nz_ins,nx_ins,ny_ins
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
  real                                 :: x2,y2,lat2,long2,z2
  Integer                              :: same_res
  Integer                              :: ascii_or_bi
  Real                                 :: waves(2),freq
  Character(len=100)                    :: title, nc_title,units,plot_title
  !
  character(len=256)                   :: ecsim_home, scatt_lib,error_str
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
  waves1=0.0
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
  Call read_uff_header_info(infilename,nscatt_types,scatt_list_names,&
       & n_gasses,gasses,scene_nx,scene_ny,scene_nz,&
       & x1,y1,lat1,long1,z1,&
       & x2,y2,lat2,long2,z2,ascii_or_bi)
  !
  !
  !-----------------------------------------------------
  ! Build the instrument resolution verticle coordinate
  !-----------------------------------------------------
  !
  !
  same_res=0
  !
  if (dz_ins.le.0) then 
     call write_error( 'Must specifiy a positive dz_ins !')
     call exit(1)
  else if (z1_ins.lt.z1) then
     z1_ins=z1
     same_res=0
  else if (z2_ins.gt.z2) then
     z2_ins=z2
     same_res=0
  endif
  !
  if ((dz_ins.gt.(z2-z1)/scene_nz).and.(same_res.eq.0)) then 
     call write_info( '*********************************')
     call write_info( 'Warning requested resolution')
     call write_info( 'Coarser than the file resolution')
     call write_info( '*********************************')
     call exit(2)
  endif
  !
  nz_ins=int((z2_ins-z1_ins)/dz_ins+0.5)
  Allocate(z_ins(1:nz_ins))
  do iz=1,nz_ins
     z_ins(iz)=(iz-1)*dz_ins+z1
  enddo
  !
  scene_nz=scene_nz-1
  !     
  !-----------------------------------
  ! Read the scattering master lists
  !-----------------------------------
  !
  Allocate(scatt_master_info(1:nscatt_types))
  !
  call write_info('=================================================================')
  call write_info( 'The Following scattering types are referenced in the input file')
  call write_info('=================================================================')
  Do i=1,nscatt_types
     write(error_str,*) i,':',trim(adjustl(scatt_list_names(i)))
     call write_info(error_str)
  enddo
  !
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
     if ((qindex==100).or.(qindex==101)) then 
        freq=94.0
     else if ((qindex==102).or.(qindex==103)) then 
        freq=32.0
     else
        freq=5.6
     endif
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
     if (waves1.le.0.0) then 
        waves1=0.353
     endif
     
     waves=waves1
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
  nx_ins=ix2-ix1+1
  ny_ins=iy2-iy1+1
  Allocate(X_grid(ix1:ix2,iy1:iy2))
  Allocate(Y_grid(ix1:ix2,iy1:iy2))
  Allocate(Quantity_grid(ix1:ix2,iy1:iy2,1:nz_ins))
  Allocate(X(ix1:ix2))
  Allocate(Y(iy1:iy2))
  !
  X_grid=0.0
  Y_grid=0.0
  Quantity_grid=0.0
  X=0.0
  Y=0.0
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
  !---------------------------------
  ! Read the uff_file and populate 
  ! Quantity_grid
  ! Process the continious
  ! part first then take any
  ! wrap around into account
  !---------------------------------
  !
  Do ix=1,scene_nx
     Do iy=1,scene_ny
        !
        Call Read_uff_column_intpol(infilename,ascii_or_bi,3,ix,iy,scene_nx,scene_ny,scene_nz,data_column,&
             & nz_ins,z_ins,n_gasses,nscatt_types,scatt_master_info)
        ! 
        ! Is the column one of the ones we want ?
        !
        if ((ix.ge.ix1).and.(ix.le.ix2).and.(iy.ge.iy1).and.(iy.le.iy2)) then
           !
           X_grid(ix,iy)=data_column(1)%x
           Y_grid(ix,iy)=data_column(1)%y
           !
           call populate(ix,iy,qindex,qindex2)
           !
        Endif
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
  !-------------------
  ! Read the buffer
  !-------------------
  !
  do ix=ix1,ix2
     !
     if (ix.lt.1) then
        iix=ix+scene_nx
     else if (ix.gt.scene_nx) then
        iix=ix-scene_nx
     else
        iix=ix
     endif
     !
     do iy=iy1,iy2
        if ((iy.lt.1).or.(iy.gt.scene_ny).or.(ix.lt.1).or.(ix.gt.scene_nx)) then
           !
           if (iy.lt.1) then
              iiy=iy+scene_ny
           else if(iy.gt.scene_ny) then
              iiy=iy-scene_ny
           else
              iiy=iy
           endif
           !
           Call Read_uff_column_intpol(infilename,ascii_or_bi,3,iix,iiy,scene_nx,scene_ny,&
                & scene_nz,data_column,&
                & nz_ins,z_ins,n_gasses,nscatt_types,scatt_master_info)
           !
           if ((ascii_or_bi.ne.2).and.(ascii_or_bi.ne.3)) then
              Call Close_uff(3)
           endif
           !
           call populate(ix,iy,qindex,qindex2)
           !
           if (ix.lt.1) then
              X_grid(ix,iy)=x_start+(ix-ix1)*hor_res+hor_res/2.0
           else if (ix.gt.scene_nx) then
              X_grid(ix,iy)=x_finish+(ix-ix2)*hor_res-hor_res/2.0
           endif
           !
           if (iy.lt.1) then
              Y_grid(ix,iy)=y_start+(iy-iy1)*hor_res+hor_res/2.0
           else if (iy.gt.scene_ny) then
              Y_grid(ix,iy)=y_finish+(iy-iy2)*hor_res-hor_res/2.0
           endif
           !
        endif
     enddo
  enddo
  Call Close_uff(3)
  !
  !-----------------------------------
  ! Write the Results to the Output file
  !-----------------------------------
  !
  !
  X(:)=X_grid(:,iy1)
  Y(:)=Y_grid(ix1,:)
  !
  !  call write_results(outfilename,ix2-ix1+1,iy2-iy1+1,nz_ins,x_grid,y_grid,z_ins,Quantity_grid,title)
  call write_results_ncdf(outfilename,nx_ins,ny_ins,nz_ins,x,y,z_ins,Quantity_grid,nc_title,title,units,plot_title,status,error_str)
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
  !
  subroutine populate(ix,iy,qindex,qindex2)
    !
    Implicit none
    !
    Integer,intent(in)             :: ix,iy,qindex,qindex2
    !
    Integer                        :: irh,il,it,iz,itheta,igass
    Real                           :: work1,work2
    Real,dimension(:),allocatable  :: Nsize
    Character(len=4)               :: waves_val
    !
    !write(*,*) ix,iy
    !write(*,*) '--------'
    !
    write(unit=waves_val, fmt='(I4)') nint(waves1*1000)
     
    Do iz=1,nz_ins
       !
       ! Populate the output arrays
       !
       work1=0.0
       work2=0.0
       !

       if ((qindex.gt.0).and.(qindex.lt.10)) then ! We want something not involving the scattering properties
          !
          work2=1.0
          if (qindex==1) then                     ! Extract temperature (K)
             work1=data_column(iz)%T       
             title="Temperature"       
             nc_title="Temperature"       
             units="K"
             plot_title="Temperature [K]"
          else if (qindex==2) then 
             title='Pressure'
             nc_title='Pressure'
             units="mb"
             plot_title="Pressure [mb]"
             work1=data_column(iz)%p               ! Pressure (mb)
          else if (qindex==3) then 
             plot_title='Density [mol/cm\u3\d]'
             nc_title='Density'
             units="mol cm^-3"
             title="Density"
             work1=data_column(iz)%rho             ! density molecules/cm^3
          else if (qindex==4) then                    
             igass=0
             title='X(H2O)'
             nc_title='H2O'
             units="-"
             plot_title="X(H\D2\UO)"
             Do i=1,n_gasses
                If ((trim(adjustl(gasses(i)))=='h2o').Or.(trim(adjustl(gasses(i))).Eq.'H2O')) Then
                   igass=i
                Endif
             Enddo
             If (igass==0) Then
                call Write_error('Water mixing ratio is not available !!')
                call write_error('Check the UFF file !')
                call exit(1)
             Endif
             !
             work1=data_column(iz)%X_vol(igass)            ! O3 volumne mixing ratio
          else if (qindex==5) then                    
             plot_title='X(O\D3\U)'
             nc_title='O3'
             units="-"
             title="X(O3)"
             igass=0
             Do i=1,n_gasses
                If ((trim(adjustl(gasses(i)))=='o3').Or.(trim(adjustl(gasses(i))).Eq.'O3')) Then
                   igass=i
                Endif
             Enddo
             If (igass==0) Then
                call Write_error('O3 mixing ratio is not available !!')
                call Write_error('Check the UFF file !')
                Call exit(1)
                return
             Endif
             !
             work1=data_column(iz)%X_vol(igass)            ! O3 volumne mixing ratio
             !
          else
             write(error_str,'(A20,i4)') 'Unknown Option ',qindex
             call Write_error(error_str)
             call write_error('Run program with no arguments to get help')
             call exit(1)
             return
             !
          endif
          !
          Quantity_grid(ix,iy,iz)=work1/work2
          !
       else if ((qindex.ge.10).and.(qindex.lt.100)) then ! We want some IR/Vis property of the scatterers
          !
          Do isc=1,nscatt_types
             if ((qindex2.lt.1).or.(qindex2==isc)) then
                If ((Associated(data_column(iz)%size_dists(isc)%N_bin)).Or.&
                     (Associated(global_size_dists(isc,iz)%N_bin))) Then
                   Allocate(Nsize(1:lid_scatt_info(isc)%n_sizes))
                   Nsize(:)=0.0
                   !
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
                      plot_title='Mass content g/m\u3\d'
                      nc_title='Mass_content'
                      units="g m^-3"
                      title="Mass content"
                   else if (qindex==11) then ! We want Reff
                      work1=work1+Sum(lid_scatt_info(isc)%Vol*Nsize)*0.75
                      work2=work2+Sum(lid_scatt_info(isc)%Ac*Nsize)
                      plot_title='R\deff\u [microns]'
                      nc_title='R_eff'
                      units="microns"
                      title="effective radius"                   
                   else if (qindex==12) then  ! Mean maximum dimension
                      work1=work1+Sum(lid_scatt_info(isc)%R*Nsize)*2.0
                      work2=work2+Sum(Nsize)
                      Quantity_grid(ix,iy,iz)=work1/work2
                      plot_title='<D> [microns]'
                      nc_title='Mean_Maximum_Dimension'
                      units="microns"
                      title="Mean Maximum Dimension"     
                   else if (qindex==13) then ! Visable extinction
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      !
                      work1=work1+Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      work2=1.0
                      plot_title='Extinction [1/m]'
                      nc_title='Extinction'
                      units="m^-1"
                      title="Extinction at "//trim(waves_val)//" nm"
                      !
                   else if (qindex==14) then ! alpha/beta ratio
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      work1=work1+Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      !
                      itheta=lid_scatt_info(isc)%n_angles
                      !
                      work2=work2+Sum(lid_scatt_info(isc)%p11(1,it,irh,:,itheta)*&
                           &lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      plot_title='Extinction/Backscatter [sr]'
                      nc_title='Extinction_Backscatter_ratio'
                      units='sr'
                      title='Extinction to Backscatter ratio at'//trim(waves_val)//' nm'
                      !
                   else if (qindex==15) then ! backscatter
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      work2=1.0
                      !
                      itheta=lid_scatt_info(isc)%n_angles
                      !
                      work1=work1+Sum(lid_scatt_info(isc)%p11(1,it,irh,:,itheta)*&
                           &lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      plot_title='Backscatter [1/m/sr]'
                      nc_title='Backscatter'
                      units='m^-1 sr^-1'
                      title='Backscatter at'//trim(waves_val)//' nm'
                   else if (qindex==16) then ! We want Rpeff
                      work1=work1+Sum((lid_scatt_info(isc)%Vol**2*Nsize))
                      work2=work2+Sum(lid_scatt_info(isc)%Ac*Nsize)
                      plot_title="R'\deff\u [microns]"
                      nc_title='Rpeff'
                      units="microns"
                      title="lidar-radar effective radius"   
                      !
                   else if (qindex==17) then ! We want No
                      work1=work1+Sum(Nsize)
                      work2=1.0 
                      plot_title="No [1/m\U3\D]"
                      nc_title='N_0'
                      units="cm^-3"
                      title="Np"
                   else if (qindex==18) then ! We want Ra
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      !
                      work2=work2+Sum(lid_scatt_info(isc)%Ac*lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)
                      work1=work1+Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)
                      plot_title="R\da\u [microns]"
                      nc_title='R_a'
                      units="microns"
                      title="R_a"
                   else
                      
                      write(error_str,'(A20,i4)') 'Unknown Option ',qindex
                      call write_error(error_str)
                      call write_error('Run program with not arguments to get help')
                      Call exit(1)
                      return
                   endif
                   !
                   DeAllocate(Nsize)
                   !
                endif
             endif
          enddo
          if (work2==0.0) then
             Quantity_grid(ix,iy,iz)=0.0
          else
             if (qindex==16) then
                Quantity_grid(ix,iy,iz)=(9.0/16.0/pi*work1/work2)**(0.25)
             else if (qindex==18) then
                Quantity_grid(ix,iy,iz)=sqrt((work2/work1)/pi)
             else
                Quantity_grid(ix,iy,iz)=work1/work2
             endif
          endif
       else if (qindex.ge.100) then
          !
          if ((qindex.ge.100).and.(qindex.le.105)) then   ! we want to calculate the reflectivity
             !
             if (qindex==100) then
                title='Unattenuated Reflectivity [94 GHz]'
                nc_title='Unattenuated_Reflectivity_94_GHz'
                units='mm^6 m^-3'
                plot_title='Unattenuated Reflectivity [mm\U6\D/m\U3\D] [94 GHz]'
             else if (qindex==101) then
                title='Radar attenuation [94 GHz]'
                nc_title='Radar_attenuation_94_GHz'
                units='m^-1'
                plot_title='Radar attenuation [1/m] [94 GHz]'
             else if (qindex==102) then
                title='Unattenuated Reflectivity [32 GHz]'
                nc_title='Unattenuated_Reflectivity_32_GHz'
                units='mm^6 m^-3'
                plot_title='Unattenuated Reflectivity [mm\U6\D/m\U3\D] [32 GHz]'
             else if (qindex==103) then
                title='Radar attenuation [32 GHz]'
                nc_title='Radar_attenuation_32_GHz'
                units='m^-1'
                plot_title='Radar attenuation [1/m] [32 GHz]'
             else if (qindex==104) then
                title='Unattenuated Reflectivity [5 GHz]'
                nc_title='Unattenuated_Reflectivity_5_GHz'
                units='mm^6 m^-3'
                plot_title='Unattenuated Reflectivity [mm\U6\D/m\U3\D] [5 GHz]'
             else if (qindex==105) then
                title='Radar attenuation [5 GHz]'
                nc_title='Radar_attenuation_5_GHz'
                units='m^-1'
                plot_title='Radar attenuation [1/m] [5 GHz]'
             endif
             !
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
                      if (qindex==100) then
                         !
                         work1=work1+Sum(Nsize*Ze_vec)
                         work2=1.0
                      else 
                         !
                         work1=work1+Sum(Nsize*ext_vec)*1.0e-4
                         work2=1.0
                         !
                      endif
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
          else
             !
             write(error_str,'(A20,i4)') 'Unknown Option ',qindex
             call write_info(error_str)
             call write_error('Run program with not arguments to get help')
             call exit(1)
             return
             !
          endif
          !
       endif
    enddo
    return
  end subroutine populate
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
    Character(len=160)                :: arg_str
    Character(len=190)               :: error_str
    Character(len=100)               :: error_str2
    !
    error_str='Usage: extract_quantity_3d inputfilename outputfilename x_start(km) x_finish(km)'//&
         &       ' y_start(km) y_finish(km) z1(km) z2(km) vert_res(km) qindex qindex2'
    error_str2='Usage: extract_quantity_3d inputfilename outputfilename x_start(km) x_finish(km)'//&
         &       ' y_start(km) y_finish(km) z1(km) z2(km) vert_res(km) qindex qindex2 wavelength'
    !
    status=0
    !
    nargs=iargc()  
    If (nargs.lt.11.or.nargs.gt.12) Then
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
       error_str2='Error in x starting value (km)'
       Read(arg_str,*,iostat=status,err=100) x_start
    Endif
    !
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str2='Error in x ending value (km)'
       Read(arg_str,*,iostat=status,err=100) x_finish
    Endif
    !
    !
    If (status == 0) Then
       Call getarg(5,arg_str)
       error_str2='Error in x starting value (km)'
       Read(arg_str,*,iostat=status,err=100) y_start
    Endif
    !
    If (status == 0) Then
       Call getarg(6,arg_str)
       error_str2='Error in x ending value (km)'
       Read(arg_str,*,iostat=status,err=100) y_finish
    Endif
    !
    If (status == 0) Then
       Call getarg(7,arg_str)
       error_str2='Error in instrument start altitude(km)'
       Read(arg_str,*,iostat=status,err=100) z1_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(8,arg_str)
       error_str2='Error in instrument stop altitude(km)'
       Read(arg_str,*,iostat=status,err=100) z2_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(9,arg_str)
       error_str2='Error in desired vert. resolution (km)'
       Read(arg_str,*,iostat=status,err=100) dz_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(10,arg_str)
       error_str2='Error in index of quantity to extract'
       Read(arg_str,*,iostat=status,err=100) qindex
    Endif
    !
    If (status == 0) Then
       Call getarg(11,arg_str)
       error_str2='Error in index of quantity to extract'
       Read(arg_str,*,iostat=status,err=100) qindex2
    Endif
    !
    if (nargs.eq.12) then 
       if ((qindex.eq.13).or.(qindex.eq.14).or.(qindex.eq.15)) then
          If (status == 0) Then
             Call getarg(12,arg_str)
             error_str2='Error in index of waves to extract'
             Read(arg_str,*,iostat=status,err=100) waves1
          Endif
          if ((waves1.lt.0.1).or.(waves1.gt.2)) then
             status=1
             error_str2=''
             error_str='Lidar wavelength should be in microns'
          endif
      else
         error_str2='The extracted quantity does not concern a lidar product. '//&
              &               '   There should be no index of wavelength'    
         status=1
       endif
    endif
    !
100 If (status.ne.0) Then 
       call Write_error(error_str)
       call Write_error(error_str2)
    Endif
    !
    If (status.eq.2) then
       call write_info( '-------------------------------')
       call write_info( 'Current Choices For qindex are:')
       call write_info( '-------------------------------')
       call write_info( '1  ===> Temperature [K]')
       call write_info( '2  ===> Pressure [mb]')
       call write_info( '3  ===> Density [molecules/cm^3]')
       call write_info( '4  ===> Xvol H2O')
       call write_info( '5  ===> Xvol O3')
       call write_info( '10 ===> mass [g/m^3]')
       call write_info( '11 ===> Reff [um]')
       call write_info( '12 ===> Mean maximum dimension [um]')
       call write_info( '13 ===> Extinction at 353 nm [1/m]')
       call write_info( '14 ===> Extinction/Backscatter Ratio at 353 nm [1/sr] ')
       call write_info( '15 ===> Backscatter at 353 nm [1/m/sr]')
       call write_info( "16 ===> R'eff [microns]")
       call write_info( "17 ===> No [#/cm^3]")
       call write_info( "18 ===> R_a [microns]")
       call write_info( '100 ==> Idealized Radar reflectivity [mm^6/m^3]')
       call write_info( '101 ==> Attenuation at 94 GHz [1/m].')
       call write_info( '102 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (32 GHz).')
       call write_info( '103 ==> Attenuation at 32 GHz [1/m].')
       call write_info( '104 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (5 GHz).')
       call write_info( '105 ==> Attenuation at 5 GHz [1/m].')
       call write_info( '---------------------------------------------------------------')
       call write_info( 'qindex2 ==> scatterer index (set to a specific scatterer index)')
       call write_info( 'set to 0 to use all scatterer types present')
       call write_info( '---------------------------------------------------------------')
       call write_info( ' In case of qindex=[13,14,15] an extra argument can be given   ')
       call write_info( ' This argument is the lidar wavelength. If non is given the    ')
       call write_info( ' default EarthCARE wavelength is given: 0.353                  ')
       call write_info( '---------------------------------------------------------------')
       call exit(1)
    endif
    !
  End Subroutine get_data
  
end Program extract_quantity_3d
!
!
 Subroutine write_results(filename,nx,ny,nz,x,y,z,Quant,title)
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: nx,ny,nz
    Real,intent(in)                :: x(nx,ny),y(nx,ny),z(nz)
    Real,intent(in)                :: Quant(nx,ny,nz)
    Character(len=*),intent(in)    :: title
    !
    Integer                        :: i,j,k
    !
    open(unit=3,file=filename,status='unknown',form='formatted',access='sequential')
    !
    write(3,*) title
    write(3,*) nx,ny,nz
    do i=1,nx
       do j=1,ny
          write(3,*) '=================='
          write(3,*) x(i,j),y(i,j)
          write(3,*) (z(k),k=1,nz)
          write(3,*) (Quant(i,j,k),k=1,nz)
       enddo
    enddo
    !
    close(3)
    !
    Return
  End Subroutine write_results
  !
!
  Subroutine write_results_ncdf(filename,nx,ny,Nz,x,y,z,Quant,nc_title,title,units,plot_title,status,error_str)
    !
    use typeSizes
    use netcdf
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: nx,ny
    Integer,intent(in)             :: Nz
    Real,intent(in)                :: x(nx),y(ny),z(Nz)
    Real,intent(in)                :: Quant(nx,ny,nz)
    Character(len=*),intent(in)    :: nc_title,title,units,plot_title
    !
    Integer                        :: i,j
    !
    Integer :: ncid, status
    CHARACTER(len=256)             ::error_str
    !
    integer :: X_dist, Y_dist, Altitude,ix,iy
    integer :: XDistId, YDistId, AlongTrackDistId, AcrossTrackDistId, HeightId, QuantId
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
    status = nf90_def_dim(ncid, "nx", nx, X_dist)    
    if (status /= 0) then 
       error_str='error in nf90_def_dim: nx'
       return
    endif
    status = nf90_def_dim(ncid, "ny", ny, Y_dist)    
    if (status /= 0) then 
       error_str='error in nf90_def_dim: ny'
       return
    endif
    status = nf90_def_dim(ncid, "nz", nZ,  Altitude)    
    if (status /= 0) then 
       error_str='error in nf90_def_dim: nz'
       return
    endif
    !
    ! Defining variables
    status = nf90_def_var(ncid, "x_scene", NF90_FLOAT, (/X_dist/), XDistId)
    if (status /= 0) then 
       error_str='error in nf90_def_var1'
       return
    endif
    status = nf90_def_var(ncid, "y_scene", NF90_FLOAT, (/Y_dist/), YDistId)
    if (status /= 0) then 
       error_str='error in nf90_def_var2'
       return
    endif
!    status = nf90_def_var(ncid, "along_track", NF90_FLOAT, (/X_dist,Y_dist/), alongtrackDistId)
!    if (status /= 0) then 
!       error_str='error in nf90_def_var2'
!       return
!    endif
!    status = nf90_def_var(ncid, "cross_track", NF90_FLOAT, (/X_dist,Y_dist/), acrosstrackDistId)
!    if (status /= 0) then 
!       error_str='error in nf90_def_var2'
!       return
!    endif
    status = nf90_def_var(ncid, "height", NF90_FLOAT, (/altitude/), HeightId)
    if (status /= 0) then 
       error_str='error in nf90_def_var3'
       return
    endif
    status = nf90_def_var(ncid, trim(nc_title), NF90_FLOAT, (/X_dist, Y_dist, Altitude/), QuantId)
    if (status /= 0) then 
       error_str='error in nf90_def_var4'
       return
    endif
    !
    ! Defining attributes
    !
    status=nf90_put_att(ncid,XDistId,"long_name","Distance along X-axis")
    if (status /= nf90_noerr) then 
       error_str='Error in nf90_atrdef1a'
       return
    endif
    status=nf90_put_att(ncid,XDistId,"units","km")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef1b'
       return
    endif
    status=nf90_put_att(ncid,XDistId,"plot_title"," Distance along X-axis [km]")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef1c'
       return
    endif
    !
    status=nf90_put_att(ncid,YDistId,"long_name","Distance along Y-axis")
    if (status /= nf90_noerr) then 
       error_str='Error in nf90_atrdef2a'
       return
    endif
    status=nf90_put_att(ncid,YDistId,"units","km")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef2b'
       return
    endif
    status=nf90_put_att(ncid,YDistId,"plot_title","Distance along Y-axis [km]")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef2c'
       return
    endif
    !
    status=nf90_put_att(ncid,HeightId,"long_name","Height")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef4a'
       return
    endif
    status=nf90_put_att(ncid,HeightId,"units","km")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef4b'
       return
    endif
    status=nf90_put_att(ncid,HeightId,"plot_title","Height [km]")
    if (status /= 0) then 
       error_str='Error in nf90_atrdef4c'
       return
    endif
    !
    status=nf90_put_att(ncid,QuantId,"long_name",title)
    if (status /= 0) then 
       error_str='Error in nf90_atrdef5a'
       return
    endif
    status=nf90_put_att(ncid,QuantId,"units",units)
    if (status /= 0) then 
       error_str='Error in nf90_atrdef5b'
       return
    endif
    status=nf90_put_att(ncid,QuantId,"plot_title",plot_title)
    if (status /= 0) then 
       error_str='Error in nf90_atrdef5c'
       return
    endif
    !
    status = nf90_enddef(ncid)
    if (status /= 0) then 
       error_str='Error in nf90_enddef'
       return
    endif
    !        
    ! Writing variables
    status = nf90_put_var(ncid, XDistId, x(1:nx))    
    if (status /= 0) then 
       error_str='Error in nf90_put_var1'
       return
    endif
    status = nf90_put_var(ncid, YDistId, y(1:ny))    
    if (status /= 0) then 
       error_str='Error in nf90_put_var2'
       return
    endif
    status = nf90_put_var(ncid, HeightId, z(1:nz))    
    if (status /= 0) then 
       error_str='Error in nf90_put_var3'
       return
    endif

    do ix=1,nx
       do iy=1,ny
          status = nf90_put_var(ncid, QuantId, Quant(ix,iy,:), &
               start=(/ ix,iy,1 /), &
               count=(/ 1,1,nz/)) 
          if (status /= 0) then 
             error_str='Error in nf90_put_var4'
             return
          endif
       enddo
    enddo


    ! Close the file
    status = nf90_close(ncid)
    if (status /= 0) then 
       error_str='Error in nf90_close'
       return
    endif
    !  
   Return
  End Subroutine write_results_ncdf
  !
  
