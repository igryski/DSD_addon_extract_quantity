!!
!! *PROGRAM* extract_quantity
!!
!! *USAGE*: extract_quantity inputfile.uff outfile.nc x1 y1 x2 y2
!!               z1(km) z2(km) vert_res(km)  
!!               qindex qindex2 (wavelength)
!! wavelength is only needed for qindex=13,14 or 15, 19 20 or 50
!!
!! @version ecsim 1.4
!!
!! *SRC_FILE*
!!
!! tools/product_tools/src/extract_quantity.f90
!!
!! *LAST CHANGES*
!! 
!! -Apr 15, 2015: I.S. Added DSD output option
!! -Feb 26, 2014: D.D. Added phase function output option
!! -Apr 03, 2013: D.D. Fixed error in generating reflectivity at 32 or 3 GHz (atteuation was output instead !!)
!! -Nov 24, 2011: D.D. Fixes added on Nov 22 were incomplete
!! -Nov 22, 2011: D.D. Added missing atmos_point deallocation statments..fixed running out of memory on large scenes.
!! -Nov 19, 2010: D.D. Fixed error in the calculation of g values
!! -Sep 23, 2010: D.D. Updated a few info and error messages to reflect new options 
!! -Jun 20, 2010: D.D. Fixed format error when wavelength is greater than 1000 nm
!! -Jun 15, 2010: D.D. Added SS_alb and g options
!! -Mar 10, 2008: D.D. Modified the header and Copyright info
!! -Feb 04, 2008, D.D. Cleaned up exit code reporting
!! -Jan  28, 2008: DD: Adjusted the along_track array (+ res/2)
!! -Dec  18, 2007, RV: Added 'header_' to UFF filename
!! -Oct  04, 2007, DD: Extended length of character strings in get_data and changed ESIM_HOME==>ECSIM_HOME
!! -Sept 24, 2007, DD: Changed prog messages and fixed wrap around problem when both ix1, and ix2 are less than one or greater than scene_nx etc..
!! -May  15, 2007, GJvZ: Added attributes for the NCDF output file
!! -May  16, 2007, GJvZ Added wavelength as 12th argument, value is passed to the NCDF output file
!!
!! *DESCRIPTION*
!!
!! Extract quantity is used to extract information from a UFF file along
!! a straight line between (x1,y1) and (x2,y2). The prgram will read the
!! UFF file as well as any relevent information from the scattering
!! libraries in order to buid-up the requested data.
!!
!! The following command line arguments are expected:
!! -inputfile.uff   : UFF file header to read from.
!! -outfile.nc      : netcdf file to write data to.
!! -x1,y1           : starting scene coordinates in km
!! -x2,y2           : ending scene coordinates in km
!! -z1              : Lower altitude bound in km.
!! -z2              : Upper altitude bound in km.
!! -vert_res        : Vertical output resolution in km.
!! -qindex          : Index of quantity to extract (see table below).
!! -qindex2         : Index of scattering type (0=> use all scattering types).
!! -wavelength      : Lidar wavelength (microns), only needed for qindex=13,14 or 15
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
!! -18 ===> R_a [microns]"
!! -19 ===> g
!! -20 ===> omega (SS_alb) 
!! -50 ===> Phase function 
!! -51 ===> DSD (Drop Size Distribution)
!! -100 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (94 GHz).
!! -101 ==> Attenuation at 94 GHz [1/m].
!! -102 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (32 GHz).
!! -103 ==> Attenuation at 32 GHz [1/m].
!! -104 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (5 GHz).
!! -105 ==> Attenuation at 5 GHz [1/m].
!! 
!! If the incorrect number of commandline arguments are given then
!! an error message is printed along with a help message.
!!
!! *OUTPUT*
!!
!! OUTFILENAME will contain the output information which may then
!! be plotted using plot_slice (except if phase function output is selected).
!! The output is in netcdf.
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
Program extract_quantity
  !
  Use scene_creator_types
  Use data_types
  Use physical_parameters
  Use read_uff
  Use Ray_SCATT
  Use write_messages
  !
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
  Real                           :: cnst,start,finish,width,waves1
  Real                           :: z1_ins,z2_ins,dz_ins,buffer,hor_res
  Integer                        :: qindex,qindex2
  !
  !----------------
  ! Main vaiables
  !----------------
  !
  Real                                      :: x_start,y_start,x_finish,y_finish
  Real                                      :: phi
  Real                                      :: h_res
  Character(len=256), Dimension (:), Pointer :: scatt_list_names
  Character(len=10), Dimension (:), Pointer :: gasses
  Integer                                   :: nscatt_types,n_gasses,scene_nx,scene_ny,scene_nz
  !
  Type(scatt_prop_master),Dimension(:),Allocatable   :: scatt_master_info ! structure containg scattering list info
  Type(scatterer_info),Dimension(:),Allocatable      :: rad_scatt_info    ! where the data is really stored
  Type(scatterer_info_pol),Dimension(:),Allocatable  :: lid_scatt_info    ! where the data is really stored
  Type(scatterer_info),Dimension(:),Allocatable  :: lid_scatt_info_nopol    ! where the data is really stored
  Type(scatterer_info),Dimension(:),Allocatable  :: lid_scatt_info_nopol_angs ! Interpolated to a standard angular grid
  Integer                                        :: nangles
  Real,dimension(:),Allocatable                  :: angles
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
  Real,Dimension(:),Allocatable              :: x,y,dist
  Real,Dimension(:,:),Allocatable            :: Quantity
  Real,Dimension(:,:,:),Allocatable            :: Quantity_ang
  !
  Real,Dimension(:,:),Allocatable            :: X_grid,Y_grid,weight_grid
  Real,Dimension(:,:,:),Allocatable          :: Quantity_grid 
  Real,Dimension(:,:,:,:),Allocatable          :: Quantity_grid_ang 
  !
  Real,Dimension(:),Allocatable              :: z_ins ! instrument resolution altitude vector (km)
  Integer                                    :: nz_ins
  !
  !-------------------------
  ! Misc working variables 
  !-------------------------
  !
  Integer                              :: status,i,ix,iy,iz,isc,j,ia
  Integer                              :: ix1,ix2,iy1,iy2,iix,iiy
  Integer                              :: ix1_tmp,iy1_tmp
  Integer                              :: ix2_tmp,iy2_tmp
  Type(scatterer_info)                 :: scatt_info_temp
  Integer                              :: ix_out_min,ix_out_max
  Integer                              :: iy_out_min,iy_out_max
  real                                 :: x1,y1,lat1,long1,z1
  real                                 :: x2,y2,lat2,long2,z2
  Integer                              :: same_res
  Integer                              :: ascii_or_bi
  Real                                 :: waves(2),freq
  Character(len=80)                    :: title, nc_title,units,plot_title
  real,dimension(:),allocatable        :: xg,yg
  integer                              :: nxg,nyg
  !
  character(len=256)                   :: ecsim_home, scatt_lib,error_str
  !
  !-------------------
  ! Weighting related
  ! variables
  !-------------------
  !
  real,dimension(:),pointer     :: xi,yi,ri
  Real                          :: xo,yo,dr,xo1,yo1
  Integer                       :: ni,nshots

  ! UFF header update
  Integer :: loc
  Character(len=256) :: char1, char2
  !
! ***************************
  character(len=256)    :: tempchar, tempname, tempname2
  integer :: charlen
! ***************************
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
     subroutine find_intercepts_2d(nx,ny,x,y,xo,yo,x1,y1,phi,xi,yi,ri,ni)
       !
       implicit none
       !
       integer,intent(in)    :: nx,ny
       real,intent(in)       :: x(nx),y(ny)
       real,intent(in)       :: xo,yo,phi,x1,y1
       !
       real,dimension(:),pointer  :: xi,yi,ri
       integer,intent(out)   :: ni
       !
     end subroutine find_intercepts_2d
  End Interface
  !
  nullify(xi)
  nullify(yi)
  nullify(ri)
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
     call write_error('Command line arguments not present or incomplete')
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
  !-----------------------------------
  ! Initalize some stuff regarding the 
  ! resolution and the path
  !-----------------------------------
  !
  if (x_start==x_finish) then
     phi=90.0
  else 
     phi=atand((y_finish-y_start)/(x_finish-x_start))
  endif
  !
  if (x_finish.lt.x_start) then
     phi=phi+180.0
  endif
  !
  !
  !---------------------------------------------------
  ! Open and read the the UFF file header information
  !---------------------------------------------------
  !
  call write_info('Reading the uff file :'//trim(adjustl(infilename)))
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
     call write_error('Must specifiy a positive dz_ins !')
     goto 200
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
     call write_warning( 'Requested resolution')
     call write_info( 'Coarser than the file resolution')
     call write_info( 'You will be able to average the ')
     call write_info( 'results later')
     call write_info( '*********************************')
     goto 200
  endif
  !
  if ((dz_ins.gt.(z2-z1)/scene_nz).and.(same_res.eq.0)) then 
     call write_info( '*********************************')
     call write_warning( 'Requested resolution')
     call write_info( 'Coarser than the file resolution')
     call write_info( 'You will be able to average the ')
     call write_info( 'results later')
     call write_info( '*********************************')
     goto 200
  endif
  !
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
     if (qindex.eq.19) then
        Allocate(lid_scatt_info_nopol(1:nscatt_types))
     else if (qindex.eq.50) then
        Allocate(lid_scatt_info_nopol(1:nscatt_types))
        Allocate(lid_scatt_info_nopol_angs(1:nscatt_types))
     else
        Allocate(lid_scatt_info(1:nscatt_types))
     endif
     !
     !
     if (waves1.le.0.0) then 
        waves1=0.353
     endif
     
     waves=waves1
     !
     if (qindex==50) then
        nangles=1801
        allocate(angles(1:nangles))
        do ia=1,nangles
           angles(ia)=(ia-1)*180.0/(nangles-1)
        enddo
     endif
     !
     Do i=1,nscatt_types
        !
        if (qindex.eq.19) then
           call find_scatt_at_waves(7,scatt_master_info(i),lid_scatt_info_nopol(i),&
                & 2,waves)
        else if (qindex.eq.50) then
           !
           call find_scatt_at_waves(7,scatt_master_info(i),lid_scatt_info_nopol(i),&
                & 2,waves)
           !
           call find_scatt_at_angle(lid_scatt_info_nopol_angs(i),lid_scatt_info_nopol(i),nangles)
           !
        else
           call find_scatt_at_waves_pol(7,scatt_master_info(i),lid_scatt_info(i),&
                & 2,waves)
        endif
        !
     Enddo
     !
  endif
  !
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
  ix1_tmp=Int((x_start-x1)/hor_res+1)
  ix2_tmp=Int((x_finish-x1)/hor_res+1)
  iy1_tmp=Int((y_start-y1)/hor_res+1)
  iy2_tmp=Int((y_finish-y1)/hor_res+1)
  !
  ix1=min(ix1_tmp,ix2_tmp)
  ix2=max(ix1_tmp,ix2_tmp)
  iy1=min(iy1_tmp,iy2_tmp)
  iy2=max(iy1_tmp,iy2_tmp)
  !
  Allocate(X_grid(ix1:ix2,iy1:iy2))
  Allocate(Y_grid(ix1:ix2,iy1:iy2))
  Allocate(Quantity_grid(ix1:ix2,iy1:iy2,1:nz_ins))
  Allocate(weight_grid(ix1:ix2,iy1:iy2))
  !
  X_grid=0.0
  Y_grid=0.0
  Quantity_grid=0.0
  weight_grid=0.0
  !
  if (qindex==50) then
     Allocate(Quantity_grid_ang(ix1:ix2,iy1:iy2,1:nz_ins,nangles))
     Quantity_grid_ang=0.0
  endif
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
  !
  !---------------------------------
  ! Read the uff_file and populate 
  ! Quantity_grid
  ! Process the continious
  ! part first then take any
  ! wrap around into account
  !---------------------------------
  !
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
           if (status.ne.0) goto 200
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
  !Call Close_uff(3)
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
           do iz=1,nz_ins
              Call deallocate_atmos_point(data_column(iz))
           enddo
           !
        endif
     enddo
  enddo
  !
  Call Close_uff(3)
  !
  !
  dr=hor_res
  !
  nshots=max(1,int(sqrt((x_finish-x_start)**2+(y_finish-y_start)**2)/dr))
  !
  !--------------------------
  ! Determine the weighting
  ! and apply 
  !--------------------------
  !
  xo=x_start
  yo=y_start
  !
  Allocate(x(1:nshots))
  Allocate(y(1:nshots))
  Allocate(dist(1:nshots))
  Allocate(Quantity(1:nz_ins,1:nshots))
  !
  if (qindex==50) then
     Allocate(Quantity_ang(1:nz_ins,1:nshots,1:nangles))
     Quantity_ang=0.0
  endif
  !
  allocate(xg(size(x_grid(:,iy1))))
  allocate(yg(size(y_grid(ix1,:))))
  xg=x_grid(:,iy1)
  yg=y_grid(ix1,:)
  nxg=size(xg)
  nyg=size(yg)
  !
  !
  do i=1,nshots
     !
     xo1=xo+dr*cosd(phi)
     yo1=yo+dr*sind(phi)
     !
     weight_grid=0.0
     !
     if (((int((xo-xg(1))/hor_res+0.5)==int((xo1-xg(1))/hor_res+0.5)).and.&
          (int((yo-yg(1))/hor_res+0.5)==int((yo1-yg(1))/hor_res+0.5))).or.&
          ((nxg==1).or.(nyg==1))) then
        !
        ni=1
        weight_grid(int((xo-x1)/Hor_res+1),int((yo-y1)/Hor_res+1))=1.0
        !
     else
        !
        call find_intercepts_2d(nxg,nyg,xg,yg,xo,yo,xo1,yo1,phi,xi,yi,ri,ni)
        !
        ri=ri/ri(ni)
        !
        do j=1,ni-1
           weight_grid(int((xi(j)-x1)/Hor_res+1),int((yi(j)-y1)/Hor_res+1))=ri(j+1)-ri(j)
        enddo
        !
     endif
     !
     x(i)=(xo+xo1)/2.
     y(i)=(yo+yo1)/2.
     !
     dist(i)=sqrt((xo-x_start)**2+(yo-y_start)**2)
     !
     xo=xo+dr*cosd(phi)
     yo=yo+dr*sind(phi)
     !
     if (qindex.ne.50) then
        do iz=1,nz_ins
           Quantity(iz,i)=Sum(Quantity_grid(:,:,iz)*weight_grid)
        enddo
     else
        do iz=1,nz_ins
           do ia=1,nangles
              Quantity_ang(iz,i,ia)=Sum(Quantity_grid_ang(:,:,iz,ia)*weight_grid)
           enddo
        enddo
     endif       
     !
  enddo
  !
  dist=dist+(dist(2)-dist(1))/2
  !
  
200 if (status.NE.0) then
     call write_error('The execution of extract_quantity stopped')
     call exit(1)
  endif

  !  call write_results(outfilename,size(x),size(x),nz_ins,x,y,dist,z_ins,Quantity,title)
  if (qindex.ne.50) then
     call write_results_ncdf(outfilename,size(x),size(x),nz_ins,x,y,dist,z_ins,Quantity,nc_title,title,units,plot_title)
  else
     call write_results_ncdf_ang(outfilename,size(x),size(x),nz_ins,nangles,x,y,dist,z_ins,angles,Quantity_ang,nc_title,title,units,plot_title)
  endif
 !
  !
  !
  call Write_info('******************FINISHED******************')
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
    Real,dimension(:),allocatable  :: work_ar                           
    Real,dimension(:),allocatable  :: Nsize
    Character(len=7)               :: waves_val
    !
    !
    !
    if (qindex==50) then 
       allocate(work_ar(1:nangles))
    endif
    ;
    write(unit=waves_val, fmt='(I7)') nint(waves1*1000)
    !
    Do iz=1,nz_ins
       !
       ! Populate the output arrays
       !
       work1=0.0
       work2=0.0
       !
       if (qindex==50) then 
          work_ar(:)=0.0
       endif
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
                If ((gasses(i)==' h2o').Or.(gasses(i).Eq.' H2O').or.(gasses(i).Eq.'H2O ')) Then
                   igass=i
                Endif
             Enddo
             If (igass==0) Then
                call Write_error('Water mixing ratio is not available !!')
                call Write_error('Check the UFF file !')
                call exit(1)
                return
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
                If ((gasses(i)==' o3 ').Or.(gasses(i).Eq.' O3 ').or.(gasses(i).Eq.'O3  ')) Then
                   igass=i
                Endif
             Enddo
             If (igass==0) Then
                call Write_error('O3 mixing ratio is not available !!')
                call Write_error('Check the UFF file !')
                call exit(1)
             Endif
             !
             work1=data_column(iz)%X_vol(igass)            ! O3 volumne mixing ratio
             !
          else
             write(error_str,'(a20,i5)') 'Unknown Option ',qindex
             call write_error(error_str)
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
          !
          Do isc=1,nscatt_types
             if ((qindex2.lt.1).or.(qindex2==isc)) then
                If ((Associated(data_column(iz)%size_dists(isc)%N_bin)).Or.&
                     (Associated(global_size_dists(isc,iz)%N_bin))) Then
                   !
                   if (qindex.eq.19) then
                      Allocate(Nsize(1:lid_scatt_info_nopol(isc)%n_sizes))
                   else if (qindex.eq.50) then
                      Allocate(Nsize(1:lid_scatt_info_nopol_angs(isc)%n_sizes))
                   else
                      Allocate(Nsize(1:lid_scatt_info(isc)%n_sizes))
                   endif
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
                   else if (qindex==20) then ! SS_alb
                      call find_irh_it_pol(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info(isc),irh,it)
                      !
                      work1=work1+Sum(lid_scatt_info(isc)%abs(1,it,irh,:)*Nsize)*1.0e-4 
                      work2=work2+Sum(lid_scatt_info(isc)%ext(1,it,irh,:)*Nsize)*1.0e-4 
                      plot_title='SS_alb'
                      nc_title='SS_alb'
                      units="-"
                      title="Single Scattering Albedo "//trim(waves_val)//" nm"
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
                   else if (qindex==19) then ! g
                      call find_irh_it(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info_nopol(isc),irh,it)
                      work1=work1+Sum(lid_scatt_info_nopol(isc)%g(1,it,irh,:)*Nsize*&
                           &(lid_scatt_info_nopol(isc)%ext(1,it,irh,:)-&
                           &lid_scatt_info_nopol(isc)%abs(1,it,irh,:)))
                      !
                      work2=work2+Sum((lid_scatt_info_nopol(isc)%ext(1,it,irh,:)-&
                           &lid_scatt_info_nopol(isc)%abs(1,it,irh,:))*Nsize)
                      !
                      plot_title='g'
                      nc_title='g'
                      units='-'
                      title='Ass. factor'//trim(waves_val)//' nm'
                      !
                   else if (qindex==50) then ! The Phase Function
                      nc_title="P11"       
                      call find_irh_it(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info_nopol_angs(isc),irh,it)
                      !
                      work1=Sum(lid_scatt_info_nopol_angs(isc)%ext(1,it,irh,:)*Nsize)
                      work2=work2+work1
                      !
                      do i=1,lid_scatt_info_nopol_angs(isc)%n_angles
                         work_ar(i)=work_ar(i)+work1*Sum(lid_scatt_info_nopol_angs(isc)%p11(1,it,irh,:,i)*Nsize)
                      enddo
                      !
                   else if (qindex==51) then ! The Drop Size Distribution
                      nc_title="DSD"       
                      call find_irh_it(data_column(iz)%T,&
                           & data_column(iz)%RH,lid_scatt_info_nopol_angs(isc),irh,it)
                      !
                      work1=Sum(Nsize)
                      work2=work2+work1
                      !
                      do i=1,lid_scatt_info_nopol_angs(isc)%n_angles
                         work_ar(i)=work_ar(i)+work1*Sum(Nsize)
                      enddo
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
                      write(error_str,'(a20,i5)')
                      call write_error( error_str)
                      call write_error('Run program with not arguments to get help')
                      call exit(1)
                      return
                   endif
                   !
                   DeAllocate(Nsize)
                   !
                endif
             endif
          enddo
          !
          if (work2==0.0) then
             if (qindex==50) then
                Quantity_grid_ang(ix,iy,iz,:)=0.0
             else 
                Quantity_grid(ix,iy,iz)=0.0
             endif
          else
             if (qindex==16) then
                Quantity_grid(ix,iy,iz)=(9.0/16.0/pi*work1/work2)**(0.25)
             else if (qindex==18) then
                Quantity_grid(ix,iy,iz)=sqrt((work2/work1)/pi)
             else if (qindex==50) then
                Quantity_grid_ang(ix,iy,iz,:)=work_ar(:)/work2
             else
                Quantity_grid(ix,iy,iz)=work1/work2
             endif
          endif
          !
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
                If ((qindex2.lt.1).or.(qindex2==isc)) then
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
!!$                      write(*,*) isc
!!$                      write(*,*) rad_scatt_info(isc)%P11(1,:,:,:,1)
!!$                      write(*,*) '-------------'
!!$                      write(*,*) Nsize/1.0e+4
                      !
                      !
                      ! D.D. April 03, 2013
                      !
                      if ((qindex==100).or.(qindex==102).or.(qindex==104)) then
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
                if (qindex.ne.20) then
                   Quantity_grid(ix,iy,iz)=work1/work2
                else
                   Quantity_grid(ix,iy,iz)=1.0-work1/work2
                endif
             endif
             !
          else
             !
             write(*,*) 'Unknown Option ',qindex
             write(*,*) 'Run program with not arguments to get help'
             call exit(1)
             return
             !
          endif
          !
       endif
    enddo
    !
    if (allocated(work_ar)) then
       DeAllocate(work_ar)
    endif
    !
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
    error_str='Usage: extract_quantity inputfilename outputfilename '//&
&       ' x1(km) y1(km) x2(km) y2(km) z1(km) z2(km) vert_res(km) qindex qindex2 wavelength'
    error_str2='Usage: extract_quantity inputfilename outputfilename '//&
&       ' x1(km) y1(km) x2(km) y2(km) z1(km) z2(km) vert_res(km) qindex qindex2 wavelength'
    !
    status=0
    start=0.0
    finish=0.0
    width=0.0
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
       error_str2='Error in x1 (km)'
       Read(arg_str,*,iostat=status,err=100) x_start
    Endif
    !
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str2='Error in x2 (km)'
       Read(arg_str,*,iostat=status,err=100) y_start
    Endif
    !
    !
    If (status == 0) Then
       Call getarg(5,arg_str)
       error_str2='Error in x2 (km)'
       Read(arg_str,*,iostat=status,err=100) x_finish
    Endif
    !
    If (status == 0) Then
       Call getarg(6,arg_str)
       error_str2='Error in y2 (km)'
       Read(arg_str,*,iostat=status,err=100) y_finish
    Endif
    !
!    If (status == 0) Then
!       Call getarg(7,arg_str)
!       error_str2='Error in width (km)'
!       Read(arg_str,*,iostat=status,err=100) width
!    Endif
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
       error_str2='Error in desired resolution (km)'
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
       if ((qindex.eq.13).or.(qindex.eq.14).or.(qindex.eq.15).or.(qindex.eq.19).or.(qindex.eq.20).or.(qindex.eq.50)) then
          If (status == 0) Then
             Call getarg(12,arg_str)
             error_str2='Error in index of waves to extract'
             Read(arg_str,*,iostat=status,err=100) waves1
          Endif
          if ((waves1.lt.0.1).or.(waves1.gt.400)) then
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
    If (start.Gt.finish) Then
       status=1
       error_str=' '
       error_str2='Start must be less than Finish !'
    Endif
    !
    If (width.Lt.0) Then
       status=1
       error_str=' '
       error_str2='Width must be greater than 0 ! '
    Endif
    !
100 If (status.ne.0) Then 
       call Write_error(error_str)
       call Write_error( error_str2)
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
       call write_info( '10 ===> Mass [g/m^3]')
       call write_info( '11 ===> Reff [um]')
       call write_info( '12 ===> Mean maximum dimension [um]')
       call write_info( '13 ===> Extinction at 353 nm [1/m]')
       call write_info( '14 ===> Extinction/Backscatter Ratio at 353 nm [1/sr] ')
       call write_info( '15 ===> Backscatter at 353 nm [1/m/sr]')
       call write_info( "16 ===> R'eff [microns]")
       call write_info( "17 ===> No [#/cm^3]")
       call write_info( "18 ===> R_a [microns]")
       call write_info( "19 ===> g")
       call write_info( "20 ===> omega (SS_alb)")
       call write_info( "50 ===> Phase Function")
       call write_info( '100 ==> Idealized Radar reflectivity [mm^6/m^3]')
       call write_info( '101 ==> Attenuation at 94 GHz [1/m].')
       call write_info( '102 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (32 GHz).')
       call write_info( '103 ==> Attenuation at 32 GHz [1/m].')
       call write_info( '104 ==> Idealized Radar reflectivity [mm$^6$/m$^3$] (5 GHz).')
       call write_info('105 ==> Attenuation at 5 GHz [1/m].')
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
  
end Program extract_quantity
!
!
 Subroutine write_results(filename,nmax,n,Nz,x,y,dist,z,Quant,title)
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: n,nmax
    Integer,intent(in)             :: Nz
    Real,intent(in)                :: x(nmax),y(nmax),dist(nmax),z(Nz)
    Real,intent(in)                :: Quant(Nz,nmax)
    Character(len=*),intent(in)    :: title
    !
    Integer                        :: i,j
    !
    open(unit=3,file=filename,status='unknown',form='formatted',access='sequential')
    !
    write(3,*) 1
    write(3,*) title
    write(3,*) 1,n
    write(3,*) Nz
    write(3,*) (x(i),i=1,n)
    write(3,*) (y(i),i=1,n)
    write(3,*) (dist(i),i=1,n)
    write(3,*) (z(i),i=1,Nz)
    !
    do i=1,n
       write(3,*) '=================='
       write(3,*) i
       write(3,*) (Quant(j,i),j=1,Nz)
    enddo
    !
    close(3)
    !
    Return
  End Subroutine write_results
  !
!
  Subroutine write_results_ncdf(filename,nmax,n,Nz,x,y,dist,z,Quant,nc_title,title,units,plot_title)
    !
    use typeSizes
    use netcdf
    Use write_messages
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: n,nmax
    Integer,intent(in)             :: Nz
    Real,intent(in)                :: x(nmax),y(nmax),dist(nmax),z(Nz)
    Real,intent(in)                :: Quant(Nz,nmax)
    Character(len=*),intent(in)    :: nc_title,title,units,plot_title
    !
    Integer                        :: i,j
    !
    Integer :: ncid, status
    !
    integer :: GroundDist, Altitude
    integer :: XDistId, YDistId, AlongTrackDistId, HeightId, QuantId
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
    status = nf90_def_dim(ncid, "nx", n, GroundDist)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: GroundDist'
       goto 400
    endif
    status = nf90_def_dim(ncid, "nz", nZ,  Altitude)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: Altitude'
       goto 400
    endif
    
    !
    ! Defining variables
    status = nf90_def_var(ncid, "x_scene", NF90_FLOAT, (/GroundDist/), XDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var1'
       goto 400
    endif
    status = nf90_def_var(ncid, "y_scene", NF90_FLOAT, (/GroundDist/), YDistId)
    if (status /= 0) then
       error_str= 'error in nf90_def_var2'
       goto 400
    endif
    status = nf90_def_var(ncid, "along_track", NF90_FLOAT, (/GroundDist/), AlongTrackDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var3'
       goto 400
    endif
    status = nf90_def_var(ncid, "height", NF90_FLOAT, (/Altitude/), HeightId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var4'
       goto 400
    endif
    status = nf90_def_var(ncid, trim(nc_title), NF90_FLOAT, (/Altitude, GroundDist/), QuantId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var5'
       goto 400
    endif
    !
    ! Defining attributes
    !
    status=nf90_put_att(ncid,XDistId,"long_name","Distance along X-scene")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"plot_title"," Distance along X-scene [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,YDistId,"long_name","Y-scene")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef2a'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2b'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"plot_title","Y-scene [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,AlongTrackDistId,"long_name","Along Track Distance")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef3a'
       goto 400
    endif
    status=nf90_put_att(ncid,AlongTrackDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef3b'
       goto 400
    endif
    status=nf90_put_att(ncid,AlongTrackDistId,"plot_title","Along Track Distance [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef3c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,HeightId,"long_name","Height")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef4a'
       goto 400
    endif
    status=nf90_put_att(ncid,HeightId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef4b'
       goto 400
    endif
    status=nf90_put_att(ncid,HeightId,"plot_title","Height [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef4c'
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
    status = nf90_put_var(ncid, XDistId, x(1:n))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var1'
       goto 400
    endif
    status = nf90_put_var(ncid, YDistId, y(1:n))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var2'
       goto 400
    endif
    status = nf90_put_var(ncid, AlongTrackDistId, dist(1:n))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var3'
       goto 400
    endif
    status = nf90_put_var(ncid, HeightId, z(1:nz))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var4'
       goto 400
    endif
    status = nf90_put_var(ncid, QuantId, Quant(1:nz,1:n))    
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
  Subroutine write_results_ncdf_ang(filename,nmax,n,Nz,na,x,y,dist,z,angles,Quant_ang,nc_title,title,units,plot_title)
    !
    use typeSizes
    use netcdf
    Use write_messages
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: n,nmax
    Integer,intent(in)             :: Nz,na
    Real,intent(in)                :: x(nmax),y(nmax),dist(nmax),z(Nz)
    Real,intent(in)                :: angles(na)
    Real,intent(in)                :: Quant_ang(Nz,nmax,na)
    Character(len=*),intent(in)    :: nc_title,title,units,plot_title
    !
    Integer                        :: i,j
    !
    Integer :: ncid, status
    !
    integer :: GroundDist, Altitude, Angle
    integer :: XDistId, YDistId, AngleId,AlongTrackDistId, HeightId, QuantId
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
    status = nf90_def_dim(ncid, "nx", n, GroundDist)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: GroundDist'
       goto 400
    endif
    status = nf90_def_dim(ncid, "nz", nZ,  Altitude)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: Altitude'
       goto 400
    endif
    status = nf90_def_dim(ncid, "na", na,  Angle)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: Angle'
       goto 400
    endif
    
    !
    ! Defining variables
    status = nf90_def_var(ncid, "x_scene", NF90_FLOAT, (/GroundDist/), XDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var1'
       goto 400
    endif
    status = nf90_def_var(ncid, "y_scene", NF90_FLOAT, (/GroundDist/), YDistId)
    if (status /= 0) then
       error_str= 'error in nf90_def_var2'
       goto 400
    endif
    status = nf90_def_var(ncid, "along_track", NF90_FLOAT, (/GroundDist/), AlongTrackDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var3'
       goto 400
    endif
    status = nf90_def_var(ncid, "angle", NF90_FLOAT, (/Angle/), AngleId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var4'
       goto 400
    endif
    status = nf90_def_var(ncid, "height", NF90_FLOAT, (/Altitude/), HeightId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var5'
       goto 400
    endif
    !
    !
    status = nf90_def_var(ncid, trim(nc_title), NF90_FLOAT, (/Altitude, GroundDist, Angle/), QuantId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var6'
       goto 400
    endif
    !
    ! Defining attributes
    !
    status=nf90_put_att(ncid,XDistId,"long_name","Distance along X-scene")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"plot_title"," Distance along X-scene [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,YDistId,"long_name","Y-scene")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef2a'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2b'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"plot_title","Y-scene [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef2c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,AlongTrackDistId,"long_name","Along Track Distance")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef3a'
       goto 400
    endif
    status=nf90_put_att(ncid,AlongTrackDistId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef3b'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,AngleId,"long_name","Angle")
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef3d'
       goto 400
    endif
    status=nf90_put_att(ncid,AngleId,"units","Deg")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef3e'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,AlongTrackDistId,"plot_title","Along Track Distance [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef3c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,HeightId,"long_name","Height")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef4a'
       goto 400
    endif
    status=nf90_put_att(ncid,HeightId,"units","km")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef4b'
       goto 400
    endif
    status=nf90_put_att(ncid,HeightId,"plot_title","Height [km]")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef4c'
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
    status = nf90_put_var(ncid, XDistId, x(1:n))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var1'
       goto 400
    endif
    status = nf90_put_var(ncid, YDistId, y(1:n))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var2'
       goto 400
    endif
    status = nf90_put_var(ncid, AlongTrackDistId, dist(1:n))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var3'
       goto 400
    endif
    status = nf90_put_var(ncid, HeightId, z(1:nz))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var4'
       goto 400
    endif
    status = nf90_put_var(ncid, AngleId, angles(1:na))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var6'
       goto 400
    endif
    status = nf90_put_var(ncid, QuantId, Quant_ang(1:nz,1:n,1:na))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var6'
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
    
  End Subroutine write_results_ncdf_ang
  !
  !
  subroutine find_intercepts_2d(nx,ny,x,y,xo,yo,x1,y1,phi,xi,yi,ri,ni)
    !
    implicit none
    !
    integer,intent(in)    :: nx,ny
    real,intent(in)       :: x(nx),y(ny)
    real,intent(in)       :: xo,yo,phi,x1,y1
    !
    real,dimension(:),pointer  :: xi,yi,ri
    integer,intent(out)   :: ni
    !
    integer               :: nix,niy
    real                  :: xinc(nx+2),xyinc(nx+2)
    real                  :: yinc(ny+2),yxinc(ny+2)
    !
    real                  :: xw,yw,rx,ry
    real                  :: res_y,res_x,r
    real                  :: sinphi,cosphi,sm
    !
    real                  :: maxx,minx,maxy,miny
    integer               :: ix,iy,i,ii,idel
    !
    res_x=x(2)-x(1)
    res_y=y(2)-y(1)
    !
    sm=1.0e-7
    !
    maxx=min(max(xo,x1),maxval(x)+res_x/2.)
    minx=max(min(xo,x1),minval(x)-res_x/2.)
    maxy=min(max(yo,y1),maxval(y)+res_y/2.)
    miny=max(min(yo,y1),minval(y)-res_y/2.)
    !
    !
    xw=xo
    yw=yo
    !
!!$    stop
!!$    if (xo.ge.x(1)) then
!!$       ix=max(1,int((xo-x(1))/res_x+0.5)+1)
!!$    else
!!$       ix=min(nx,int((xo-x(1))/res_x-0.5)+1)
!!$    endif
!!$    !
!!$    if (yo.ge.y(1)) then
!!$       iy=max(1,int((yo-y(1))/res_y+0.5)+1)
!!$    else
!!$       iy=min(ny,int((yo-y(1))/res_y-0.5)+1)
!!$    endif
    !
    ix=max(1,int((xo-x(1))/res_x+0.5)+1)
    iy=max(1,int((yo-y(1))/res_y+0.5)+1)
    !
    sinphi=sind(phi)
    cosphi=cosd(phi)
    !
    xw=xo
    yw=yo
    !
    if (associated(xi)) then
       deallocate(xi)
       nullify(xi)
    endif
    if (associated(yi)) then
       deallocate(yi)
       nullify(yi)
    endif
    if (associated(ri)) then
       deallocate(ri)
       nullify(ri)
    endif
    !
    if (cosphi.ne.0.0) then 
       !
       ! Find the x crossing points
       !
       i=1
       do while (&
            & (xw.le.maxx).and.(xw.ge.minx).and.&
            & (yw.le.maxy).and.(yw.ge.miny))
          !
          xinc(i)=xw
          xyinc(i)=yw
          !
          if (cosphi.ge.0) then ! increasing x
             if (ix+i-1.gt.nx) then
                xw=maxx+1
             else
                xw=x(ix+i-1)+res_x/2.0
             endif
          else
             if (ix-i+1.lt.1) then
                xw=minx-1
             else
                xw=x(ix-i+1)-res_x/2.0 ! decreasing x
             endif
          endif
          !
          r=abs((xw-xo)/cosphi)
          yw=abs(r*sinphi)+yo
          !
          i=i+1
       enddo
       nix=i-1
       !
    else
       nix=1
       xinc(1)=xo
       xyinc(1)=yo
    endif
    !
    xw=xo
    yw=yo
    !
    if (sinphi.ne.0.0) then 
       !
       ! Find the y crossing points
       !
       i=1
       do while (&
            & (xw.le.maxx).and.(xw.ge.minx).and.&
            & (yw.le.maxy).and.(yw.ge.miny))
          !
          yxinc(i)=xw
          yinc(i)=yw
          !
          if (sinphi.gt.0.0) then ! Increasing y
             if (iy+i-1.gt.ny) then
                yw=maxy+1
             else
                yw=y(iy+i-1)+res_y/2.0
             endif
          else
             if (iy-i+1.lt.1) then
                yw=miny-1
             else
                yw=y(iy-i+1)-res_y/2.0 ! Decreasing y
             endif
          endif
          !
          r=abs((yw-yo)/sinphi)
          xw=abs(r)*cosphi+xo
          !
          i=i+1
       enddo
       niy=i-1
    else
       niy=1
       yxinc(1)=xo
       yinc(1)=yo
    endif
    !
    !----------------------------
    ! Now combine the results
    !----------------------------
    !
    ni=nix+niy-1
    !
!!$    write(*,*) xo,x1
!!$    write(*,*) yo,y1
!!$    !
!!$    write(*,*) minx,maxx
!!$    write(*,*) miny,maxy
!!$    !
!!$    write(*,*) xinc
!!$    write(*,*) xyinc
!!$    !
!!$    write(*,*)
!!$    !
!!$    write(*,*) yinc
!!$    write(*,*) yxinc
    !
    !
    allocate(xi(1:max(ni,2)))
    allocate(yi(1:max(ni,2)))
    allocate(ri(1:max(ni,2)))
    !
    ix=2
    iy=2
    !
    xi(1)=xo
    yi(1)=yo
    ri(1)=0.0
    !
    ii=2
    do i=2,ni-1
       !
       if (ix.gt.nix) then
          rx=9e+20
       else
          rx=(xinc(ix)-xo)**2+(xyinc(ix)-yo)**2
       endif
       if (iy.gt.niy) then
          ry=9e+20
       else
          ry=(yxinc(iy)-xo)**2+(yinc(iy)-yo)**2
       endif
       !
       if (rx.le.ry) then
          xi(ii)=xinc(ix)
          yi(ii)=xyinc(ix)
          ri(ii)=sqrt(rx)
          ix=ix+1
       else 
          xi(ii)=yxinc(iy)
          yi(ii)=yinc(iy)
          ri(ii)=sqrt(ry)
          iy=iy+1
       endif
       !
       if ((abs(xi(ii)-xi(ii-1)).gt.sm).or.&
            & (abs(yi(ii)-yi(ii-1)).gt.sm)) then
          ii=ii+1
       endif
    enddo
    !
    if (ni.gt.1) then
       xi(ii)=x1
       yi(ii)=y1
       ri(ii)=sqrt((x1-xo)**2+(y1-yo)**2)
       !
    else
       !
       ni=2
       xi(2)=x1
       yi(2)=y1
       ri(2)=sqrt((x1-xo)**2+(y1-yo)**2)
       !
    endif
!!$
!!$    rx=(xinc(nix)-xo)**2+(xyinc(nix)-yo)**2
!!$    ry=(yxinc(niy)-xo)**2+(yinc(niy)-yo)**2
!!$    !
!!$    if (ni.gt.1) then
!!$       if (rx.ge.ry) then
!!$          xi(ii)=xinc(nix)
!!$          yi(ii)=xyinc(nix)
!!$          ri(ii)=sqrt(rx)
!!$       else 
!!$          xi(ii)=yxinc(niy)
!!$          yi(ii)=yinc(niy)
!!$          ri(ii)=sqrt(ry)
!!$       endif
!!$       ni=ii
!!$    else
!!$       xi(1)=xo
!!$       yi(1)=yo
!!$       ri(1)=1.0
!!$    endif
!!$    !

    !
    return
    !
  end subroutine find_intercepts_2d
