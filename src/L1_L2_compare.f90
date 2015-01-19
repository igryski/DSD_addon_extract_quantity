!!
!! *PROGRAM* L1_L2_compare
!!
!! USAGE: L1_L2_compare input1.nc title1 input2.nc title2 output.nc
!!
!! @version ecsim 1.2.3
!!
!! *SRC_FILE*
!!
!! product_tools/tools/src/L1_L2_compare.f90
!!
!! *LAST CHANGES*
!! 
!! -Nov 24, 2008: D.D. Various ncdf attribute character strings were too short ! This created propblems with the output of merger_lid_rad_msi
!!                     These attribute characters now have a length of 256. ALSO. Fixed Oct 17 fix. It was only applied to the case of
!!                     horizontal data files.
!! -Oct 17, 2008: D.D. Fixed it so that relative diff =0 when abs diff=0 even when data2 = 0.
!! -Mar 12, 2008: D.D. Fixed problem when croaa_track was not present in the input ncdf files
!! -Mar 10, 2008: D.D. Modified the header and Copyright info
!! -Jan 14, 2008: D.D. Fixed it so it works well with 2d horizontal nc files
!! -Oct 04, 2007: D.D. Extended length of character strings in get_data and changed ESIM_HOME==>ECSIM_HOME
!!
!! * DESCRIPTION*
!!
!! L1_L2_compare is used to compare the quantities in two different 2-D ncdf files. The program will read the
!! data from the netcdf files and return the diffences (absolute and relative of file1-file2) in the output netcdf file..
!! It reads in data from the L1_L2_rebin tool. This means that there will be a predefined slab in the 
!! data stream of the input file and this argument is not needed here.
!!
!! The following command line arguments are expected:
!! 
!!
!! -input1.nc : netcdf file 1 to read from.
!! -title1    : data stream to read in file 1 
!! -input2.nc : netcdf file 2 to read from.
!! -title2    : data stream to read in file 2  
!! -output.nc : netcdf file to write the data to
!! 
!! If the incorrect number of commandline arguments are given then
!! an error message is printed along with a help message.
!!
!! *OUTPUT*
!!
!! output.nc will contain the output information which may then
!! be plotted using plot_slice. The output is in netcdf.
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
Program L1_L2_compare
  ! 
  Use write_messages
  Use netcdf
  Use ncdf_utilities
  !
  Implicit None
  !
  !
  !------------------
  ! Input variables
  !------------------
  !
  Character(len=256)             :: outfilename
  Character(len=256)             :: infilename1,infilename2
  Real                           :: cnst,start,finish,width,waves1
  Real                           :: z1_ins,z2_ins,dz_ins,buffer,hor_res
  Integer                        :: qindex,qindex2,ifile
  !
  !----------------
  ! Main vaiables
  !----------------
  !
  Real                                      :: phi
  Real                                      :: x_max_diff,y_max_diff,z_max_diff
  Integer                                   :: scene_nx1,scene_nz1,scene_nx2,scene_nz2
  !
  Integer                                   :: nx1,ny1,nz1,nx2,ny2,nz2
  Real,Dimension(:),Allocatable             :: along_track1,cross_track1
  Real,Dimension(:),Allocatable             :: z1,z2
  Real,Dimension(:),Allocatable             :: along_track2,cross_track2
  Real,Dimension(:),Allocatable             :: x_scene1,y_scene1
  Real,Dimension(:),Allocatable             :: x_scene2,y_scene2
  !

  Real,Dimension(:,:),Allocatable            :: data1,data2,data_diff_abs,data_diff_rel
  !
  !--------------------
  ! Attribute strings
  !--------------------
  !
  character(len=256) :: title_x,units_x,title_plot_x
  character(len=256) :: title_y,units_y,title_plot_y
  character(len=256) :: title_z,units_z,title_plot_z
  character(len=256) :: title_xscene,title_yscene
  character(len=256) :: units_xscene,units_yscene
  character(len=256) :: title_plot_xscene,title_plot_yscene
  !
  !-------------------------
  ! Misc working variables 
  !-------------------------
  !
  Integer                              :: status,i,ix,iy,iz,isc,j,hor
  Integer                              :: same_res
  Character(len=256)                    :: title, nc_title, units,plot_title
  Character(len=256)                    :: title_1,title_plot_1,units_1,title1,out_title
  Character(len=256)                    :: title_2,title_plot_2,units_2,title2,out_plot_title
  Character(len=256)                    :: title_x1,units_x1,title_plot_x1
  Character(len=256)                    :: title_y1,units_y1,title_plot_y1
  Character(len=256)                    :: title_z1,units_z1,title_plot_z1
  Character(len=256)                    :: title_x2,units_x2,title_plot_x2
  Character(len=256)                    :: title_y2,units_y2,title_plot_y2
  Character(len=256)                    :: title_z2,units_z2,title_plot_z2  
  character(len=256)                   :: ecsim_home, scatt_lib,error_str
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
     call write_error('Command line arguments not present or incomplete')
     goto 200
  Endif
  !
  !----------------------------------
  ! Add ECSIM_HOME to the input and 
  ! output files
  !----------------------------------
  !
  infilename1=trim(adjustl(ecsim_home))//infilename1
  infilename2=trim(adjustl(ecsim_home))//infilename2
  outfilename=trim(adjustl(ecsim_home))//outfilename
  !
  !
  !---------------------------------------------------
  ! Open and read the the NETCDF file header information
  !---------------------------------------------------
  !
  hor=0
  IFILE=1
  !
  call read_data_ncdf(11,infilename1,status,title1,title_plot_1,units_1,data1,&
       & along_track1,cross_track1,x_scene1,y_scene1,z1,hor)
  if (status /= 0) then
     call write_error('Error in read_data_ncdf file 1')
     call exit(status)
  endif
  !
  title_x1=title_x
  units_x1=units_x
  title_plot_x1=title_plot_x
  title_y1=title_y
  units_y1=units_y
  title_plot_y1=title_plot_y
  title_z1=title_z
  units_z1=units_z
  title_plot_z1=title_plot_z
  !
  IFILE=2
  call read_data_ncdf(12,infilename2,status,title2,title_plot_2,units_2,data2,&
       & along_track2,cross_track2,x_scene2,y_scene2,z2,hor)
  if (status /= 0) then
     call write_error('Error in read_data_ncdf file 2')
     call exit(status)
  endif
  !
  !
  title_x2=title_x
  units_x2=units_x
  title_plot_x2=title_plot_x
  title_y2=title_y
  units_y2=units_y
  title_plot_y2=title_plot_y
  title_z2=title_z
  units_z2=units_z
  title_plot_z2=title_plot_z
  !
  if (hor==0) then
     !
     ! Check if the dimensions and x and z -values are the same
     !
     if ((nx1.ne.nx2).or.(nz1.ne.nz2)) then
        status=1
        call write_error('The along_track or z dimensions of the two files are not the same')
        call write_error('Use l1_l2_rebin to get the same formats') 
        goto 200
     endif
     !
     ! Check for the maximum difference between the 2 along_track and height grids
     !
     !
     z_max_diff=maxval(abs(z2-z1))
     x_max_diff=maxval(abs(along_track2-along_track1))
     !
     if (x_max_diff.gt.1.e-3) then 
        status=1
        call write_error('The two along_track grids are not the same')
        call write_error('Use l1_l2_rebin to get the same formats') 
        goto 200
     endif
     if (z_max_diff.gt.1.e-3) then 
        status=1
        call write_error('The two z-grids are not the same')
        call write_error('Use l1_l2_rebin to get the same formats') 
        goto 200
     endif
     allocate(data_diff_abs(1:nx1,1:nz1),data_diff_rel(1:nx1,1:nz1))
     
     data_diff_abs=data1-data2
     data_diff_rel(:,:)=-9999.0
     
     where(data2.ne.0.0)
        data_diff_rel=(data1-data2)/data2*100.0
     end where
     !
     where(data_diff_abs.eq.0.0)
        data_diff_rel=0.0
     end where
     !
     call write_results_ncdf(trim(adjustl(outfilename)))
     !
  else
     !
     ! Check if the dimensions and x and y -values are the same
     !
     if ((nx1.ne.nx2).or.(ny1.ne.ny2)) then
        status=1
        call write_error('The along_track or cross_track dimensions of the two files are not the same')
        call write_error('Use l1_l2_rebin to get the same formats') 
        goto 200
     endif
     !
     ! Check for the maximum difference between the 2 along_track and height grids
     !
     !
     y_max_diff=maxval(abs(cross_track2-cross_track1))
     x_max_diff=maxval(abs(along_track2-along_track1))
     !
     if (x_max_diff.gt.1.e-3) then 
        status=1
        call write_error('The two along_track grids are not the same')
        call write_error('Use l1_l2_rebin to get the same formats') 
        goto 200
     endif
     if (y_max_diff.gt.1.e-3) then 
        status=1
        call write_error('The two cross_track grids are not the same')
        call write_error('Use l1_l2_rebin to get the same formats') 
        goto 200
     endif
     allocate(data_diff_abs(1:nx1,1:ny1),data_diff_rel(1:nx1,1:ny1))
     
     data_diff_abs=data1-data2
     data_diff_rel(:,:)=-9999.0
     
     where(data2.ne.0.0)
        data_diff_rel=(data1-data2)/data2*100.0
     end where
     !
     ! D.D. Oct 17, 2008
     !
     where(data_diff_abs.eq.0.0)
        data_diff_rel=0.0
     end where
     !
     nz1=1
     nz2=1
     !
     call write_results_ncdf(outfilename)
     !
  end if
  !
  !
  !
  call Write_info('******************FINISHED******************')
  call exit(0)
200 call Write_error('An error in L1_L2_compare occured')
  call exit(status)
  ! 
Contains
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
    error_str='Usage:L1_L2_compare INPUTFILENAME1 title1 INPUTFILENAME2 title2  OUTPUTFILENAME '

    error_str2='Usage:L1_L2_compare INPUTFILENAME1 title1 INPUTFILENAME2 title2  OUTPUTFILENAME '
    !
    status=0
    !
    nargs=iargc() 
 
    If (nargs.ne.5) Then
       status=1
       error_str2='Wrong number of arguments'
    Endif
    !
    If (status == 0) Then
       Call getarg(1,arg_str)
       error_str2='Error in infilename1 specification'
       Read(arg_str,'(a160)',iostat=status,err=100) infilename1
    Endif
    If (status == 0) Then
       Call getarg(2,arg_str)
       error_str2='Error in title1 specification'
       Read(arg_str,'(a80)',iostat=status,err=100) title1
    Endif
    !
    If (status == 0) Then
       Call getarg(3,arg_str)
       error_str2='Error in infilename2 specification'
       Read(arg_str,'(a160)',iostat=status,err=100) infilename2
    Endif
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str2='Error in title2 specification'
       Read(arg_str,'(a80)',iostat=status,err=100) title2
    Endif
    If (status == 0) Then
       Call getarg(5,arg_str)
       error_str2='Error in outfilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) outfilename
    Endif
    !If (status == 0) Then
    !   Call getarg(5,arg_str)
    !   error_str2='Error in out title specification'
    !   Read(arg_str,'(a160)',iostat=status,err=100) out_title
    !Endif
    !
    !If (status == 0) Then
    !   Call getarg(7,arg_str)
    !   error_str2='Error in out plot title specification'
    !   Read(arg_str,'(a160)',iostat=status,err=100) out_plot_title
    !Endif
    !
    !
100 If (status.ne.0) Then 
       call Write_error(error_str)
       call Write_error( error_str2)
    Endif
    
    !
  End Subroutine get_data
  !
  Subroutine read_data_ncdf(id,inputfile,exit_status,title_q,plot_title,units,data,&
       &along_track,cross_track,x_scene,y_scene,z,hor)
    !
    use typeSizes
    use netcdf
    !
    integer        :: ncid,ndims,nvars,id
    Character(len=*),intent(in) :: inputfile,title_q
    integer, intent(out) :: exit_status
    Character(len=*),intent(out)  :: plot_title,units
    !
    ! NEW variables
    Integer :: status, io_mode, iDim ,iVar, length 
    character(len=nf90_max_name) :: name, err_msg
    integer :: found,iwork1,iwork2,iwork3(10),nx,ny,nz
    real,Dimension(:,:),Allocatable :: data
    real,Dimension(:),Allocatable   :: along_track
    real,Dimension(:),Allocatable   :: cross_track
    real,Dimension(:),Allocatable   :: x_scene
    real,Dimension(:),Allocatable   :: y_scene
    real,Dimension(:),Allocatable   :: z
    integer, allocatable, dimension(:) :: DimLengths
    character(len=80), allocatable, dimension(:) :: Dimnames
    Integer                            :: hor
    !
    exit_status = 0
    
    !
    io_mode = 0  ! default: nf90_nowrite
    !
    status = nf90_open(inputfile, io_mode, ncid)
    !
    err_msg = 'Error in nf90_open '//inputfile
    if (status/=0) then
       call write_error(err_msg)
       exit_status=1 
       goto 999
    endif
    !
    ! Get info on
    ! - nDims
    ! - nVars
    status  = nf90_inquire(ncid, nDims, nVars)   
    if (status /= 0) then
       call write_error('Error in nf90_inquire in L1_L2_compare read_data_ncdf')
       call exit(1)
    endif
    !
    ! Create arrays that contain Dimension Info   
    allocate(DimNames(1:nDims))
    allocate(DimLengths(1:nDims))
    !
    nz=-99
    !
    ! Get info on the dimensions:
    Do iDim = 1,nDims
       Call my_nf90_inquire_dimension(ncid, idim, status, name=name, len=length) 
       if (status/=0) then
          call write_error(err_msg)
          exit_status=1 
          goto 999
       endif
       dimNames(iDim)   = name
       dimLengths(iDim) = length
       Select Case (name)
       Case ('nx') ;  nx = length
       Case ('ny') ;  ny = length
       Case ('nz') ;  nz = length
       End Select
    End Do
    !
    if (nz.lt.0) then
       hor=1
       nz=1
    else
       hor=0
    endif
    !
    found=0
    !
    title_y=' '
    !
    do iVar=1,nVars
       !
       status = nf90_inquire_variable(ncid, iVar, title,iwork1,iwork2,iwork3)
       if (status/=0) goto 999
       status =nf90_get_att(ncid, iVar, "plot_title",plot_title)
       if (status/=0) then
          plot_title='--'
          status=0
       endif
       status =nf90_get_att(ncid, iVar, "units",units)
       if (status/=0) then
          units='--'
          status=0
       endif
       if (title.eq.'along_track') then
          allocate(along_track(1:nx))
          title_x=trim(title)
          units_x=trim(units)
          title_plot_x=trim(plot_title)
          status = nf90_get_var(ncid, iVar, along_track)
          !
       elseif (title.eq.'cross_track') then
          allocate(cross_track(1:ny))
          title_y=trim(title)
          units_y=trim(units)
          title_plot_y=trim(plot_title)
          status = nf90_get_var(ncid, iVar, cross_track)
          !
       elseif (title.eq.'x_scene') then
          allocate(x_scene(1:nx))
          title_xscene=trim(title)
          units_xscene=trim(units)
          title_plot_xscene=trim(plot_title)
          status = nf90_get_var(ncid, iVar, x_scene)
          !
       elseif (title.eq.'y_scene') then
          allocate(y_scene(1:nx))
          title_yscene=trim(title)
          units_yscene=trim(units)
          title_plot_yscene=trim(plot_title)
          status = nf90_get_var(ncid, iVar, y_scene)
          !
       elseif (title.eq.'height') then
          allocate(z(1:nz))
          title_z=trim(title)
          units_z=trim(units)
          title_plot_z=trim(plot_title)
          status = nf90_get_var(ncid, iVar, z)
          !
       elseif (trim(title).eq.trim(title_q)) then
          if (iwork2.ne.2) then
             call write_error('The array to compare is not 2 dimensional')
             exit_status=1
             goto 999
          endif
          !
          if (hor==1) then
             allocate(data(1:nx,1:ny))
          else
             allocate(data(1:nx,1:nz))
          endif
          !
          units=trim(units)
          title=trim(plot_title)
          status = nf90_get_var(ncid, iVar, data)
          found=1
       endif
       
       if (status/=0) then
          call write_error(err_msg)
          exit_status=1 
          goto 999
       endif
    enddo
    !
    if (len(trim(title_y))==0) then ! cross_track has not be defined in the ncdf file
       allocate(cross_track(1:ny))
       title_y="cross_track"
       units_y="km"
       title_plot_y="Across Track Distance [km]"
       cross_track=0.0
    endif
    !
    if (found /= 1) then 
       err_msg = 'Variable not found: '// trim(title_q)
       !
       call write_error(err_msg)
       call write_info('The 2-D variables in this file are:')
       !
       do iVar=1,nVars
          status = nf90_inquire_variable(ncid, iVar, title,iwork1,iwork2)
          if (iwork2==2) then
             call write_info(title)
          endif
       enddo
       !
       exit_status=1 
       goto 999
    endif
    !
    if (IFILE.EQ.1) THEN
    !   allocate(x1(1:nx),z1(1:nz),data1(1:nx,1:nz))
       nx1=nx
       nz1=nz
       ny1=ny
    elseif (IFILE.EQ.2) THEN
    !   allocate(x2(1:nx),z2(1:nz),data2(1:nx,1:nz))
       nx2=nx
       nz2=nz
       ny2=ny
    else
       call write_error('Incorrect file number')
       exit_status=1
       goto 999
    endif
    !
    !
999 return
  end Subroutine read_data_ncdf
  !

  Subroutine write_results_ncdf(filename)
    !
    use typeSizes
    use netcdf
    Use write_messages
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename
    !
    Integer                        :: i,j
    !
    Integer :: ncid, status
    !
    integer :: nx,ny,nz
    integer :: XDistId, YDistId,HeightID, QuantId1,QuantId2
    integer :: xsceneid,ysceneid
    Character(len=190)               :: error_str,title_q1,title_q2
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
    status = nf90_def_dim(ncid, "nx", nx1, nx)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: nx'
       goto 400
    endif
    !
    status = nf90_def_dim(ncid, "ny", ny1, ny)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: ny'
       goto 400
    endif
    !
    status = nf90_def_dim(ncid, "nz", nz1,  nz)    
    if (status /= 0) then 
       error_str= 'error in nf90_def_dim: nz'
       goto 400
    endif    
    !
    !
    ! Defining variables
    status = nf90_def_var(ncid, title_x1, NF90_FLOAT, (/nx/), XDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var1'
       goto 400
    endif
    !
    status = nf90_def_var(ncid, title_y1, NF90_FLOAT, (/ny/), YDistId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var2'
       goto 400
    endif
    !
    if (hor==0) then
       status = nf90_def_var(ncid, title_z1, NF90_FLOAT, (/nz/), HeightId)
       if (status /= 0) then 
          error_str= 'error in nf90_def_var3'
          goto 400
       endif
    endif
    !
    title_q1='abs_diff_'//trim(adjustl(title1))
    title_q2='rel_diff_'//trim(adjustl(title1))
    !
    !
    status = nf90_def_var(ncid, title_xscene, NF90_FLOAT, (/nx/), xsceneId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var3'
       goto 400
    endif
    !
    status = nf90_def_var(ncid, title_yscene, NF90_FLOAT, (/nx/), ysceneId)
    if (status /= 0) then 
       error_str= 'error in nf90_def_var4'
       goto 400
    endif
    !
    if (hor==0) then
       !
       status = nf90_def_var(ncid, title_q1, NF90_FLOAT, (/nx, nz/), QuantId1)
       if (status /= 0) then 
          error_str= 'error in nf90_def_var5'
          goto 400
       endif
       !
       status = nf90_def_var(ncid, title_q2, NF90_FLOAT, (/nx, nz/), QuantId2)
       if (status /= 0) then 
          error_str= 'error in nf90_def_var6'
          goto 400
       endif
       !
    else
       !
       status = nf90_def_var(ncid, title_q1, NF90_FLOAT, (/nx, ny/), QuantId1)
       if (status /= 0) then 
          error_str= 'error in nf90_def_var5'
          goto 400
       endif
       !
       status = nf90_def_var(ncid, title_q2, NF90_FLOAT, (/nx, ny/), QuantId2)
       if (status /= 0) then 
          error_str= 'error in nf90_def_var6'
          goto 400
       endif
       !
    endif
    !
    ! Defining attributes
    !
    status=nf90_put_att(ncid,XDistId,"long_name",trim(title_x1))
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"units",trim(units_x1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,XDistId,"plot_title",trim(title_plot_x1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,YDistId,"long_name",trim(title_y1))
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"units",trim(units_y1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,YDistId,"plot_title",trim(title_plot_y1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,XsceneId,"long_name",trim(title_xscene))
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,XsceneId,"units",trim(units_xscene))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,XsceneId,"plot_title",trim(title_plot_xscene))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    status=nf90_put_att(ncid,YsceneId,"long_name",trim(title_yscene))
    if (status /= nf90_noerr) then 
       error_str= 'Error in nf90_atrdef1a'
       goto 400
    endif
    status=nf90_put_att(ncid,YsceneId,"units",trim(units_yscene))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1b'
       goto 400
    endif
    status=nf90_put_att(ncid,YsceneId,"plot_title",trim(title_plot_yscene))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef1c'
       goto 400
    endif
    !
    if (hor==0) then
       status=nf90_put_att(ncid,HeightId,"long_name",trim(title_z1))
       if (status /= nf90_noerr) then 
          error_str= 'Error in nf90_atrdef2a'
          goto 400
       endif
       status=nf90_put_att(ncid,HeightId,"units",trim(units_z1))
       if (status /= 0) then 
          error_str= 'Error in nf90_atrdef2b'
          goto 400
       endif
       status=nf90_put_att(ncid,HeightId,"plot_title",trim(title_plot_z1))
       if (status /= 0) then 
          error_str= 'Error in nf90_atrdef2c'
          goto 400
       endif
    endif
    !
    status=nf90_put_att(ncid,QuantId1,"long_name","absolute difference "//trim(title1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5a'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId1,"units",trim(units_1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5b'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId1,"plot_title","difference "//trim(title1)//'['//trim(units_1)//']')
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5c'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId2,"long_name","relative difference "//trim(title1))
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5a'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId2,"units","%")
    if (status /= 0) then 
       error_str= 'Error in nf90_atrdef5b'
       goto 400
    endif
    status=nf90_put_att(ncid,QuantId2,"plot_title","relative difference "//trim(title1)//'[%]')
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
    !
    status = nf90_put_var(ncid, XDistId, along_track1(1:nx1))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var1'
       goto 400
    endif
    !
    status = nf90_put_var(ncid, YDistId, cross_track1(1:ny1))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var2'
       goto 400
    endif
    !
    status = nf90_put_var(ncid, XsceneId, x_scene1(1:nx1))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var3'
       goto 400
    endif
    !
    status = nf90_put_var(ncid, YsceneId, y_scene1(1:nx1))    
    if (status /= 0) then 
       error_str= 'Error in nf90_put_var4'
       goto 400
    endif
    !
    if (hor==0) then
       status = nf90_put_var(ncid, HeightId, z1(1:nz1))    
       if (status /= 0) then 
          error_str= 'Error in nf90_put_var5'
          goto 400
       endif
    endif
    !
    if (hor==0) then
       !
       status = nf90_put_var(ncid, QuantId1, data_diff_abs(1:nx1,1:nz1))    
       if (status /= 0) then 
          error_str= 'Error in nf90_put_var6'
          goto 400
       endif
       status = nf90_put_var(ncid, QuantId2, data_diff_rel(1:nx1,1:nz1))    
       if (status /= 0) then 
          error_str= 'Error in nf90_put_var7'
          goto 400
       endif
       !
    else
       !
       status = nf90_put_var(ncid, QuantId1, data_diff_abs(1:nx1,1:ny1))    
       if (status /= 0) then 
          error_str= 'Error in nf90_put_var6'
          goto 400
       endif
       status = nf90_put_var(ncid, QuantId2, data_diff_rel(1:nx1,1:ny1))    
       if (status /= 0) then 
          error_str= 'Error in nf90_put_var7'
          goto 400
       endif
       !
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
end Program L1_L2_compare
