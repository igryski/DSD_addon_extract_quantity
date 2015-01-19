!!
!! *PROGRAM* L1_L2_rebin
!!
!! USAGE: L1_L2_rebin input.nc output.nc dx(km) dz(dy)(km) 
!!
!! @version ecsim1.3
!!
!! *SRC_FILE*
!!
!! product_tools/tools/src/L1_L2_rebin.f90
!!
!! *LAST CHANGES*
!! 
!! -Oct 12, 2008: D.D. : Fixed so it now handels the output of extract_quantity_hor properly
!! -Mar 12, 2008: D.D. : Fixed so it now processes x_scene and y_scene in 1-d cases. 
!! -Mar 10, 2008: D.D. : Modified the header and Copyright info
!! -Feb 04, 2008: D.D. : Minor changes to error messages
!! -Oct 04, 2007: D.D. : Extended length of character strings in get_data and changed ESIM_HOME==>ECSIM_HOME
!!
!! *DESCRIPTION*
!!
!! L1_L2_rebin is used to change the resolution of the lidar and radar output files or the output of extract_quantity 
!! or extract_quantity_hor. The program will read the
!! data from the netcdf file and return the data in the new resolution in the output netcdf file. It will retrieve 
!! the along_track and height for vertical slabs (made with extract_quantity) and the x-scene and y-scene for the 
!! horizontal slabs (made with extract_quantity_hor).
!! 
!! The following command line arguments are expected:
!! 
!! -input.nc  : netcdf file 1 to read from.
!! -output.nc : netcdf file to write the data to
!! -dx              : horizontal resolution of the rebinned file.  A smaller than the 
!!                    original resolution or negative value uses the original resolution
!! -dz(dy)          : vertical resolution of the rebinned file (or 'y' resolution in the 
!!                    case of horizontal data).  A smaller than the 
!!                    original resolution or negative value uses the original resolution.
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
Program L1_L2_rebin
  !!
  Use average
  Use physical_parameters
  Use write_messages
  Use netcdf
  !
  Implicit None
  !
  !------------------
  ! Input variables
  !------------------
  !
  Character(len=200)             :: outfilename
  Character(len=200)             :: infilename
  real                           :: dx_ins,dz_ins
  !
  !----------------
  ! Main vaiables
  !----------------
  !
  Real,Dimension(:),Allocatable              :: x,y,dist,z,x_ori,y_ori,z_ori
  Real,Dimension(:),Allocatable              :: x_scene,y_scene
  Real,Dimension(:),Allocatable              :: x_scene_orig,y_scene_orig
  Real,Dimension(:,:),Allocatable              :: x_scene_hor,y_scene_hor
  Real,Dimension(:,:),Allocatable              :: x_scene_hor_orig,y_scene_hor_orig
  Real,Dimension(:),Allocatable              :: data11,data12
  ! dummy arrays to be rebinned
  Real,Dimension(:,:),Allocatable            :: data21,data22
  Real,Dimension(:,:,:),Allocatable          :: data31,data32
  Real,Dimension(:,:,:,:),Allocatable        :: data41,data42
  !
  Integer,Dimension(:),Allocatable           :: Var_dims 
  !
  Real,Dimension(:),Allocatable              :: z_ins,x_ins,y_ins ! output resolution  (km)
  character(len=nf90_max_name), allocatable, dimension(:) :: DimNames, VarNames
  integer, allocatable, dimension(:) :: DimLengths,DimLengths_out
  Integer                                    :: nz_ins,nx_ins,ixx2,ixx3,idim
  character(len=nf90_max_name) :: name, err_msg
  !
  !-------------------------
  ! Misc working variables 
  !-------------------------
  !
  Integer                              :: ncid,ncid2
  integer                              :: XDistId(30), YDistId(30)
  Integer                              :: status,i,ix,iy,iz,isc,j,i1,i2, nDims, nVars,ivar
  Integer                              :: ix1,ix2,iy1,iy2,iix,iyy,ixx,izz,itel1,itel2
  Integer                              :: ixs,iys
  Integer                              :: ix_out_min,ix_out_max,nx,ny,nz
  Integer                              :: iy_out_min,iy_out_max,nx_t,ny_t
  Integer                              :: n1_i,n1_o,n2_i,n2_o,n3_i,n3_o,n4_i,n4_o
  Integer                              :: dimid(10),idimxx,idimzz,idimyy
  Integer                              :: vert_hor
  real                                 :: x1,y1,z1
  real                                 :: x2,y2,z2
  Character(len=2)                     :: d_under
  Character(len=80)                    :: title, nc_title,units,plot_title
  Character(len=80),Dimension(30)      :: title_1d,title_plot_1d,units_1d
  Character(len=80),Dimension(30)      :: title_2d,title_plot_2d,units_2d
  integer,dimension(30,4)              :: dims_2d
  Real                                 :: dx_ins_orig
  Real                                 :: dz_ins_orig
  !
  character(len=200)                   :: ecsim_home,error_str
  !
  !-------------------
  ! Weighting related
  ! variables
  !-------------------
  !
!  real,dimension(:),pointer     :: xi,yi,ri
!  Real                          :: xo,yo,dr,xo1,yo1,dx_ins
!  Integer                       :: ni,nshots

! ***************************
  character(len=200)    :: tempchar, tempname, tempname2
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
       Integer,Intent(in)                     :: n
       Character(len=*),Intent(out)           :: arg
     End Subroutine getarg
     !   
  End Interface
  !
!  nullify(xi)
!  nullify(yi)
!  nullify(ri)
  d_under='d_'
  !
  !-----------------------------------------
  ! GET  environment variable ECSIM_HOME 
  !-----------------------------------------
  !
  call getenv('ECSIM_HOME', ecsim_home)
  !ecsim_home = trim(adjustl(esim_home))//'/models/'
  !call getenv('SCATT_LIB', scatt_lib)
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
  infilename=trim(adjustl(ecsim_home))//infilename
  outfilename=trim(adjustl(ecsim_home))//outfilename
  !
  !-----------------------------------
  ! Initalize some stuff regarding the 
  ! resolution and the path
  !-----------------------------------
  !
  call read_data_ncdf(11,infilename,status,plot_title)
  if (status /= 0) then
     call write_error('Error in read_data_ncdf')
     call exit(status)
  endif
  status = nf90_create(outfilename, 0, ncid2)
  if (status /= 0) then 
     error_str='error in nf90_create'
     goto 200
  endif

  !
  !-----------------------------------------------------
  ! Build the instrument resolution verticle coordinate
  !-----------------------------------------------------
  !
  !
  z1=min(z_ori(1),z_ori(nz))
  z2=max(z_ori(1),z_ori(nz))
  x1=min(x_ori(1),x_ori(nx))
  x2=max(x_ori(1),x_ori(nx))
  !
  dx_ins_orig=x_ori(2)-x_ori(1)
  dz_ins_orig=z_ori(2)-z_ori(1)
  !
  !
  if ((dx_ins.le.(x2-x1)/nx)) then 
     call write_info( '*************************************************')
     call write_info( 'Warning using the original along_track resolution')
     call write_info( '*************************************************')
     nx_ins=nx
     Allocate(x_ins(1:nx_ins))
     x_ins=x_ori
     dx_ins=(x2-x1)/nx
  else
     nx_ins=int((x2-x1)/dx_ins)
     if ((x2-x1)/dx_ins.gt.nx_ins) nx_ins=nx_ins+1
     Allocate(x_ins(1:nx_ins))
     do ix=1,nx_ins
        x_ins(ix)=(ix-1)*dx_ins+x1+dx_ins/2.0-dx_ins_orig/2.
     enddo
  endif
  if ((dz_ins.le.(z2-z1)/nz)) then
     call write_info( '********************************************')
     call write_info( 'Warning using the original height resolution')
     call write_info( '********************************************')
     nz_ins=nz
     dz_ins=(z2-z1)/nz
     Allocate(z_ins(1:nz_ins))
     z_ins=z_ori
  else
     nz_ins=int((z2-z1)/dz_ins)+1
     Allocate(z_ins(1:nz_ins))
     do iz=1,nz_ins
        z_ins(iz)=(iz-1)*dz_ins+z1+dz_ins/2.0-dz_ins_orig/2.
     enddo
  endif
  !
  ! Open the output file and start the basic write out
  !
  call write_results_ncdf(outfilename,nx_ins,nz_ins)
  !
  ! All the variables will now be read in, rebinned when neccesary and written to the output file.
  !
  do ivar=1,nvars
     title=title_1d(ivar)
     ! First the along_track and height are written
     if (ivar.eq.ixx) then
        status = nf90_put_var(ncid2, XDistId(ixx), x_ins(1:nx_ins))
     elseif (ivar.eq.izz) then
        status = nf90_put_var(ncid2, XDistId(izz), z_ins(1:nz_ins)) 
     elseif (ivar.eq.iyy) then
        status = nf90_put_var(ncid2, XDistId(iyy), z_ins(1:nz_ins)) 
     elseif (ivar.eq.ixs) then
        if (vert_hor==1) then
           allocate(data11(1:nx_ins))
           call linpol(x_ori,x_scene_orig,nx,x_ins,data11,nx_ins)
           status = nf90_put_var(ncid2, XDistId(ixs), data11) 
           deallocate(data11)
        else
           allocate(data22(1:nx_ins,1:nz_ins))
           do j=1,nx_ins
              call linpol(x_scene_hor_orig(j,:),nx,x_ins,data22(j,:),nx_ins)
           enddo
           status = nf90_put_var(ncid2, XDistId(ixs), data22) 
           deallocate(data11)
        endif
     elseif (ivar.eq.iys) then
        !
        if (vert_hor==1) then
           allocate(data11(1:nx_ins))
           call linpol(x_ori,y_scene_orig,nx,x_ins,data11,nx_ins)
           status = nf90_put_var(ncid2, XDistId(iys), data11) 
           deallocate(data11)
        else
           allocate(data22(1:nx_ins,1:nz_ins))
           do j=1,nz_ins
              call linpol(y_scene_hor_orig(:,j),nz,z_ins,data22(:,j),nz_ins)
           enddo
           status = nf90_put_var(ncid2, XDistId(iys), data22) 
           deallocate(data22)
        endif
     else
        !
        ! All one dimensional arrays are read in and written
        !
        if (var_dims(ivar).eq.1) then
           n1_i=Dimlengths_out(dims_2d(ivar,1))
           n1_o=Dimlengths(dims_2d(ivar,1))
           allocate(data11(1:n1_o))
           allocate(data12(1:n1_o))
           status = nf90_get_var(ncid, iVar, data11)
           data12=data11
           status = nf90_put_var(ncid2, XDistId(ivar), data12(1:n1_o))  
           deallocate(data11)
           deallocate(data12)
           !
           ! Two dimensional arrays are read in, checked what the dimensions are and rebinned when they 
           ! are related to <x> or <z>. A check is made if the order of the dimensions is <x,z> or <z,x>
           !
        elseif (var_dims(ivar).eq.2) then
           n1_o=Dimlengths(dims_2d(ivar,1))
           n2_o=Dimlengths(dims_2d(ivar,2))
           allocate(data21(1:n1_o,1:n2_o))
           status = nf90_get_var(ncid, iVar, data21)
           !

           !           print*,'xxxxx',n1_o,n2_o
           if ((idimxx.eq.dims_2d(ivar,1))) then 
              If (idimzz.eq.dims_2d(ivar,2)) then 
                 Allocate(data22(1:nx_ins,1:nz_ins))
                 if (trim(title(1:2)).eq.d_under) data21=data21**2
                 Call average_2d('a',n1_o,n2_o,x_ori,z_ori,data21,nx_ins,nz_ins,x_ins,z_ins,data22)
                 if (trim(title(1:2)).eq.d_under) data22=sqrt(data22)
                 status = nf90_put_var(ncid2, XDistId(ivar), data22(1:nx_ins,1:nz_ins))  
                 deallocate(data22)
              else
                 Allocate(data22(1:nx_ins,1:n2_o),data12(1:nx_ins),data11(1:n1_o))
                 do i1=1,n2_o
                    data11(:)=data21(:,i1)
                     if (trim(title(1:2)).eq.d_under) data11=data11**2
                     Call average_1d('a',n1_o,x_ori,data11,nx_ins,x_ins,data12)
                     if (trim(title(1:2)).eq.d_under) data12=sqrt(data12)
                     data22(:,i1)=data12
                  enddo
                 status = nf90_put_var(ncid2, XDistId(ivar), data22(1:nx_ins,1:n2_o))
                 deallocate(data11,data12,data22)
              endif
           elseif ((idimxx.eq.dims_2d(ivar,2)).and.(idimzz.eq.dims_2d(ivar,1))) then 
              Allocate(data22(1:nz_ins,1:nx_ins))
              if (trim(title(1:2)).eq.d_under) data21=data21**2
!              print*,n2_o,n1_o
!              print*,nz_ins,nx_ins
              Call average_2d('a',n1_o,n2_o,z_ori,x_ori,data21,nz_ins,nx_ins,z_ins,x_ins,data22)
              if (trim(title(1:2)).eq.d_under) data22=sqrt(data22)
              status = nf90_put_var(ncid2, XDistId(ivar), data22(1:nz_ins,1:nx_ins))  
              deallocate(data22)
           else
              call write_error('The along_track or height dimension are not the first dimension')
              call write_error('of the '//trim(title_1d(ivar))//' variable')
              call exit(1)
           endif
           !
           deallocate(data21)
           !
           ! Three dimensional arrays are read in, checked what the dimensions are and rebinned when they 
           ! are related to <x> or <z>. A check is made if the order of the dimensions is <x,z,i> or <z,x,i>
           !
        elseif  (var_dims(ivar).eq.3) then
           n1_o=Dimlengths(dims_2d(ivar,1))
           n2_o=Dimlengths(dims_2d(ivar,2))
           n3_o=Dimlengths(dims_2d(ivar,3))
           allocate(data31(1:n1_o,1:n2_o,1:n3_o))
           status = nf90_get_var(ncid, iVar, data31)
           !
           if ((idimxx.eq.dims_2d(ivar,1)).and.(idimzz.eq.dims_2d(ivar,2))) then 
              allocate(data32(1:nx_ins,1:nz_ins,1:n3_o),data21(1:n1_o,1:n2_o),data22(1:nx_ins,1:nz_ins))
              do i1=1,n3_o
                 data21(:,:)=data31(:,:,i1)
                 if (trim(title(1:2)).eq.d_under) data21=data21**2
                 Call average_2d('a',n1_o,n2_o,x_ori,z_ori,data21,nx_ins,nz_ins,x_ins,z_ins,data22)
                 if (trim(title(1:2)).eq.d_under) data22=sqrt(data22)
                 data32(:,:,i1)=data22(:,:)
                 data21=0.0
                 data22=0.0
              enddo
              status = nf90_put_var(ncid2, XDistId(ivar), data32(1:nx_ins,1:nz_ins,1:n3_o))
              deallocate(data21,data22,data32)
           elseif ((idimxx.eq.dims_2d(ivar,2)).and.(idimzz.eq.dims_2d(ivar,1))) then 
              allocate(data32(1:nz_ins,1:nx_ins,1:n3_o),data21(1:n1_o,1:n2_o),data22(1:nz_ins,1:nx_ins))
              do i1=1,n3_o
                 data21(:,:)=data31(:,:,i1)
                 if (trim(title(1:2)).eq.d_under) data21=data21**2
                 Call average_2d('a',n1_o,n2_o,z_ori,x_ori,data21,nz_ins,nx_ins,z_ins,x_ins,data22)
                 if (trim(title(1:2)).eq.d_under) data22=sqrt(data22)
                 data32(:,:,i1)=data22(:,:)
                 data21=0.0
                 data22=0.0
              enddo
              status = nf90_put_var(ncid2, XDistId(ivar), data32(1:nz_ins,1:nx_ins,1:n3_o))
              deallocate(data21,data22,data32)
           else
              call write_error('The along_track or height dimension are not the first dimension')
              call write_error('of the '//trim(title_1d(ivar))//' variable')
              call exit(1)
           endif
           !
           deallocate(data31)
        endif
     endif
  enddo

  status = nf90_close(ncid2)
  if (status /= 0) then 
     error_str= 'Error in nf90_close'
     goto 200
  endif

  status = nf90_close(ncid)
  err_msg =  'Error in nf90_close'
  if (status.gt.0) goto 200
  
  call Write_info('******************L1_L2_REBIN FINISHED******************')
  call exit(0)
200 call Write_error('*****************An error in L1_L2_rebin occured******************')
  call Write_error(err_msg)
  call exit(1)
  !  
  ! 
Contains
  !
  Subroutine read_data_ncdf(id,inputfile, exit_status,plot_title)
    !
    use typeSizes
  !  use netcdf
    !
    integer,intent(in)          :: id
    Character(len=*),intent(in) :: inputfile
    integer, intent(out) :: exit_status
    character(len=80),intent(out):: plot_title
    !
    ! NEW variables
    Integer :: status, io_mode, iDim ,iVar, length 
    character(len=nf90_max_name) :: name, err_msg
    integer :: found,iwork1,iwork2,iwork3(10)
    !
    exit_status = 0
    !
    io_mode = 0  ! default: nf90_nowrite
    status = nf90_open(inputfile, io_mode, ncid)
    err_msg = 'Error in nf90_open '//inputfile
    if (status/=0) then
       call write_error(err_msg)
       exit_status=1 
       goto 999
    endif
    ! Get info on
    ! - nDims
    ! - nVars
    status  = nf90_inquire(ncid, nDims, nVars)   
    if (status /= 0) stop 'Error in nf90_inquire'
    !
    ! Create arrays that contain Dimension Info   
    allocate(DimNames(1:nDims))
    allocate(DimLengths(1:nDims))

    ! Get info on the dimensions:
    do iDim = 1,nDims
       !
       status=nf90_inquire_dimension(ncid, iDim, name, length)
       !       
       if (status/=0) then
          err_msg = 'Error in nf90_inquire_dimension '
          call write_error(err_msg)
          exit_status=1 
          goto 999
       endif
       !
       DimNames(iDim)   = name
       DimLengths(iDim) = length
       Select Case (name)
       Case ('nx') ;  nx = length
       Case ('ny') ;  ny = length
       Case ('nz') ;  nz = length
       End Select
       !
    end do
    allocate(var_dims(nvars))
    !
    iyy=-1
    izz=-1
    ixx=-1  
    ixs=-1
    iys=-1
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
       var_dims(ivar)=iwork2
       title_1d(ivar)=trim(title)
       title_plot_1d(ivar)=trim(plot_title)
       units_1d(ivar)=trim(units)
       dims_2d(ivar,1:iwork2)=iwork3(1:iwork2)
       
       if (title.eq.'along_track') then
          allocate(x_ori(1:nx))
          ixx=iVar
          idimxx=iwork3(1)
          status = nf90_get_var(ncid, iVar, x_ori)
       elseif (title.eq.'cross_track') then
          allocate(y_ori(1:ny))
          idimzz=iwork3(1)
          iyy=iVar
          status = nf90_get_var(ncid, iVar, y_ori)
       elseif (title.eq.'height') then
          allocate(z_ori(1:nz))
          izz=iVar
          idimzz=iwork3(1)
          status = nf90_get_var(ncid, iVar, z_ori)
       endif
       !
       if (status/=0) then
          call write_error(err_msg)
          exit_status=1 
          goto 999
       endif
    enddo
    !
    vert_hor=1
    !
    if (.not.allocated(z_ori)) then    ! Assume this file is x-y
       if (allocated(y_ori)) then
          nz=ny
          allocate(z_ori(1:ny))
          z_ori=y_ori
          vert_hor=0
       else
          call write_error('No z or y coordinate in the input ncdf file')
          call exit(1)
       endif
    endif
    !
    if (vert_hor==1) then
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
          var_dims(ivar)=iwork2
          title_1d(ivar)=trim(title)
          title_plot_1d(ivar)=trim(plot_title)
          units_1d(ivar)=trim(units)
          dims_2d(ivar,1:iwork2)=iwork3(1:iwork2)
          !
          if (title.eq.'x_scene') then
             allocate(x_scene_orig(1:nx))
             ixs=iVar
             idimxx=iwork3(1)
             status = nf90_get_var(ncid, iVar, x_scene_orig)
          elseif (title.eq.'y_scene') then
             allocate(y_scene_orig(1:nx))
             iys=iVar
             idimyy=iwork3(1) 
             status = nf90_get_var(ncid, iVar, y_scene_orig)
          endif
       enddo
       !
    else
       do iVar=1,nVars
          if (title.eq.'x_scene') then
             allocate(x_scene_hor_orig(1:nx,1:ny))
             ixs=iVar
             idimxx=iwork3(1)
             status = nf90_get_var(ncid, iVar, x_scene_hor_orig)
          elseif (title.eq.'y_scene') then
             allocate(y_scene_hor_orig(1:nx,1:ny))
             iys=iVar
             idimzz=iwork3(1)
             status = nf90_get_var(ncid, iVar, y_scene_hor_orig)
          endif
       enddo
       !
    endif
    !
999 return
  end Subroutine read_data_ncdf


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
    error_str='Usage:L1_L2_rebin INPUTFILENAME1 OUTPUTFILENAME dx dz '
    error_str2=error_str
    !
    status=0
    !
    nargs=iargc() 
!    print*,nargs
    If (nargs.lt.4.or.nargs.gt.5) Then
       status=1
       error_str2='Wrong number of arguments'
       goto 100
    Endif
    !
    If (status == 0) Then
       Call getarg(1,arg_str)
       error_str2='Error in infilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) infilename
    Endif
    !
    !
    If (status == 0) Then
       Call getarg(2,arg_str)
       error_str2='Error in outfilename specification'
       Read(arg_str,'(a160)',iostat=status,err=100) outfilename
    Endif
    !
    If (status == 0) Then
       Call getarg(3,arg_str)
       error_str2='Error in horizontal resolution specification'
       Read(arg_str,*,iostat=status,err=100) dx_ins
    Endif
    !
    If (status == 0) Then
       Call getarg(4,arg_str)
       error_str2='Error in vertical resolution specification'
       Read(arg_str,*,iostat=status,err=100) dz_ins
    Endif
    !
    vert_hor=0
!    if (nargs.eq.5.and.status==0) then
!       Call getarg(5,arg_str)
!       error_str2='Error in vertical_horizontal specification'
!       Read(arg_str,*,iostat=status,err=100) vert_hor
!    endif
    !
100 If (status.ne.0) Then 
       call Write_error(error_str)
       call Write_error( error_str2)
    Endif
    
    !
  End Subroutine get_data
  !
  !
  !
  Subroutine write_results_ncdf(filename,nx_new,nz_new)
    !
    use typeSizes
    use netcdf
    Use write_messages
    !
    implicit none
    !
    Character(len=*),intent(in)    :: filename
    Integer,intent(in)             :: nx_new
    Integer,intent(in)             :: Nz_new
    Real                :: x(nx_new),y(nx_new),dist(nx_new),z(nz_new)
    Real                :: Quant(Nz_new,nx_new)
    !Character(len=*),intent(in)    :: nc_title,title,units,plot_title
    !
    Integer                        :: i,j,ivar,xidim,zidim
    !
    Integer :: status,n,idim
    !
    integer :: GroundDist, Altitude,crossdist
    integer :: HeightId, QuantId
    Character(len=190)               :: error_str
    !
    n=nx_new
    !
    !
    ! Assume the file does not exist;
    ! If you want to check, you need to use the nf90_open function
    allocate(dimlengths_out(1:ndims))
    
    do idim=1,ndims
       if (DimNames(iDim)=='nx') then
          xidim=idim
          dimlengths_out(idim)=nx_new
!          print*,nx_new,nz_new
          status = nf90_def_dim(ncid2, "nx", nx_new, GroundDist)    
          if (status /= 0) then 
             error_str= 'error in nf90_def_dim: nx'
             goto 400
          endif
       elseif (DimNames(iDim)=='ny') then
          xidim=idim
          dimlengths_out(idim)=nx_new
!          print*,nx_new,nz_new
          status = nf90_def_dim(ncid2, "ny", nx_new, crossDist)    
          if (status /= 0) then 
             error_str= 'error in nf90_def_dim: ny'
             goto 400
          endif
       elseif (DimNames(iDim)=='nz') then
          zidim=idim
          dimlengths_out(idim)=nz_new
          status = nf90_def_dim(ncid2, "nz", nz_new,  Altitude)    
          if (status /= 0) then 
             error_str= 'error in nf90_def_dim: nz'
             goto 400
          endif
       else     
          dimlengths_out(idim)=dimlengths(idim)
          status = nf90_def_dim(ncid2,trim(DimNames(idim)), dimlengths(idim), dimid(idim))
           if (status /= 0) then 
             error_str= 'error in nf90_def_dim remainder'
             goto 400
          endif
      endif
    enddo
    !
    ! Defining the variables
    !
    do ivar=1,nvars
       if (var_dims(ivar).eq.1) then
          status = nf90_def_var(ncid2, trim(title_1d(ivar)), NF90_FLOAT,(/ dims_2d(ivar,1)/), XDistId(ivar))
       elseif (var_dims(ivar).eq.2) then
          status = nf90_def_var(ncid2, trim(title_1d(ivar)), NF90_FLOAT,(/dims_2d(ivar,1),dims_2d(ivar,2) /), XDistId(ivar))
        
       elseif (var_dims(ivar).eq.3) then
          status = nf90_def_var(ncid2, trim(title_1d(ivar)), NF90_FLOAT,(/dims_2d(ivar,1),&
               dims_2d(ivar,2),dims_2d(ivar,3) /),XDistId(ivar))
          
       elseif (var_dims(ivar).eq.4) then
          status = nf90_def_var(ncid2, trim(title_1d(ivar)), NF90_FLOAT,(/dims_2d(ivar,1),&
               dims_2d(ivar,2),dims_2d(ivar,3),dims_2d(ivar,4) /),XDistId(ivar))
       endif
       
       if (status /= 0) then 
          error_str= 'error in nf90_def_var1'
          goto 400
       endif
    enddo
    !
    ! Defining attributes
    !
    do ivar=1,nvars
       status=nf90_put_att(ncid2,XDistId(ivar),"long_name",title_1d(ivar))
       status=nf90_put_att(ncid2,XDistId(ivar),"units",units_1d(ivar))
       status=nf90_put_att(ncid2,XDistId(ivar),"plot_title",title_plot_1d(ivar))
    enddo
    if (status /= 0) then 
       error_str= 'error in nf90_def_atr'
       goto 400
    endif
    !    
    status = nf90_enddef(ncid2)
    if (status /= 0) then 
       error_str= 'Error in nf90_enddef'
       goto 400
    endif

400 if (status.ne.0) then
       call write_error(error_str)
       call exit(1)
    endif
    !
    Return
    
  End Subroutine write_results_ncdf
  !
end Program L1_L2_rebin





  
