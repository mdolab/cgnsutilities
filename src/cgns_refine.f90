program cgns_refine

  !use cgnsGrid

  implicit none
  include 'cgnslib_f.h'

  integer  CellDim, PhysDim
  integer nbases, nzones, in_dims(9),out_dims(9)

  integer ier, zonetype
  integer n,nn,mm,i,j,k,idim

  integer cg_in, cg_out, base, zone, zoneCounter

  integer :: ptSetType, normalDataType
  integer :: sizes(9),coordID,blockStart(3),blockEnd(3)

  integer donor_range(6),transform(3)
  integer nbocos,n1to1,bocotype,BCout
  integer NormalIndex(3), NormalListFlag, ndataset,datatype
  integer ptset_type, npnts, pnts(3,2),pnts_donor(3,2),ncon

  character*32 in_filename,out_filename
  character*32 basename, zonename,boconame
  character*32 connectname,donorname, famname
  real(kind=8), allocatable,dimension(:,:,:,:) :: data3d_in,data3d_out
  real(kind=8), allocatable,dimension(:,:,:) :: data3d_temp

  real(kind=8)::  MI(64,64)
  real(kind=8):: data_double(6)

  N = IARGC ()

  if (N .ne. 2) then
     print *,'Errors: cgns_refine must be called with TWO arguments:'
     print *,'./cgns_refine cgns_in_file.cgns cgns_out_file.cgns'
     stop
  end if

  call getarg(1, in_filename)
  call getarg(2, out_filename)

  call cg_open_f(in_filename,CG_MODE_READ, cg_in, ier)

  if (ier .eq. CG_ERROR) call cg_error_exit_f
  print *,'Reading input file: ',in_filename

  call cg_nbases_f(cg_in, nbases, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  if (nbases .lt. 1) then
     print *, 'Error: No bases found in CGNS file'
  else if(nbases .gt. 1) then
     print *, 'Warning: More than one base found. Using only the first base'
  end if

  ! Explicitly use ONLY the first base
  base = 1

  call cg_base_read_f(cg_in, base, basename, CellDim, PhysDim, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f
  if (cellDim .ne. 3 .or. PhysDim .ne. 3) then
     print *,'The Cells must be hexahedreal in 3 dimensions'
     stop
  end if

  ! Goto Base Node
  call cg_goto_f(cg_in, base, ier, 'end')
  call cg_nzones_f(cg_in, base, nzones, ier)

  ! Open the output file:
  call cg_open_f(out_filename,CG_MODE_WRITE, cg_out, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  call cg_base_write_f(cg_out,"Base#1", Celldim,Physdim, base, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  ! Get hte inverse interpolation matrix:
  call Minv(MI)

  do nn=1,nzones
     print *,'Processing zone...',nn
     call cg_zone_read_f(cg_in, base, nn, zonename, in_dims, ier)

     ! Double up the dimensions:
     do i=1,3
        out_dims(i) = (in_dims(i)-1)*2 + 1
     end do
     do i=4,6
        out_dims(i) = in_dims(i)*2
     end do

     ! Write new grid zone
     call cg_zone_write_f(cg_out,base,zonename,out_dims,Structured,zoneCounter,ier)

     ! Allocate memory for incomming grid
     allocate(data3d_temp(in_dims(1),in_dims(2),in_dims(3)))
     allocate(data3d_in (3,in_dims(1),in_dims(2),in_dims(3)))
     allocate(data3d_out(3,out_dims(1),out_dims(2),out_dims(3)))

     ! Initial data arrays
     data3d_in(:,:,:,:) = 0.0
     data3d_out(:,:,:,:) = 0.0

     blockStart = (/1,1,1/)
     blockEnd   = in_dims(1:3)

     ! Read/Write each Coordinate:

     ! X coordinate Read
     call cg_coord_read_f(cg_in,base,nn,'CoordinateX',RealDouble,&
          blockStart,blockEnd,data3d_temp,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     call addData(data3d_temp,data3d_in,in_dims(1),in_dims(2),in_dims(3),1)

     ! Y coordinate Read
     call cg_coord_read_f(cg_in,base,nn,'CoordinateY',RealDouble,&
          blockStart,blockEnd,data3d_temp,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     call addData(data3d_temp,data3d_in,in_dims(1),in_dims(2),in_dims(3),2)

     ! Z coordinate Read
     call cg_coord_read_f(cg_in,base,nn,'CoordinateZ',RealDouble,&
          blockStart,blockEnd,data3d_temp,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     call addData(data3d_temp,data3d_in,in_dims(1),in_dims(2),in_dims(3),3)

     ! Interpolate all the Data:
     call interp2(data3d_in,data3d_out,in_dims(1),in_dims(2),in_dims(3),MI)

     ! Reallocate temp array to the size of the output
     deallocate(data3d_temp)
     allocate(data3d_temp(out_dims(1),out_dims(2),out_dims(3)))

     ! Write out coordinates:
     
     call convertOutput(data3d_out,data3d_temp,out_dims(1),out_dims(2),&
          out_dims(3),1)

     call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
          'CoordinateX',data3d_temp, coordID,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     call convertOutput(data3d_out,data3d_temp,out_dims(1),out_dims(2),&
          out_dims(3),2)

     call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
          'CoordinateY',data3d_temp, coordID,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     call convertOutput(data3d_out,data3d_temp,out_dims(1),out_dims(2),&
          out_dims(3),3)
     
     call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
          'CoordinateZ',data3d_temp, coordID,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     ! Ditch the allocated data
     deallocate(data3d_in,data3d_out,data3d_temp)

     ! Next to the BC's if there are any:

     call cg_nbocos_f(cg_in, base, nn, nBocos,ier)

     do mm=1,nBocos
        ! Get Boundary Condition Info
        call cg_boco_info_f(cg_in, base, nn, mm , boconame,bocotype,&
             ptset_type,npnts,NormalIndex,NormalListFlag,datatype,ndataset,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        call cg_boco_read_f(cg_in, base, nn, mm, pnts,data_double, ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        ! Process the point range:
        do i=1,3
           do j=1,2
              pnts(i,j) = (pnts(i,j)-1)*2+1
           end do
        end do

        ! Write Boundary Condition Info:
        call cg_boco_write_f(cg_out,base,nn,boconame, bocotype, PointRange, &
             2, pnts, BCout ,  ier )
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        ! Get/Write the Family Info
        call cg_goto_f(cg_in,base,ier,"Zone_t",nn,"ZoneBC_t",1,"BC_t",mm,"end")
        if (ier == 0) then ! Node exits
           call cg_famname_read_f(famName, ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f

           call cg_goto_f(cg_out,base,ier,'Zone_t', nn,"ZoneBC_t", 1,&
                "BC_t", BCOut, "end")
           call cg_famname_write_f(famName, ier)
           if (ier .eq. CG_ERROR) call cg_error_exit_f
        end if
     end do

     ! Next read the 1to1 Connectivity if there are any:
     call cg_n1to1_f(cg_in, base, nn, n1to1, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f
  
     do mm=1,n1to1
        call cg_1to1_read_f(cg_in, base, nn, mm, connectName, donorname, &
             pnts,pnts_donor, transform, ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        ! Process the point range:
        do i=1,3
           do j=1,2
              pnts(i,j) = (pnts(i,j)-1)*2+1
              pnts_donor(i,j) = (pnts_donor(i,j)-1)*2+1
           end do
        end do

        call cg_1to1_write_f(cg_out,base,nn,connectName,donorname,&
             pnts,pnts_donor,transform,nCon,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f
     end do ! 1ot1 Loop
  end do ! Zone Loop

  ! Close both CGNS files
  call cg_close_f(cg_in, ier)
  call cg_close_f(cg_out,ier)

end program cgns_refine

subroutine addData(temp,data_in,nx,ny,nz,idx)
  implicit none

  ! Input/Output
  integer :: nx,ny,nz ! These are the GRID dimensions:
  real(kind=8),dimension(nx,ny,nz), intent(in) :: temp
  real(kind=8),dimension(3,nx,ny,nz),intent(inout) :: data_in
  integer, intent(in) :: idx

  ! Working
  integer :: i,j,k

  do k=1,nz
     do j=1,ny
        do i=1,nx
           data_in(idx,i,j,k) = temp(i,j,k)
        end do
     end do
  end do

end subroutine addData


subroutine convertOutput(data_out,temp,nx,ny,nz,idx)
  implicit none

  ! Input/Output
  integer :: nx,ny,nz ! These are the GRID dimensions:
  real(kind=8),dimension(3,nx,ny,nz),intent(in) :: data_out
  real(kind=8),dimension(nx,ny,nz), intent(inout) :: temp

  integer, intent(in) :: idx

  ! Working
  integer :: i,j,k

  do k=1,nz
     do j=1,ny
        do i=1,nx
           temp(i,j,k) = data_out(idx,i,j,k) 
        end do
     end do
  end do
  
end subroutine convertOutput
  
subroutine interp2(input,output,nx,ny,nz,MI)

  implicit none

  ! Input/Output
  integer :: nx,ny,nz ! These are the GRID dimensions:
  real(kind=8), dimension(3,nx,ny,nz),intent(in) :: input
  real(kind=8), dimension(3,2*nx-1,2*ny-1,2*nz-1),intent(out) :: output
  real(kind=8),dimension(:,:,:,:,:),allocatable :: fcn
  real(kind=8),dimension(:,:,:,:),allocatable :: paraS
  
  real(kind=8) :: corners(8),s(3),shp(8)
  real(kind=8) :: MI(64,64)
  real(kind=8) :: coef(3,64)
  real(kind=8) :: b(64)
  ! Working
  integer :: i,j,k,ii,jj,kk,ndim

  ! The purpose of this routine is to interpolate from the coarse
  ! "input" to the refined "output" using a higher order cubic
  ! interpolation. This will attempt to capture the the curvature
  ! information of (especially) the surface which is required to
  ! obtain the desired order of convergence. 

  ! Now we calculate the following partial derivatives (in this
  ! order): F,Fx,Fy,Fz,Fxy,Fxz,Fyz, and Fxyz
  ! Index : 1,2 ,3 ,4, 5  ,6  ,7  ,     8

  ! Run the parameterization:
  ndim = 3
  allocate(fcn(ndim,nx,ny,nz,8),paras(ndim,nx,ny,nz))

  call para3d(input,nx,ny,nz,ndim,paraS)
    
  ! Values:
  do k=1,nz
     do j=1,ny
        do i=1,nx
           do ii=1,3
              fcn(ii,i,j,k,1) = input(ii,i,j,k)
           end do
        end do
     end do
  end do
  
  ! X Derivatives
  do k=1,nz
     do j=1,ny
        do ii=1,3
           call splineDERIV(fcn(ii,:,j,k,1),paraS(1,:,j,k),fcn(ii,:,j,k,2),nx) !Fx
        end do
     end do
  end do

  ! Y Derivatives
  do k=1,nz
     do i=1,nx
        do ii=1,3
           call splineDERIV(fcn(ii,i,:,k,1),paraS(2,i,:,k),fcn(ii,i,:,k,3),ny) !Fy
           call splineDERIV(fcn(ii,i,:,k,2),paraS(2,i,:,k),fcn(ii,i,:,k,5),ny) !(Fx)y
        end do
     end do
  end do

  ! Z Derivatives
  do j=1,ny
     do i=1,nx
        do ii=1,3
           call splineDERIV(fcn(ii,i,j,:,1),paraS(3,i,j,:),fcn(ii,i,j,:,4),nz) ! Fz
           call splineDERIV(fcn(ii,i,j,:,2),paraS(3,i,j,:),fcn(ii,i,j,:,6),nz) !(Fx)z
           call splineDERIV(fcn(ii,i,j,:,3),paraS(3,i,j,:),fcn(ii,i,j,:,7),nz) !(Fy)z
           call splineDERIV(fcn(ii,i,j,:,5),paraS(3,i,j,:),fcn(ii,i,j,:,8),nz) !(fxy)z
        end do
     end do
  end do

  ! Loop over each CELL:
  do k=1,nz-1
     do j=1,ny-1
        do i=1,nx-1
           do ii = 1,3
              ! Get the coefficients for this cell
              do jj=1,8
                 b((jj-1)*8+1) = fcn(ii,i  ,j  ,k  ,jj)
                 b((jj-1)*8+2) = fcn(ii,i+1,j  ,k  ,jj)
                 b((jj-1)*8+3) = fcn(ii,i  ,j+1,k  ,jj)
                 b((jj-1)*8+4) = fcn(ii,i+1,j+1,k  ,jj)
                 b((jj-1)*8+5) = fcn(ii,i  ,j  ,k+1,jj)
                 b((jj-1)*8+6) = fcn(ii,i+1,j  ,k+1,jj)
                 b((jj-1)*8+7) = fcn(ii,i  ,j+1,k+1,jj)
                 b((jj-1)*8+8) = fcn(ii,i+1,j+1,k+1,jj)
              end do

              ! Compute the coefficients with the Matult
              coef(ii,:) =  matmul(MI,b)
           end do

           ! Fill in the 27 values for this cell
           do kk=0,2
              do jj=0,2
                 do ii=0,2
                    s = (/dble(ii)/2.0,dble(jj)/2.0,dble(kk)/2.0/)

                    call splineEval(s,coef,&
                         output(:,(i-1)*2+ii+1,(j-1)*2+jj+1,(k-1)*2+kk+1))
                 end do
              end do
           end do
        end do
     end do
  end do
  deallocate(paraS,fcn)
end subroutine interp2

subroutine Minv(M)

  implicit none

  real(kind=8), intent(out) :: M(64,64)

  M = 0.0

  M	(1,1)	=	1.0
  M	(3,1)	=	-3.0
  M	(4,1)	=	2.0
  M	(9,1)	=	-3.0
  M	(11,1)	=	9.0
  M	(12,1)	=	-6.0
  M	(13,1)	=	2.0
  M	(15,1)	=	-6.0
  M	(16,1)	=	4.0
  M	(33,1)	=	-3.0
  M	(35,1)	=	9.0
  M	(36,1)	=	-6.0
  M	(41,1)	=	9.0
  M	(43,1)	=	-27.0
  M	(44,1)	=	18.0
  M	(45,1)	=	-6.0
  M	(47,1)	=	18.0
  M	(48,1)	=	-12.0
  M	(49,1)	=	2.0
  M	(51,1)	=	-6.0
  M	(52,1)	=	4.0
  M	(57,1)	=	-6.0
  M	(59,1)	=	18.0
  M	(60,1)	=	-12.0
  M	(61,1)	=	4.0
  M	(63,1)	=	-12.0
  M	(64,1)	=	8.0
  M	(3,2)	=	3.0
  M	(4,2)	=	-2.0
  M	(11,2)	=	-9.0
  M	(12,2)	=	6.0
  M	(15,2)	=	6.0
  M	(16,2)	=	-4.0
  M	(35,2)	=	-9.0
  M	(36,2)	=	6.0
  M	(43,2)	=	27.0
  M	(44,2)	=	-18.0
  M	(47,2)	=	-18.0
  M	(48,2)	=	12.0
  M	(51,2)	=	6.0
  M	(52,2)	=	-4.0
  M	(59,2)	=	-18.0
  M	(60,2)	=	12.0
  M	(63,2)	=	12.0
  M	(64,2)	=	-8.0
  M	(9,3)	=	3.0
  M	(11,3)	=	-9.0
  M	(12,3)	=	6.0
  M	(13,3)	=	-2.0
  M	(15,3)	=	6.0
  M	(16,3)	=	-4.0
  M	(41,3)	=	-9.0
  M	(43,3)	=	27.0
  M	(44,3)	=	-18.0
  M	(45,3)	=	6.0
  M	(47,3)	=	-18.0
  M	(48,3)	=	12.0
  M	(57,3)	=	6.0
  M	(59,3)	=	-18.0
  M	(60,3)	=	12.0
  M	(61,3)	=	-4.0
  M	(63,3)	=	12.0
  M	(64,3)	=	-8.0
  M	(11,4)	=	9.0
  M	(12,4)	=	-6.0
  M	(15,4)	=	-6.0
  M	(16,4)	=	4.0
  M	(43,4)	=	-27.0
  M	(44,4)	=	18.0
  M	(47,4)	=	18.0
  M	(48,4)	=	-12.0
  M	(59,4)	=	18.0
  M	(60,4)	=	-12.0
  M	(63,4)	=	-12.0
  M	(64,4)	=	8.0
  M	(33,5)	=	3.0
  M	(35,5)	=	-9.0
  M	(36,5)	=	6.0
  M	(41,5)	=	-9.0
  M	(43,5)	=	27.0
  M	(44,5)	=	-18.0
  M	(45,5)	=	6.0
  M	(47,5)	=	-18.0
  M	(48,5)	=	12.0
  M	(49,5)	=	-2.0
  M	(51,5)	=	6.0
  M	(52,5)	=	-4.0
  M	(57,5)	=	6.0
  M	(59,5)	=	-18.0
  M	(60,5)	=	12.0
  M	(61,5)	=	-4.0
  M	(63,5)	=	12.0
  M	(64,5)	=	-8.0
  M	(35,6)	=	9.0
  M	(36,6)	=	-6.0
  M	(43,6)	=	-27.0
  M	(44,6)	=	18.0
  M	(47,6)	=	18.0
  M	(48,6)	=	-12.0
  M	(51,6)	=	-6.0
  M	(52,6)	=	4.0
  M	(59,6)	=	18.0
  M	(60,6)	=	-12.0
  M	(63,6)	=	-12.0
  M	(64,6)	=	8.0
  M	(41,7)	=	9.0
  M	(43,7)	=	-27.0
  M	(44,7)	=	18.0
  M	(45,7)	=	-6.0
  M	(47,7)	=	18.0
  M	(48,7)	=	-12.0
  M	(57,7)	=	-6.0
  M	(59,7)	=	18.0
  M	(60,7)	=	-12.0
  M	(61,7)	=	4.0
  M	(63,7)	=	-12.0
  M	(64,7)	=	8.0
  M	(43,8)	=	27.0
  M	(44,8)	=	-18.0
  M	(47,8)	=	-18.0
  M	(48,8)	=	12.0
  M	(59,8)	=	-18.0
  M	(60,8)	=	12.0
  M	(63,8)	=	12.0
  M	(64,8)	=	-8.0
  M	(2,9)	=	1.0
  M	(3,9)	=	-2.0
  M	(4,9)	=	1.0
  M	(10,9)	=	-3.0
  M	(11,9)	=	6.0
  M	(12,9)	=	-3.0
  M	(14,9)	=	2.0
  M	(15,9)	=	-4.0
  M	(16,9)	=	2.0
  M	(34,9)	=	-3.0
  M	(35,9)	=	6.0
  M	(36,9)	=	-3.0
  M	(42,9)	=	9.0
  M	(43,9)	=	-18.0
  M	(44,9)	=	9.0
  M	(46,9)	=	-6.0
  M	(47,9)	=	12.0
  M	(48,9)	=	-6.0
  M	(50,9)	=	2.0
  M	(51,9)	=	-4.0
  M	(52,9)	=	2.0
  M	(58,9)	=	-6.0
  M	(59,9)	=	12.0
  M	(60,9)	=	-6.0
  M	(62,9)	=	4.0
  M	(63,9)	=	-8.0
  M	(64,9)	=	4.0
  M	(3,10)	=	-1.0
  M	(4,10)	=	1.0
  M	(11,10)	=	3.0
  M	(12,10)	=	-3.0
  M	(15,10)	=	-2.0
  M	(16,10)	=	2.0
  M	(35,10)	=	3.0
  M	(36,10)	=	-3.0
  M	(43,10)	=	-9.0
  M	(44,10)	=	9.0
  M	(47,10)	=	6.0
  M	(48,10)	=	-6.0
  M	(51,10)	=	-2.0
  M	(52,10)	=	2.0
  M	(59,10)	=	6.0
  M	(60,10)	=	-6.0
  M	(63,10)	=	-4.0
  M	(64,10)	=	4.0
  M	(10,11)	=	3.0
  M	(11,11)	=	-6.0
  M	(12,11)	=	3.0
  M	(14,11)	=	-2.0
  M	(15,11)	=	4.0
  M	(16,11)	=	-2.0
  M	(42,11)	=	-9.0
  M	(43,11)	=	18.0
  M	(44,11)	=	-9.0
  M	(46,11)	=	6.0
  M	(47,11)	=	-12.0
  M	(48,11)	=	6.0
  M	(58,11)	=	6.0
  M	(59,11)	=	-12.0
  M	(60,11)	=	6.0
  M	(62,11)	=	-4.0
  M	(63,11)	=	8.0
  M	(64,11)	=	-4.0
  M	(11,12)	=	-3.0
  M	(12,12)	=	3.0
  M	(15,12)	=	2.0
  M	(16,12)	=	-2.0
  M	(43,12)	=	9.0
  M	(44,12)	=	-9.0
  M	(47,12)	=	-6.0
  M	(48,12)	=	6.0
  M	(59,12)	=	-6.0
  M	(60,12)	=	6.0
  M	(63,12)	=	4.0
  M	(64,12)	=	-4.0
  M	(34,13)	=	3.0
  M	(35,13)	=	-6.0
  M	(36,13)	=	3.0
  M	(42,13)	=	-9.0
  M	(43,13)	=	18.0
  M	(44,13)	=	-9.0
  M	(46,13)	=	6.0
  M	(47,13)	=	-12.0
  M	(48,13)	=	6.0
  M	(50,13)	=	-2.0
  M	(51,13)	=	4.0
  M	(52,13)	=	-2.0
  M	(58,13)	=	6.0
  M	(59,13)	=	-12.0
  M	(60,13)	=	6.0
  M	(62,13)	=	-4.0
  M	(63,13)	=	8.0
  M	(64,13)	=	-4.0
  M	(35,14)	=	-3.0
  M	(36,14)	=	3.0
  M	(43,14)	=	9.0
  M	(44,14)	=	-9.0
  M	(47,14)	=	-6.0
  M	(48,14)	=	6.0
  M	(51,14)	=	2.0
  M	(52,14)	=	-2.0
  M	(59,14)	=	-6.0
  M	(60,14)	=	6.0
  M	(63,14)	=	4.0
  M	(64,14)	=	-4.0
  M	(42,15)	=	9.0
  M	(43,15)	=	-18.0
  M	(44,15)	=	9.0
  M	(46,15)	=	-6.0
  M	(47,15)	=	12.0
  M	(48,15)	=	-6.0
  M	(58,15)	=	-6.0
  M	(59,15)	=	12.0
  M	(60,15)	=	-6.0
  M	(62,15)	=	4.0
  M	(63,15)	=	-8.0
  M	(64,15)	=	4.0
  M	(43,16)	=	-9.0
  M	(44,16)	=	9.0
  M	(47,16)	=	6.0
  M	(48,16)	=	-6.0
  M	(59,16)	=	6.0
  M	(60,16)	=	-6.0
  M	(63,16)	=	-4.0
  M	(64,16)	=	4.0
  M	(5,17)	=	1.0
  M	(7,17)	=	-3.0
  M	(8,17)	=	2.0
  M	(9,17)	=	-2.0
  M	(11,17)	=	6.0
  M	(12,17)	=	-4.0
  M	(13,17)	=	1.0
  M	(15,17)	=	-3.0
  M	(16,17)	=	2.0
  M	(37,17)	=	-3.0
  M	(39,17)	=	9.0
  M	(40,17)	=	-6.0
  M	(41,17)	=	6.0
  M	(43,17)	=	-18.0
  M	(44,17)	=	12.0
  M	(45,17)	=	-3.0
  M	(47,17)	=	9.0
  M	(48,17)	=	-6.0
  M	(53,17)	=	2.0
  M	(55,17)	=	-6.0
  M	(56,17)	=	4.0
  M	(57,17)	=	-4.0
  M	(59,17)	=	12.0
  M	(60,17)	=	-8.0
  M	(61,17)	=	2.0
  M	(63,17)	=	-6.0
  M	(64,17)	=	4.0
  M	(7,18)	=	3.0
  M	(8,18)	=	-2.0
  M	(11,18)	=	-6.0
  M	(12,18)	=	4.0
  M	(15,18)	=	3.0
  M	(16,18)	=	-2.0
  M	(39,18)	=	-9.0
  M	(40,18)	=	6.0
  M	(43,18)	=	18.0
  M	(44,18)	=	-12.0
  M	(47,18)	=	-9.0
  M	(48,18)	=	6.0
  M	(55,18)	=	6.0
  M	(56,18)	=	-4.0
  M	(59,18)	=	-12.0
  M	(60,18)	=	8.0
  M	(63,18)	=	6.0
  M	(64,18)	=	-4.0
  M	(9,19)	=	-1.0
  M	(11,19)	=	3.0
  M	(12,19)	=	-2.0
  M	(13,19)	=	1.0
  M	(15,19)	=	-3.0
  M	(16,19)	=	2.0
  M	(41,19)	=	3.0
  M	(43,19)	=	-9.0
  M	(44,19)	=	6.0
  M	(45,19)	=	-3.0
  M	(47,19)	=	9.0
  M	(48,19)	=	-6.0
  M	(57,19)	=	-2.0
  M	(59,19)	=	6.0
  M	(60,19)	=	-4.0
  M	(61,19)	=	2.0
  M	(63,19)	=	-6.0
  M	(64,19)	=	4.0
  M	(11,20)	=	-3.0
  M	(12,20)	=	2.0
  M	(15,20)	=	3.0
  M	(16,20)	=	-2.0
  M	(43,20)	=	9.0
  M	(44,20)	=	-6.0
  M	(47,20)	=	-9.0
  M	(48,20)	=	6.0
  M	(59,20)	=	-6.0
  M	(60,20)	=	4.0
  M	(63,20)	=	6.0
  M	(64,20)	=	-4.0
  M	(37,21)	=	3.0
  M	(39,21)	=	-9.0
  M	(40,21)	=	6.0
  M	(41,21)	=	-6.0
  M	(43,21)	=	18.0
  M	(44,21)	=	-12.0
  M	(45,21)	=	3.0
  M	(47,21)	=	-9.0
  M	(48,21)	=	6.0
  M	(53,21)	=	-2.0
  M	(55,21)	=	6.0
  M	(56,21)	=	-4.0
  M	(57,21)	=	4.0
  M	(59,21)	=	-12.0
  M	(60,21)	=	8.0
  M	(61,21)	=	-2.0
  M	(63,21)	=	6.0
  M	(64,21)	=	-4.0
  M	(39,22)	=	9.0
  M	(40,22)	=	-6.0
  M	(43,22)	=	-18.0
  M	(44,22)	=	12.0
  M	(47,22)	=	9.0
  M	(48,22)	=	-6.0
  M	(55,22)	=	-6.0
  M	(56,22)	=	4.0
  M	(59,22)	=	12.0
  M	(60,22)	=	-8.0
  M	(63,22)	=	-6.0
  M	(64,22)	=	4.0
  M	(41,23)	=	-3.0
  M	(43,23)	=	9.0
  M	(44,23)	=	-6.0
  M	(45,23)	=	3.0
  M	(47,23)	=	-9.0
  M	(48,23)	=	6.0
  M	(57,23)	=	2.0
  M	(59,23)	=	-6.0
  M	(60,23)	=	4.0
  M	(61,23)	=	-2.0
  M	(63,23)	=	6.0
  M	(64,23)	=	-4.0
  M	(43,24)	=	-9.0
  M	(44,24)	=	6.0
  M	(47,24)	=	9.0
  M	(48,24)	=	-6.0
  M	(59,24)	=	6.0
  M	(60,24)	=	-4.0
  M	(63,24)	=	-6.0
  M	(64,24)	=	4.0
  M	(17,25)	=	1.0
  M	(19,25)	=	-3.0
  M	(20,25)	=	2.0
  M	(25,25)	=	-3.0
  M	(27,25)	=	9.0
  M	(28,25)	=	-6.0
  M	(29,25)	=	2.0
  M	(31,25)	=	-6.0
  M	(32,25)	=	4.0
  M	(33,25)	=	-2.0
  M	(35,25)	=	6.0
  M	(36,25)	=	-4.0
  M	(41,25)	=	6.0
  M	(43,25)	=	-18.0
  M	(44,25)	=	12.0
  M	(45,25)	=	-4.0
  M	(47,25)	=	12.0
  M	(48,25)	=	-8.0
  M	(49,25)	=	1.0
  M	(51,25)	=	-3.0
  M	(52,25)	=	2.0
  M	(57,25)	=	-3.0
  M	(59,25)	=	9.0
  M	(60,25)	=	-6.0
  M	(61,25)	=	2.0
  M	(63,25)	=	-6.0
  M	(64,25)	=	4.0
  M	(19,26)	=	3.0
  M	(20,26)	=	-2.0
  M	(27,26)	=	-9.0
  M	(28,26)	=	6.0
  M	(31,26)	=	6.0
  M	(32,26)	=	-4.0
  M	(35,26)	=	-6.0
  M	(36,26)	=	4.0
  M	(43,26)	=	18.0
  M	(44,26)	=	-12.0
  M	(47,26)	=	-12.0
  M	(48,26)	=	8.0
  M	(51,26)	=	3.0
  M	(52,26)	=	-2.0
  M	(59,26)	=	-9.0
  M	(60,26)	=	6.0
  M	(63,26)	=	6.0
  M	(64,26)	=	-4.0
  M	(25,27)	=	3.0
  M	(27,27)	=	-9.0
  M	(28,27)	=	6.0
  M	(29,27)	=	-2.0
  M	(31,27)	=	6.0
  M	(32,27)	=	-4.0
  M	(41,27)	=	-6.0
  M	(43,27)	=	18.0
  M	(44,27)	=	-12.0
  M	(45,27)	=	4.0
  M	(47,27)	=	-12.0
  M	(48,27)	=	8.0
  M	(57,27)	=	3.0
  M	(59,27)	=	-9.0
  M	(60,27)	=	6.0
  M	(61,27)	=	-2.0
  M	(63,27)	=	6.0
  M	(64,27)	=	-4.0
  M	(27,28)	=	9.0
  M	(28,28)	=	-6.0
  M	(31,28)	=	-6.0
  M	(32,28)	=	4.0
  M	(43,28)	=	-18.0
  M	(44,28)	=	12.0
  M	(47,28)	=	12.0
  M	(48,28)	=	-8.0
  M	(59,28)	=	9.0
  M	(60,28)	=	-6.0
  M	(63,28)	=	-6.0
  M	(64,28)	=	4.0
  M	(33,29)	=	-1.0
  M	(35,29)	=	3.0
  M	(36,29)	=	-2.0
  M	(41,29)	=	3.0
  M	(43,29)	=	-9.0
  M	(44,29)	=	6.0
  M	(45,29)	=	-2.0
  M	(47,29)	=	6.0
  M	(48,29)	=	-4.0
  M	(49,29)	=	1.0
  M	(51,29)	=	-3.0
  M	(52,29)	=	2.0
  M	(57,29)	=	-3.0
  M	(59,29)	=	9.0
  M	(60,29)	=	-6.0
  M	(61,29)	=	2.0
  M	(63,29)	=	-6.0
  M	(64,29)	=	4.0
  M	(35,30)	=	-3.0
  M	(36,30)	=	2.0
  M	(43,30)	=	9.0
  M	(44,30)	=	-6.0
  M	(47,30)	=	-6.0
  M	(48,30)	=	4.0
  M	(51,30)	=	3.0
  M	(52,30)	=	-2.0
  M	(59,30)	=	-9.0
  M	(60,30)	=	6.0
  M	(63,30)	=	6.0
  M	(64,30)	=	-4.0
  M	(41,31)	=	-3.0
  M	(43,31)	=	9.0
  M	(44,31)	=	-6.0
  M	(45,31)	=	2.0
  M	(47,31)	=	-6.0
  M	(48,31)	=	4.0
  M	(57,31)	=	3.0
  M	(59,31)	=	-9.0
  M	(60,31)	=	6.0
  M	(61,31)	=	-2.0
  M	(63,31)	=	6.0
  M	(64,31)	=	-4.0
  M	(43,32)	=	-9.0
  M	(44,32)	=	6.0
  M	(47,32)	=	6.0
  M	(48,32)	=	-4.0
  M	(59,32)	=	9.0
  M	(60,32)	=	-6.0
  M	(63,32)	=	-6.0
  M	(64,32)	=	4.0
  M	(6,33)	=	1.0
  M	(7,33)	=	-2.0
  M	(8,33)	=	1.0
  M	(10,33)	=	-2.0
  M	(11,33)	=	4.0
  M	(12,33)	=	-2.0
  M	(14,33)	=	1.0
  M	(15,33)	=	-2.0
  M	(16,33)	=	1.0
  M	(38,33)	=	-3.0
  M	(39,33)	=	6.0
  M	(40,33)	=	-3.0
  M	(42,33)	=	6.0
  M	(43,33)	=	-12.0
  M	(44,33)	=	6.0
  M	(46,33)	=	-3.0
  M	(47,33)	=	6.0
  M	(48,33)	=	-3.0
  M	(54,33)	=	2.0
  M	(55,33)	=	-4.0
  M	(56,33)	=	2.0
  M	(58,33)	=	-4.0
  M	(59,33)	=	8.0
  M	(60,33)	=	-4.0
  M	(62,33)	=	2.0
  M	(63,33)	=	-4.0
  M	(64,33)	=	2.0
  M	(7,34)	=	-1.0
  M	(8,34)	=	1.0
  M	(11,34)	=	2.0
  M	(12,34)	=	-2.0
  M	(15,34)	=	-1.0
  M	(16,34)	=	1.0
  M	(39,34)	=	3.0
  M	(40,34)	=	-3.0
  M	(43,34)	=	-6.0
  M	(44,34)	=	6.0
  M	(47,34)	=	3.0
  M	(48,34)	=	-3.0
  M	(55,34)	=	-2.0
  M	(56,34)	=	2.0
  M	(59,34)	=	4.0
  M	(60,34)	=	-4.0
  M	(63,34)	=	-2.0
  M	(64,34)	=	2.0
  M	(10,35)	=	-1.0
  M	(11,35)	=	2.0
  M	(12,35)	=	-1.0
  M	(14,35)	=	1.0
  M	(15,35)	=	-2.0
  M	(16,35)	=	1.0
  M	(42,35)	=	3.0
  M	(43,35)	=	-6.0
  M	(44,35)	=	3.0
  M	(46,35)	=	-3.0
  M	(47,35)	=	6.0
  M	(48,35)	=	-3.0
  M	(58,35)	=	-2.0
  M	(59,35)	=	4.0
  M	(60,35)	=	-2.0
  M	(62,35)	=	2.0
  M	(63,35)	=	-4.0
  M	(64,35)	=	2.0
  M	(11,36)	=	1.0
  M	(12,36)	=	-1.0
  M	(15,36)	=	-1.0
  M	(16,36)	=	1.0
  M	(43,36)	=	-3.0
  M	(44,36)	=	3.0
  M	(47,36)	=	3.0
  M	(48,36)	=	-3.0
  M	(59,36)	=	2.0
  M	(60,36)	=	-2.0
  M	(63,36)	=	-2.0
  M	(64,36)	=	2.0
  M	(38,37)	=	3.0
  M	(39,37)	=	-6.0
  M	(40,37)	=	3.0
  M	(42,37)	=	-6.0
  M	(43,37)	=	12.0
  M	(44,37)	=	-6.0
  M	(46,37)	=	3.0
  M	(47,37)	=	-6.0
  M	(48,37)	=	3.0
  M	(54,37)	=	-2.0
  M	(55,37)	=	4.0
  M	(56,37)	=	-2.0
  M	(58,37)	=	4.0
  M	(59,37)	=	-8.0
  M	(60,37)	=	4.0
  M	(62,37)	=	-2.0
  M	(63,37)	=	4.0
  M	(64,37)	=	-2.0
  M	(39,38)	=	-3.0
  M	(40,38)	=	3.0
  M	(43,38)	=	6.0
  M	(44,38)	=	-6.0
  M	(47,38)	=	-3.0
  M	(48,38)	=	3.0
  M	(55,38)	=	2.0
  M	(56,38)	=	-2.0
  M	(59,38)	=	-4.0
  M	(60,38)	=	4.0
  M	(63,38)	=	2.0
  M	(64,38)	=	-2.0
  M	(42,39)	=	-3.0
  M	(43,39)	=	6.0
  M	(44,39)	=	-3.0
  M	(46,39)	=	3.0
  M	(47,39)	=	-6.0
  M	(48,39)	=	3.0
  M	(58,39)	=	2.0
  M	(59,39)	=	-4.0
  M	(60,39)	=	2.0
  M	(62,39)	=	-2.0
  M	(63,39)	=	4.0
  M	(64,39)	=	-2.0
  M	(43,40)	=	3.0
  M	(44,40)	=	-3.0
  M	(47,40)	=	-3.0
  M	(48,40)	=	3.0
  M	(59,40)	=	-2.0
  M	(60,40)	=	2.0
  M	(63,40)	=	2.0
  M	(64,40)	=	-2.0
  M	(18,41)	=	1.0
  M	(19,41)	=	-2.0
  M	(20,41)	=	1.0
  M	(26,41)	=	-3.0
  M	(27,41)	=	6.0
  M	(28,41)	=	-3.0
  M	(30,41)	=	2.0
  M	(31,41)	=	-4.0
  M	(32,41)	=	2.0
  M	(34,41)	=	-2.0
  M	(35,41)	=	4.0
  M	(36,41)	=	-2.0
  M	(42,41)	=	6.0
  M	(43,41)	=	-12.0
  M	(44,41)	=	6.0
  M	(46,41)	=	-4.0
  M	(47,41)	=	8.0
  M	(48,41)	=	-4.0
  M	(50,41)	=	1.0
  M	(51,41)	=	-2.0
  M	(52,41)	=	1.0
  M	(58,41)	=	-3.0
  M	(59,41)	=	6.0
  M	(60,41)	=	-3.0
  M	(62,41)	=	2.0
  M	(63,41)	=	-4.0
  M	(64,41)	=	2.0
  M	(19,42)	=	-1.0
  M	(20,42)	=	1.0
  M	(27,42)	=	3.0
  M	(28,42)	=	-3.0
  M	(31,42)	=	-2.0
  M	(32,42)	=	2.0
  M	(35,42)	=	2.0
  M	(36,42)	=	-2.0
  M	(43,42)	=	-6.0
  M	(44,42)	=	6.0
  M	(47,42)	=	4.0
  M	(48,42)	=	-4.0
  M	(51,42)	=	-1.0
  M	(52,42)	=	1.0
  M	(59,42)	=	3.0
  M	(60,42)	=	-3.0
  M	(63,42)	=	-2.0
  M	(64,42)	=	2.0
  M	(26,43)	=	3.0
  M	(27,43)	=	-6.0
  M	(28,43)	=	3.0
  M	(30,43)	=	-2.0
  M	(31,43)	=	4.0
  M	(32,43)	=	-2.0
  M	(42,43)	=	-6.0
  M	(43,43)	=	12.0
  M	(44,43)	=	-6.0
  M	(46,43)	=	4.0
  M	(47,43)	=	-8.0
  M	(48,43)	=	4.0
  M	(58,43)	=	3.0
  M	(59,43)	=	-6.0
  M	(60,43)	=	3.0
  M	(62,43)	=	-2.0
  M	(63,43)	=	4.0
  M	(64,43)	=	-2.0
  M	(27,44)	=	-3.0
  M	(28,44)	=	3.0
  M	(31,44)	=	2.0
  M	(32,44)	=	-2.0
  M	(43,44)	=	6.0
  M	(44,44)	=	-6.0
  M	(47,44)	=	-4.0
  M	(48,44)	=	4.0
  M	(59,44)	=	-3.0
  M	(60,44)	=	3.0
  M	(63,44)	=	2.0
  M	(64,44)	=	-2.0
  M	(34,45)	=	-1.0
  M	(35,45)	=	2.0
  M	(36,45)	=	-1.0
  M	(42,45)	=	3.0
  M	(43,45)	=	-6.0
  M	(44,45)	=	3.0
  M	(46,45)	=	-2.0
  M	(47,45)	=	4.0
  M	(48,45)	=	-2.0
  M	(50,45)	=	1.0
  M	(51,45)	=	-2.0
  M	(52,45)	=	1.0
  M	(58,45)	=	-3.0
  M	(59,45)	=	6.0
  M	(60,45)	=	-3.0
  M	(62,45)	=	2.0
  M	(63,45)	=	-4.0
  M	(64,45)	=	2.0
  M	(35,46)	=	1.0
  M	(36,46)	=	-1.0
  M	(43,46)	=	-3.0
  M	(44,46)	=	3.0
  M	(47,46)	=	2.0
  M	(48,46)	=	-2.0
  M	(51,46)	=	-1.0
  M	(52,46)	=	1.0
  M	(59,46)	=	3.0
  M	(60,46)	=	-3.0
  M	(63,46)	=	-2.0
  M	(64,46)	=	2.0
  M	(42,47)	=	-3.0
  M	(43,47)	=	6.0
  M	(44,47)	=	-3.0
  M	(46,47)	=	2.0
  M	(47,47)	=	-4.0
  M	(48,47)	=	2.0
  M	(58,47)	=	3.0
  M	(59,47)	=	-6.0
  M	(60,47)	=	3.0
  M	(62,47)	=	-2.0
  M	(63,47)	=	4.0
  M	(64,47)	=	-2.0
  M	(43,48)	=	3.0
  M	(44,48)	=	-3.0
  M	(47,48)	=	-2.0
  M	(48,48)	=	2.0
  M	(59,48)	=	-3.0
  M	(60,48)	=	3.0
  M	(63,48)	=	2.0
  M	(64,48)	=	-2.0
  M	(21,49)	=	1.0
  M	(23,49)	=	-3.0
  M	(24,49)	=	2.0
  M	(25,49)	=	-2.0
  M	(27,49)	=	6.0
  M	(28,49)	=	-4.0
  M	(29,49)	=	1.0
  M	(31,49)	=	-3.0
  M	(32,49)	=	2.0
  M	(37,49)	=	-2.0
  M	(39,49)	=	6.0
  M	(40,49)	=	-4.0
  M	(41,49)	=	4.0
  M	(43,49)	=	-12.0
  M	(44,49)	=	8.0
  M	(45,49)	=	-2.0
  M	(47,49)	=	6.0
  M	(48,49)	=	-4.0
  M	(53,49)	=	1.0
  M	(55,49)	=	-3.0
  M	(56,49)	=	2.0
  M	(57,49)	=	-2.0
  M	(59,49)	=	6.0
  M	(60,49)	=	-4.0
  M	(61,49)	=	1.0
  M	(63,49)	=	-3.0
  M	(64,49)	=	2.0
  M	(23,50)	=	3.0
  M	(24,50)	=	-2.0
  M	(27,50)	=	-6.0
  M	(28,50)	=	4.0
  M	(31,50)	=	3.0
  M	(32,50)	=	-2.0
  M	(39,50)	=	-6.0
  M	(40,50)	=	4.0
  M	(43,50)	=	12.0
  M	(44,50)	=	-8.0
  M	(47,50)	=	-6.0
  M	(48,50)	=	4.0
  M	(55,50)	=	3.0
  M	(56,50)	=	-2.0
  M	(59,50)	=	-6.0
  M	(60,50)	=	4.0
  M	(63,50)	=	3.0
  M	(64,50)	=	-2.0
  M	(25,51)	=	-1.0
  M	(27,51)	=	3.0
  M	(28,51)	=	-2.0
  M	(29,51)	=	1.0
  M	(31,51)	=	-3.0
  M	(32,51)	=	2.0
  M	(41,51)	=	2.0
  M	(43,51)	=	-6.0
  M	(44,51)	=	4.0
  M	(45,51)	=	-2.0
  M	(47,51)	=	6.0
  M	(48,51)	=	-4.0
  M	(57,51)	=	-1.0
  M	(59,51)	=	3.0
  M	(60,51)	=	-2.0
  M	(61,51)	=	1.0
  M	(63,51)	=	-3.0
  M	(64,51)	=	2.0
  M	(27,52)	=	-3.0
  M	(28,52)	=	2.0
  M	(31,52)	=	3.0
  M	(32,52)	=	-2.0
  M	(43,52)	=	6.0
  M	(44,52)	=	-4.0
  M	(47,52)	=	-6.0
  M	(48,52)	=	4.0
  M	(59,52)	=	-3.0
  M	(60,52)	=	2.0
  M	(63,52)	=	3.0
  M	(64,52)	=	-2.0
  M	(37,53)	=	-1.0
  M	(39,53)	=	3.0
  M	(40,53)	=	-2.0
  M	(41,53)	=	2.0
  M	(43,53)	=	-6.0
  M	(44,53)	=	4.0
  M	(45,53)	=	-1.0
  M	(47,53)	=	3.0
  M	(48,53)	=	-2.0
  M	(53,53)	=	1.0
  M	(55,53)	=	-3.0
  M	(56,53)	=	2.0
  M	(57,53)	=	-2.0
  M	(59,53)	=	6.0
  M	(60,53)	=	-4.0
  M	(61,53)	=	1.0
  M	(63,53)	=	-3.0
  M	(64,53)	=	2.0
  M	(39,54)	=	-3.0
  M	(40,54)	=	2.0
  M	(43,54)	=	6.0
  M	(44,54)	=	-4.0
  M	(47,54)	=	-3.0
  M	(48,54)	=	2.0
  M	(55,54)	=	3.0
  M	(56,54)	=	-2.0
  M	(59,54)	=	-6.0
  M	(60,54)	=	4.0
  M	(63,54)	=	3.0
  M	(64,54)	=	-2.0
  M	(41,55)	=	1.0
  M	(43,55)	=	-3.0
  M	(44,55)	=	2.0
  M	(45,55)	=	-1.0
  M	(47,55)	=	3.0
  M	(48,55)	=	-2.0
  M	(57,55)	=	-1.0
  M	(59,55)	=	3.0
  M	(60,55)	=	-2.0
  M	(61,55)	=	1.0
  M	(63,55)	=	-3.0
  M	(64,55)	=	2.0
  M	(43,56)	=	3.0
  M	(44,56)	=	-2.0
  M	(47,56)	=	-3.0
  M	(48,56)	=	2.0
  M	(59,56)	=	-3.0
  M	(60,56)	=	2.0
  M	(63,56)	=	3.0
  M	(64,56)	=	-2.0
  M	(22,57)	=	1.0
  M	(23,57)	=	-2.0
  M	(24,57)	=	1.0
  M	(26,57)	=	-2.0
  M	(27,57)	=	4.0
  M	(28,57)	=	-2.0
  M	(30,57)	=	1.0
  M	(31,57)	=	-2.0
  M	(32,57)	=	1.0
  M	(38,57)	=	-2.0
  M	(39,57)	=	4.0
  M	(40,57)	=	-2.0
  M	(42,57)	=	4.0
  M	(43,57)	=	-8.0
  M	(44,57)	=	4.0
  M	(46,57)	=	-2.0
  M	(47,57)	=	4.0
  M	(48,57)	=	-2.0
  M	(54,57)	=	1.0
  M	(55,57)	=	-2.0
  M	(56,57)	=	1.0
  M	(58,57)	=	-2.0
  M	(59,57)	=	4.0
  M	(60,57)	=	-2.0
  M	(62,57)	=	1.0
  M	(63,57)	=	-2.0
  M	(64,57)	=	1.0
  M	(23,58)	=	-1.0
  M	(24,58)	=	1.0
  M	(27,58)	=	2.0
  M	(28,58)	=	-2.0
  M	(31,58)	=	-1.0
  M	(32,58)	=	1.0
  M	(39,58)	=	2.0
  M	(40,58)	=	-2.0
  M	(43,58)	=	-4.0
  M	(44,58)	=	4.0
  M	(47,58)	=	2.0
  M	(48,58)	=	-2.0
  M	(55,58)	=	-1.0
  M	(56,58)	=	1.0
  M	(59,58)	=	2.0
  M	(60,58)	=	-2.0
  M	(63,58)	=	-1.0
  M	(64,58)	=	1.0
  M	(26,59)	=	-1.0
  M	(27,59)	=	2.0
  M	(28,59)	=	-1.0
  M	(30,59)	=	1.0
  M	(31,59)	=	-2.0
  M	(32,59)	=	1.0
  M	(42,59)	=	2.0
  M	(43,59)	=	-4.0
  M	(44,59)	=	2.0
  M	(46,59)	=	-2.0
  M	(47,59)	=	4.0
  M	(48,59)	=	-2.0
  M	(58,59)	=	-1.0
  M	(59,59)	=	2.0
  M	(60,59)	=	-1.0
  M	(62,59)	=	1.0
  M	(63,59)	=	-2.0
  M	(64,59)	=	1.0
  M	(27,60)	=	1.0
  M	(28,60)	=	-1.0
  M	(31,60)	=	-1.0
  M	(32,60)	=	1.0
  M	(43,60)	=	-2.0
  M	(44,60)	=	2.0
  M	(47,60)	=	2.0
  M	(48,60)	=	-2.0
  M	(59,60)	=	1.0
  M	(60,60)	=	-1.0
  M	(63,60)	=	-1.0
  M	(64,60)	=	1.0
  M	(38,61)	=	-1.0
  M	(39,61)	=	2.0
  M	(40,61)	=	-1.0
  M	(42,61)	=	2.0
  M	(43,61)	=	-4.0
  M	(44,61)	=	2.0
  M	(46,61)	=	-1.0
  M	(47,61)	=	2.0
  M	(48,61)	=	-1.0
  M	(54,61)	=	1.0
  M	(55,61)	=	-2.0
  M	(56,61)	=	1.0
  M	(58,61)	=	-2.0
  M	(59,61)	=	4.0
  M	(60,61)	=	-2.0
  M	(62,61)	=	1.0
  M	(63,61)	=	-2.0
  M	(64,61)	=	1.0
  M	(39,62)	=	1.0
  M	(40,62)	=	-1.0
  M	(43,62)	=	-2.0
  M	(44,62)	=	2.0
  M	(47,62)	=	1.0
  M	(48,62)	=	-1.0
  M	(55,62)	=	-1.0
  M	(56,62)	=	1.0
  M	(59,62)	=	2.0
  M	(60,62)	=	-2.0
  M	(63,62)	=	-1.0
  M	(64,62)	=	1.0
  M	(42,63)	=	1.0
  M	(43,63)	=	-2.0
  M	(44,63)	=	1.0
  M	(46,63)	=	-1.0
  M	(47,63)	=	2.0
  M	(48,63)	=	-1.0
  M	(58,63)	=	-1.0
  M	(59,63)	=	2.0
  M	(60,63)	=	-1.0
  M	(62,63)	=	1.0
  M	(63,63)	=	-2.0
  M	(64,63)	=	1.0
  M	(43,64)	=	-1.0
  M	(44,64)	=	1.0
  M	(47,64)	=	1.0
  M	(48,64)	=	-1.0
  M	(59,64)	=	1.0
  M	(60,64)	=	-1.0
  M	(63,64)	=	-1.0
  M	(64,64)	=	1.0

end subroutine Minv

subroutine interp(input,output,nx,ny,nz)

  implicit none

  ! Input/Output
  integer :: nx,ny,nz ! These are the GRID dimensions:
  real(kind=8), dimension(nx,ny,nz) :: input
  real(kind=8), dimension(2*nx-1,2*ny-1,2*nz-1) :: output
  real(kind=8) :: corners(8),s(3),shp(8)
  ! Working
  integer :: i,j,k,ii,jj,kk

  ! The purpose of this routine is to interpolate from the coarse
  ! "input" to the refined "output"

  ! Loop over each CELL:
  do k=1,nz-1
     do j=1,ny-1
        do i=1,nx-1
           ! Here we will produce the new 27 points for this cell in
           ! lieu of the former 8. This does more work than necessary
           ! currently, but we're going for simplicity not speed.

           corners(1) = input(i  ,j  ,k  )
           corners(2) = input(i+1,j  ,k  )
           corners(3) = input(i  ,j+1,k  )
           corners(4) = input(i+1,j+1,k  )
           corners(5) = input(i  ,j  ,k+1)
           corners(6) = input(i+1,j  ,k+1)
           corners(7) = input(i  ,j+1,k+1)
           corners(8) = input(i+1,j+1,k+1)

           do kk=0,2
              do jj=0,2
                 do ii=0,2
                    s = (/dble(ii)/2.0,dble(jj)/2.0,dble(kk)/2.0/)

                    call shape_function(s,shp)

                    output((i-1)*2+ii+1,(j-1)*2+jj+1,(k-1)*2+kk+1) = &
                         dot_product(shp,corners)
                 end do
              end do
           end do
        end do
     end do
  end do


end subroutine interp

subroutine shape_function(s,shp)
 
  real(kind=8) :: s(3),shp(8)
  real(kind=8) ::  nns(2),nnr(2),nnt(2)


  nnr(1) = 1.0 - s(1)
  nns(1) = 1.0 - s(2)
  nnt(1) = 1.0 - s(3)

  shp(1) = nnr(1)*nns(1)*nnt(1)
  shp(2) = s(1)*nns(1)*nnt(1)
  shp(3) = nnr(1)*s(2)*nnt(1)
  shp(4) = s(1)*s(2)*nnt(1)
  shp(5) = nnr(1)*nns(1)*s(3)
  shp(6) = s(1)*nns(1)*s(3)
  shp(7) = nnr(1)*s(2)*s(3)
  shp(8) = s(1)*s(2)*s(3)

end subroutine shape_function

subroutine splineDERIV(X,s,dXds,n)

  implicit none

  ! Input/Output
  integer ,intent(in) :: n
  real(kind=8)   ,intent(in) :: X(n)
  real(kind=8)   ,intent(in) :: s(n)
  real(kind=8)   ,intent(out) :: dXds(n)

  ! Local
  integer            :: i
  
  ! Functions
  real(kind=8)  :: alpha
  alpha(i) = (s(i)-s(i-1))/(s(i+1) - s(i-1))

  dXds(1) = 2*(X(2)-X(1)) - ( (1-alpha(2))*(X(2)-X(1)) + alpha(2)*(X(3)-X(2))*(s(2)-s(1))/(s(3)-s(2)))
  dXds(n) = 2*(X(n)-X(n-1)) - (  (1-alpha(n-1))*(X(n-1)-X(n-2))*(s(n)-s(n-1))/(s(n-1)-s(n-2)) +alpha(n-1)*(X(n)-X(n-1)))

  do i=2,n-1
     dXds(i) = alpha(i)*(X(i+1)-X(i))+ (1-alpha(i))*(X(i)-X(i-1))
  end do

!dxds = 0
end subroutine splineDERIV

subroutine splineEVAL(s,coef,values)

  implicit none
  ! Input
  real(kind=8),intent(in) :: s(3)
  real(kind=8),intent(in) :: coef(3,64)
  real(kind=8)            :: values(3)

  ! Local
  integer :: i,j,k
  real(kind=8) :: pows(4,3)

 !  ! No loop for efficiency
  
  pows(1,1) = 1.0
  pows(2,1) = s(1)
  pows(3,1) = s(1)*s(1)
  pows(4,1) = s(1)*s(1)*s(1)

  pows(1,2) = 1.0
  pows(2,2) = s(2)
  pows(3,2) = s(2)*s(2)
  pows(4,2) = s(2)*s(2)*s(2)

  pows(1,3) = 1.0
  pows(2,3) = s(3)
  pows(3,3) = s(3)*s(3)
  pows(4,3) = s(3)*s(3)*s(3)

  values = 0.0

  do k=1,4
     do j=1,4
        do i=1,4
           values = values + coef(:,i+4*j+k*16-20)* pows(i,1)*pows(j,2)*pows(k,3)
        end do
     end do
  end do

end subroutine splineEVAL

subroutine para3d(X,n,m,l,ndim,S)
  !***DESCRIPTION
  !
  !     Written by Gaetan Kenway
  !
  !     Abstract para3d calculates the parametric locations for a 3d block
  !
  !     Description of Arguments
  !     Input
  !     X       - Real,size(3,n,m,l): Coordiantes
  !     Output
  !     S       - Real,size(3,n,m,l): The u,v,w parametric positions

  implicit none

  ! Input
  integer   , intent(in)   :: n,m,l,ndim
  real(kind=8), intent(in) :: X(ndim,n,m,l)


  ! Output
  real(kind=8),  intent(out)  :: S(ndim,n,m,l)

  ! Working 
  integer :: i,j,k
  real(kind=8)   :: eps = 1.0e-10
  real(kind=8) :: DELI,DELJ,DELK

  DELI(I,J,K) = sqrt ((X(1,I,J,K) - X(1,I-1,J,K)) ** 2 + &
       (X(2,I,J,K) - X(2,I-1,J,K)) ** 2 + &
       (X(3,I,J,K) - X(3,I-1,J,K)) ** 2)

  DELJ(I,J,K) = sqrt ((X(1,I,J,K) - X(1,I,J-1,K)) ** 2 + &
       (X(2,I,J,K) - X(2,I,J-1,K)) ** 2 + &
       (X(3,I,J,K) - X(3,I,J-1,K)) ** 2)

  DELK(I,J,K) = sqrt ((X(1,I,J,K) - X(1,I,J,K-1)) ** 2 + &
       (X(2,I,J,K) - X(2,I,J,K-1)) ** 2 + &
       (X(3,I,J,K) - X(3,I,J,K-1)) ** 2)


  ! Zero the three low-end faces (or edges if one plane is specified).

  do K = 1, l
     do J = 1, m
        S(1,1,J,K) = 0.0
     end do

     do I = 1, n
        S(2,I,1,K) = 0.0
     end do
  end do

  do J = 1, m
     do I = 1, n
        S(3,I,J,1) = 0.0
     end do
  end do

  !     Set up the low-end edge lines because they are missed by the
  !     following loops over most of the low-end faces:

  do I = 2, n
     S(1,I,1,1) = S(1,I-1,1,1) + DELI(I,1,1)
  end do

  do J = 2, m
     S(2,1,J,1) = S(2,1,J-1,1) + DELJ(1,J,1)
  end do

  do K = 2, l
     S(3,1,1,K) = S(3,1,1,K-1) + DELK(1,1,K)
  end do

  !     Set up the rest of the low-end face lines because they are
  !     missed by the the main loop over most of the volume.

  do K = 2, l
     do J = 2, m
        S(2,1,J,K) = S(2,1,J-1,K) + DELJ(1,J,K)
        S(3,1,J,K) = S(3,1,J,K-1) + DELK(1,J,K)
     end do
     do I = 2, n
        S(1,I,1,K) = S(1,I-1,1,K) + DELI(I,1,K)
        S(3,I,1,K) = S(3,I,1,K-1) + DELK(I,1,K)
     end do

  end do

  do J = 2, m
     do I = 2, n
        S(1,I,J,1) = S(1,I-1,J,1) + DELI(I,J,1)
        S(2,I,J,1) = S(2,I,J-1,1) + DELJ(I,J,1)
     end do
  end do

  !     Traverse the block just once for all lines except those within
  !     the low-end faces.

  do K = 2, l
     do J = 2, m
        do I = 2, n
           S(1,I,J,K) = S(1,I-1,J,K) + DELI(I,J,K)
           S(2,I,J,K) = S(2,I,J-1,K) + DELJ(I,J,K)
           S(3,I,J,K) = S(3,I,J,K-1) + DELK(I,J,K)
        end do
     end do
  end do

  !     Normalizing requires another pass through the volume.
  !     Handle lines of zero length first by inserting uniform
  !     distributions.  Then the standard normalization can be
  !     applied safely everywhere.

  do K = 1, l

     !        Zero-length lines in the I direction?

     do J = 1, m
        if (S(1,n,J,K) < eps) then
           do I = 2, n
              S(1,I,J,K) = I - 1
           end do
        end if
     end do

     !        Zero-length lines in the J direction?

     do I = 1, n
        if (S(2,I,m,K) < eps) then
           do J = 2, m
              S(2,I,J,K) = J - 1
           end do
        end if
     end do
  end do

  !     Zero-length lines in the K direction?

  do J = 1, m
     do I = 1, n
        if (S(3,I,J,l) < eps) then
           do K = 2, l
              S(3,I,J,K) = K - 1
           end do
        end if
     end do
  end do

  !     Normalize:

  do K = 1, l
     do J = 1, m
        do I = 1, n
           S(1,I,J,K) = S(1,I,J,K) / S(1,n,J,K)
           S(2,I,J,K) = S(2,I,J,K) / S(2,I,m,K)
           S(3,I,J,K) = S(3,I,J,K) / S(3,I,J,l)
        end do
     end do
  end do

  !     Finally, precise 1s for the three high-end faces:

  do K = 1, l
     do J = 1, m
        S(1,n,J,K) = 1.0
     end do

     do I = 1, n
        S(2,I,m,K) = 1.0
     end do
  end do

  do J = 1, m
     do I = 1, n
        S(3,I,J,l) = 1.0
     end do
  end do

end subroutine para3d
