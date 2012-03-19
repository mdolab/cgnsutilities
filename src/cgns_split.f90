module cgnsGrid
  implicit none
  save
  ! Maximum number of splits per edge...exceeding this maximum is not
  ! checked so make sure it is large enough!
  integer, parameter :: MAX_SPLITS = 100
  integer, parameter :: MAX_SUB_BLOCKS = 100

  type cgns1to1Conn
     ! Represents a 1to1 point connection on original grid
     character*32 :: connectName
     character*32 :: donorName
     integer :: range(6)
     integer :: donor_range(6)
     integer :: transform(3)
     integer :: faceID
  end type cgns1to1Conn
  
  type cgnsBC
     ! Represents a Boundary Condition on original grid
     integer :: range(6)
     integer :: faceID
     character*32 :: familyName
  end type cgnsBC

  type cgnsBlock

     ! Represents an original CGNS block in unsplit grid
     integer :: il,jl,kl
     integer :: iSplit(MAX_SPLITS)
     integer :: jSplit(MAX_SPLITS)
     integer :: kSplit(MAX_SPLITS)

     integer :: nISplit, nJSplit, nKSplit
     integer :: n1to1, nbcs
     character*32 :: zoneName
     type(cgns1to1Conn) , dimension(:) , allocatable :: b2bConn
     type(cgnsBC) , dimension(:) , allocatable :: bcs

  end type cgnsBlock


end module cgnsGrid

program cgns_split

  use cgnsGrid

  implicit none
  include 'cgnslib_f.h'

  type(cgnsBlock), dimension(:),allocatable :: blocks
  integer  CellDim, PhysDim
  integer nbases, nzones, zonesize(9)

  integer ier, zonetype
  integer n,nn,mm,i,j,k,idim

  integer cg_in, cg_out, base, zone, zoneCounter

  integer :: ptSetType, normalDataType
  integer :: sizes(9),coordID,blockStart(3),blockEnd(3)

  integer donor_range(6),transform(3)
  integer nbocos,n1to1,bocotype
  integer NormalIndex(3), NormalListFlag, ndataset,datatype
  integer ptset_type, npnts, points(6)
  integer n_user_split
  integer, allocatable, dimension(:,:) ::  user_splits
  character*32 in_filename,out_filename, split_file
  character*32 basename, zonename,boconame
  character*32 connectname,donorname, famname

  double precision data_double(6)
  double precision, allocatable, dimension(:,:,:,:) :: tempx
  double precision :: time(3)

  N = IARGC ()
  print *,'N:',N
  if ((N .ne. 2) .and. (N .ne. 3)) then
     print *,'Error: cgns_split must be called with the following arguments:'
     print *,'./cgns_split cgns_in_file.cgns cgns_out_file.cgns [split_file]'
     print *,'The optional split file contains information regarding additional'
     print *,'splits the user may want. The syntax is as follows:'
     print *,'3'
     print *,'1  1  17'
     print *,'5  2   9'
     print *,'18 3  33'
     print *,'The first number, n=3, is the number of splits. The remaining n lines contain'
     print *,'A zone number, the direction (1=i, 2=j, 3=k) and the index on which to perform'
     print *,'the split'
     stop
  end if

  call getarg(1, in_filename)
  call getarg(2, out_filename)

  if (N == 3) then
     call getArg(3, split_file)
     ! Read split file:
     open(11, FILE=split_file,status='old')
     read(11,*) n_user_split
     allocate(user_splits(n_user_split,3))
     print *,'nsplit:',n_user_split
     do i=1,n_user_split
        read(11,*) user_splits(i,1),user_splits(i,2),user_splits(i,3)
     end do

     ! Write out a check to sceen:
     print *,'The additional user-supplied splits will be used:'
     print *,'Block ID    Index    Direction'
     do i=1,n_user_split
        write(*,10) user_splits(i,1), user_splits(i,2), user_splits(i,3)
     end do

10   format(I9, I9, I5)
  else
     n_user_split = 0
     allocate(user_splits(n_user_split,3))
  end if

  call cg_open_f(in_filename,CG_MODE_READ, cg_in, ier)

  call cpu_time(time(1))
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

  ! Allocate memory for the array of blocks
  allocate(blocks(nzones))

  ! This first loop reads the bcs and the 1to1 block connections that
  ! will be used  for the basis of hte splits
  
  do nn=1,nzones
     call cg_zone_read_f(cg_in, base, nn, zonename, zonesize, ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     blocks(nn)%il = zonesize(1)
     blocks(nn)%jl = zonesize(2)
     blocks(nn)%kl = zonesize(3)
     blocks(nn)%nISplit = 0
     blocks(nn)%nJSplit = 0
     blocks(nn)%nKSplit = 0
     blocks(nn)%zoneName = zoneName
     ! Add the corners as "required" splits
     call addReqSplit(nn,(/1,1,1,zonesize(1),zonesize(2),zonesize(3)/))
     
     ! Determine number of bcs for this zone
     call cg_nbocos_f(cg_in, base, nn, nBocos,ier)
     if (ier .eq. CG_ERROR) call cg_error_exit_f

     allocate(blocks(nn)%bcs(nBocos))
     blocks(nn)%nbcs = nBocos
     do mm=1,nBocos
        ! Get Boundary Condition Info
        call cg_boco_info_f(cg_in, base, nn, mm , boconame,bocotype,&
             ptset_type,npnts,NormalIndex,NormalListFlag,datatype,ndataset,ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f

        call cg_boco_read_f(cg_in, base, nn, mm, points,data_double, ier)
        if (ier .eq. CG_ERROR) call cg_error_exit_f


        blocks(nn)%bcs(mm)%range = points
        if (points(1) == points(4) .and. points(1) == 1) then
           blocks(nn)%bcs(mm)%faceID = 1
        else if (points(1) == points(4) .and. points(1) == blocks(nn)%il) then
           blocks(nn)%bcs(mm)%faceID = 2
        else if (points(2) == points(5) .and. points(2) == 1) then
           blocks(nn)%bcs(mm)%faceID = 3
        else if (points(2) == points(5) .and. points(2) == blocks(nn)%jl) then
           blocks(nn)%bcs(mm)%faceID = 4
        else if (points(3) == points(6) .and. points(3) == 1) then
           blocks(nn)%bcs(mm)%faceID = 5
        else if (points(3) == points(6) .and. points(3) == blocks(nn)%kl) then
           blocks(nn)%bcs(mm)%faceID = 6
        end if

        ! Get the Family Info
        call cg_goto_f(cg_in,base,ier,"Zone_t",nn,"ZoneBC_t",1,"BC_t",mm,"end")
        if (ier == 0) then ! Node exits
           call cg_famname_read_f(famName, ier)
           blocks(nn)%bcs(mm)%familyName = famName
        else
           blocks(nn)%bcs(mm)%familyName = "wall"
        end if

     end do

     ! Determine number of 1to1 Conns for this zone
     call cg_n1to1_f(cg_in, base, nn, n1to1, ier)
     allocate(blocks(nn)%b2bConn(n1to1))
     blocks(nn)%n1to1 = n1to1
     do mm=1,n1to1
        call cg_1to1_read_f(cg_in, base, nn, mm, connectName, donorname, points,&
             donor_range, transform, ier)
        blocks(nn)%b2bConn(mm)%connectName = connectName
        blocks(nn)%b2bConn(mm)%donorName = donorName
        blocks(nn)%b2bConn(mm)%range = points
        blocks(nn)%b2bConn(mm)%donor_range = donor_range
        blocks(nn)%b2bConn(mm)%transform = transform
        
        if (points(1) == points(4) .and. points(1) == 1) then
           blocks(nn)%b2bconn(mm)%faceID = 1
        else if (points(1) == points(4) .and. points(1) == blocks(nn)%il) then
           blocks(nn)%b2bconn(mm)%faceID = 2
        else if (points(2) == points(5) .and. points(2) == 1) then
           blocks(nn)%b2bconn(mm)%faceID = 3
        else if (points(2) == points(5) .and. points(2) == blocks(nn)%jl) then
           blocks(nn)%b2bconn(mm)%faceID = 4
        else if (points(3) == points(6) .and. points(3) == 1) then
           blocks(nn)%b2bconn(mm)%faceID = 5
        else if (points(3) == points(6) .and. points(3) == blocks(nn)%kl) then
           blocks(nn)%b2bconn(mm)%faceID = 6
        end if

     end do ! 1ot1 Loop
  end do ! Master Zone Loop

  ! Now for the fun part the recursive propogation of each split to
  ! the rest of the blocks. This is by far the most tricky part of the
  ! entire operation

  ! Loop over each zone and BC, if we find a BC or n1to1 that needs to
  ! propogated we will do so

  do nn=1,nzones

     do mm=1,blocks(nn)%nbcs
        ! I face so j and k edges are split
        if (blocks(nn)%bcs(mm)%faceID == 1 .or. &
            blocks(nn)%bcs(mm)%faceID == 2) then
           call addSplit(nn,2,blocks(nn)%bcs(mm)%range(2))
           call addSplit(nn,2,blocks(nn)%bcs(mm)%range(5))
           call addSplit(nn,3,blocks(nn)%bcs(mm)%range(3))
           call addSplit(nn,3,blocks(nn)%bcs(mm)%range(6))
           ! J face so i and k edges are split
        else if (blocks(nn)%bcs(mm)%faceID == 3 .or. &
                 blocks(nn)%bcs(mm)%faceID == 4) then
           call addSplit(nn,1,blocks(nn)%bcs(mm)%range(1))
           call addSplit(nn,1,blocks(nn)%bcs(mm)%range(4))
           call addSplit(nn,3,blocks(nn)%bcs(mm)%range(3))
           call addSplit(nn,3,blocks(nn)%bcs(mm)%range(6))
           ! K face so i and j edges are split
        else if (blocks(nn)%bcs(mm)%faceID == 5 .or. &
                 blocks(nn)%bcs(mm)%faceID == 6) then
           call addSplit(nn,1,blocks(nn)%bcs(mm)%range(1))
           call addSplit(nn,1,blocks(nn)%bcs(mm)%range(4))
           call addSplit(nn,2,blocks(nn)%bcs(mm)%range(2))
           call addSplit(nn,2,blocks(nn)%bcs(mm)%range(5))
        end if
     end do
     do mm=1,blocks(nn)%n1to1
        ! I face so j and k edges are split
        if (blocks(nn)%b2bconn(mm)%faceID == 1 .or. &
            blocks(nn)%b2bconn(mm)%faceID == 2) then
           call addSplit(nn,2,blocks(nn)%b2bconn(mm)%range(2))
           call addSplit(nn,2,blocks(nn)%b2bconn(mm)%range(5))
           call addSplit(nn,3,blocks(nn)%b2bconn(mm)%range(3))
           call addSplit(nn,3,blocks(nn)%b2bconn(mm)%range(6))
           ! J face so i and k edges are split
        else if (blocks(nn)%b2bconn(mm)%faceID == 3 .or. &
                 blocks(nn)%b2bconn(mm)%faceID == 4) then
           call addSplit(nn,1,blocks(nn)%b2bconn(mm)%range(1))
           call addSplit(nn,1,blocks(nn)%b2bconn(mm)%range(4))
           call addSplit(nn,3,blocks(nn)%b2bconn(mm)%range(3))
           call addSplit(nn,3,blocks(nn)%b2bconn(mm)%range(6))
           ! K face so i and j edges are split
        else if (blocks(nn)%b2bconn(mm)%faceID == 5 .or. &
                 blocks(nn)%b2bconn(mm)%faceID == 6) then
           call addSplit(nn,1,blocks(nn)%b2bconn(mm)%range(1))
           call addSplit(nn,1,blocks(nn)%b2bconn(mm)%range(4))
           call addSplit(nn,2,blocks(nn)%b2bconn(mm)%range(2))
           call addSplit(nn,2,blocks(nn)%b2bconn(mm)%range(5))
        end if
     end do
  end do


  ! Finally loop over the number of user specified splits:
  do mm=1,n_user_split
     call addSplit(user_splits(mm,1),user_splits(mm,2),user_splits(mm,3))
  end do

  ! The Last step is to take all the blocks splits we have and to
  ! write a new CGNS File
  call cg_open_f(out_filename,CG_MODE_WRITE, cg_out, ier)
  if (ier .eq. CG_ERROR) call cg_error_exit_f

  print *,'Writing output file:',out_filename

  call cg_base_write_f(cg_out,"Base#1", Celldim,Physdim, base, ier)
  zoneCounter = 0
  j = 0
  do nn=1,nZones

     i = (blocks(nn)%nKSplit - 1)* &
          (blocks(nn)%nJSplit - 1)* &
          (blocks(nn)%nISplit - 1)
     j = j + i
     print *,'Block ',nn, 'is split into ',i,' pieces'
  end do
  print *,'There are ',j,' zones in the new file.'

  call cpu_time(time(2))     
  do nn=1,nZones
     do k=1,blocks(nn)%nKSplit - 1
        do j=1,blocks(nn)%nJSplit - 1
           do i=1,blocks(nn)%nISplit - 1             

              blockStart = (/blocks(nn)%iSplit(i),&
                             blocks(nn)%jSplit(j),&
                             blocks(nn)%kSplit(k)/)

              blockEnd   = (/blocks(nn)%iSplit(i+1),&
                             blocks(nn)%jSplit(j+1),&
                             blocks(nn)%kSplit(k+1)/)

              sizes(1:3) = abs(blockEnd - blockStart) + 1
              sizes(4:6) = abs(blockEnd - blockStart)
              sizes(7) = 0
              sizes(8) = 0
              sizes(9) = 0

999           FORMAT('zone',I4.4)
              write(zonename,999) zoneCounter
              zoneCounter = zoneCounter + 1           
              call cg_zone_write_f(cg_out,base,zonename,sizes,Structured,zone,ier)
              if (ier .eq. CG_ERROR)  call cg_error_exit_f
             
              allocate(tempx(sizes(1),sizes(2),sizes(3),3))

              ! Read Grid Section
              call cg_coord_read_f(cg_in,base,nn,'CoordinateX',RealDouble,&
                   blockStart,blockEnd,tempx(:,:,:,1),ier)
              call cg_coord_read_f(cg_in,base,nn,'CoordinateY',RealDouble,&
                   blockStart,blockEnd,tempx(:,:,:,2),ier)
              call cg_coord_read_f(cg_in,base,nn,'CoordinateZ',RealDouble,&
                   blockStart,blockEnd,tempx(:,:,:,3),ier)

              ! Write Grid Section
              call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
                   'CoordinateX',tempx(:,:,:,1), coordID,ier)
              call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
                   'CoordinateY',tempx(:,:,:,2), coordID,ier)
              call cg_coord_write_f(cg_out,base,zoneCounter,realDouble,&
                   'CoordinateZ',tempx(:,:,:,3), coordID,ier)

              deallocate(tempx)
           end do ! I loop 
        end do ! J loop
     end do !K loop
  end do ! Original Zone Loop

  ! Close both CGNS files
  call cg_close_f(cg_in, ier)
  call cg_close_f(cg_out,ier)

  ! Memory Cleanup
  deallocate(blocks)
  !deallocate(user_splits)
  call cpu_time(time(3))       

  ! Output some timing info
  print *,'Logic Time:',time(2)-time(1), ' seconds'
  print *,'File Write Time:',time(3)-time(2), ' seconds'

contains

  recursive subroutine addSplit(nn,idim,index)
    ! This is the master recusive subroutine. This subroutien attempts
    ! to add a split at index index, on block index, idim, (i, j or k). 

    implicit none
    integer :: nn,idim,index,mm,min,max,low,high
    logical :: added
    integer :: nn_new,idim_new, index_new
    integer :: this_block_size,offset
    integer :: donor_high,donor_low,abs_idim

    if (idim == 1) then
       added = addISplit(nn,index)
    else if (idim == 2) then
       added = addJSplit(nn,index)
    else if (idim == 3) then
       added = addKSplit(nn,index)
    end if

    if (.not. added) then
       return ! This split was already in so we just return
    else ! A new split was added...we need to keep going
       
       ! Loop over the b2b match conditions
          
       do mm=1,blocks(nn)%n1to1

          low =min(blocks(nn)%b2bconn(mm)%range(idim),&
               blocks(nn)%b2bconn(mm)%range(idim+3))
                
          high= max(blocks(nn)%b2bconn(mm)%range(idim),&
               blocks(nn)%b2bconn(mm)%range(idim+3))

          ! It must be fully contained in range low->high
          if (index > low .and. index < high)  then
             ! Now, we are ready to call addSplit again,
             ! recursively, on the connected block

             ! We must determine, nn,idim and index for the NEW block
             nn_new = findZoneName(blocks(nn)%b2bconn(mm)%donorName)
             
             idim_new = blocks(nn)%b2bconn(mm)%transform(idim)

             offset = index - low

             abs_idim = abs(idim_new)

             donor_high = max(blocks(nn)%b2bconn(mm)%donor_range(abs_idim), &
                              blocks(nn)%b2bconn(mm)%donor_range(abs_idim+3))

             donor_low = min(blocks(nn)%b2bconn(mm)%donor_range(abs_idim), &
                              blocks(nn)%b2bconn(mm)%donor_range(abs_idim+3))

             if (idim_new > 0) then
                index_new = donor_low + offset
             else 
                index_new = donor_high - offset
             end if

             call addSplit(nn_new,abs_idim,index_new)

          end if
       end do
    end if
    
  end subroutine addSplit

  subroutine addReqSplit(nn,points)
    implicit none
    integer :: nn, points(6)
    logical :: res
    if ( (abs(points(4)-points(1)) > 0 .and. abs(points(5)-points(2)) >0 ) .or. &
         (abs(points(4)-points(1)) > 0 .and. abs(points(6)-points(3)) >0 ) .or. &
         (abs(points(5)-points(2)) > 0 .and. abs(points(6)-points(3)) >0 )) then

       res = addISplit(nn,points(1))
       res = addISplit(nn,points(4))
       res = addJSplit(nn,points(2))
       res = addJSplit(nn,points(5))
       res = addKSplit(nn,points(3))
       res = addKSplit(nn,points(6))

    end if
  end subroutine addReqSplit
  
  function addISplit(nn,index)
    ! Add an split to I direction at Index to block nn
    implicit none
    integer :: nn,index
    logical :: alreadyInList,addISplit
    addISplit = .False.
    alreadyInList = inList(blocks(nn)%iSplit(1:blocks(nn)%nISplit), &
         blocks(nn)%nIsplit,index)
    if (.not. alreadyInList) then
       blocks(nn)%nISplit = blocks(nn)%nIsplit + 1
       blocks(nn)%iSplit(blocks(nn)%nISplit) = index
       call bubble_sort(blocks(nn)%iSplit(1:blocks(nn)%nISplit),blocks(nn)%nIsplit)
       addISplit = .True. 
    end if
  end function addISplit

  function addJSplit(nn,index)
    ! Add an split to J direction at Index to block nn
    implicit none
    integer :: nn,index
    logical :: alreadyInList,addJSplit
    addJSplit = .False.
    alreadyInList = inList(blocks(nn)%jSplit(1:blocks(nn)%nJSplit), &
         blocks(nn)%nJsplit,index)
    if (.not. alreadyInList) then
       blocks(nn)%nJSplit = blocks(nn)%nJsplit + 1
       blocks(nn)%jSplit(blocks(nn)%nJSplit) = index
       call bubble_sort(blocks(nn)%jSplit(1:blocks(nn)%nJSplit),blocks(nn)%nJsplit)
       addJSplit = .True. 
    end if
  end function addJSplit

  function addKSplit(nn,index)
    ! Add an split to K direction at Index to block nn
    implicit none
    integer :: nn,index
    logical :: alreadyInList,addKSplit
    addKSplit = .False.
    alreadyInList = inList(blocks(nn)%kSplit(1:blocks(nn)%nKSplit), &
         blocks(nn)%nKsplit,index)
    if (.not. alreadyInList) then
       blocks(nn)%nKSplit = blocks(nn)%nKsplit + 1
       blocks(nn)%kSplit(blocks(nn)%nKSplit) = index
       call bubble_sort(blocks(nn)%kSplit(1:blocks(nn)%nKSplit),blocks(nn)%nKsplit)
       addKSplit = .True. 
    end if
  end function addKSplit

  function inList(list,n,checkVal)
    ! Check if checkVal is in list of length n. Once again we use an
    ! inefficient linear search
    implicit none
    integer :: n,list(n),checkVal
    logical :: inList
    integer :: i
    inList = .False.
    do i=1,n
       if (checkVal == list(i)) then
          inList = .True. 
       end if
    end do
  end function inList

  function findZoneName(zoneName)
    ! Find the zone number cooresponding to zoneName'
    ! We use an (inefficient) linear search
    implicit none
    character*32 :: zoneName
    integer :: findZoneName, n

    zoneloop: do n=1,nZones
       if (zoneName == blocks(n)%zoneName) then
          findZoneName = n
          exit zoneloop
       end if
    end do zoneloop
  end function findZoneName

  subroutine bubble_sort (x,n)
    ! Sort x in ASCENDING order

    IMPLICIT NONE

    integer :: n, x(n)
    integer :: temp, i,j
    
    do i=1,n-1
       do j=1,n-1
          if (X(J+1) .lt. X(J)) then 
             temp = x(j)
             x(j) = x(j+1)
             x(j+1) = temp
          end if
       end do
    end do
  end subroutine bubble_sort

end program cgns_split
