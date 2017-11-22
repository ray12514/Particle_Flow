c------------------------------------------------------------------------c
c                   Post process data                                    c
c------------------------------------------------------------------------c 
      subroutine read_stats
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'AVG'
      include 'PARTICLES'

      character*80 filename(9999)
      character*3 prefix (3)
      character*128   command(2)
      character*128  char_filename
      
      real work1(lx1,ly1,lz1,lelt)
      real work2(lx1,ly1,lz1,lelt,ldim)
      integer icalld
      save    icalld
      data    icalld  /0/
      
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelv
      
      if (param(11).nid.eq.0) then 
        write(6,*)'Set NSTEPS to zero !!' 
        call exitt
      endif  
      
      if (icalld.eq.0) then
         icalld = icalld + 1
         call rzero(STAT,ntot*n_var) !may change this but for now this works
         call rzero(CSTAT,ntot*80) 
      endif

!      ifreguo = .true.   ! dump on regular (uniform) grid instead of GLL
!      nrg     = 16   ! dimension of regular grid (nrg**ndim)
          
  
      !Will do it this way first  
         call system("ls -1 f*elbow0.f* | wc -l > file.list")
         call system("ls f*elbow0.f* >> file.list")
        

      if (nid.eq.0) then
        open(unit=199,file='file.list',
     $    form='formatted',status='old')
          read(199,*) nfiles
          read(199,'(A80)') (filename(i),i=1,nfiles)
          close(199)
        nperfile=nfiles/n_st_f
      endif
      ! read file-list
      call bcast(nfiles,isize)
      call bcast(nperfile,isize)
      call bcast(filename,nfiles*80)
      
      
      do i = 1,nfiles
        do j = 1,nperfile
          call load_fld(filename(i))
        if(nid.eq.0) write(6,*) 'Loaded file ', filename(i)
          
          !Write data to big array
          skip=4*(i-1))
          call opadd2 (STAT(1,1+skip),STAT(1,2+skip),
     $       STAT(1,3+skip),vx,vy,vz)
          
          call add2(STAT(1,4+skip),t(1,1,1,1,1),ntot)    
        
        enddo
      enddo
        
        scale=1./nperfile
        !Adjust to the stats according to weight  
        call cmult(STAT,scale,ntot*nvar)
        
      
       return
       end
       
c----------------------------------------------------------------------------------       
      subroutine MEAN_RSTENSOR_
      include 'SIZE'
      include 'TOTAL'
      include  'PARTICLES'
      real work1(lx1*ly1*lz1*lelt)
      real work2(lx1*ly1*lz1*lelt,ldim)
      real work3(lx1*ly1*lz1*lelt,ldim)

      
      
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelv
      
      !Move the mean to the outpost array CSTAT
      call opcopy(CSTAT(1,1),CSTAT(1,2),CSTAT(1,3)
     $            STAT(1,1),STAT(1,2),STAT(1,3))
      

      call opcopy(work2(1,1),work2(1,2),work2(1,3)
     $            STAT(1,1),STAT(1,2),STAT(1,3))
  
      call col3(work2(1,1),work2(1,1),work2(1,1),ntot)
      call col3(work2(1,2),work2(1,2),work2(1,2),ntot)
      call col3(work2(1,3),work2(1,3),work2(1,3),ntot)

      !Reynolds Stress tensor
!      __
!      U^2 
      call sub3(CSTAT(1,4),STAT(1,5),work2(1,1),ntot)
!      __
!      V^2 
      call sub3(CSTAT(1,5),STAT(1,6),work2(1,2),ntot)
!      __
!      W^2 
      call sub3(CSTAT(1,6),STAT(1,7),work2(1,3),ntot)

      call opcopy(work2(1,1),work2(1,2),work2(1,3)
     $            STAT(1,1),STAT(1,2),STAT(1,3))
      
      !Copy means back to work array
      call opcopy(work2(1,1),work2(1,2),work2(1,3)
     $            STAT(1,1),STAT(1,2),STAT(1,3))

      !<U><V> 
      call col3(work3(1,1),work2(1,1),work2(1,2),ntot)
      !<U><W> 
      call col3(work3(1,2),work2(1,1),work2(1,3),ntot)
      !<V><W>
      call col3(work3(1,3),work2(1,2),work2(1,3),ntot)
!      __
!      UV 
      call sub3(CSTAT(1,7),STAT(1,9),work3(1,1),ntot)
!      __
!      UW 
      call sub3(CSTAT(1,8),STAT(1,10),work3(1,2),ntot)
!      __
!      VW 
      call sub3(CSTAT(1,9),STAT(1,11),work3(1,3),ntot)
      
                  !d<U>/dx     d<U>/dy     d<U>/dz    
      call gradm1(CSTAT(1,12),CSTAT(1,13),CSTAT(1,14),work2(1,1))
                  !d<V>/dx     d<V>/dy     d<V>/dz  
      call gradm1(CSTAT(1,15),CSTAT(1,16),CSTAT(1,17),work2(1,2))
                  !d<W>/dx     d<W>/dy     d<W>/dz  
      call gradm1(CSTAT(1,18),CSTAT(1,19),CSTAT(1,20),work2(1,3))

      
      
      

         
      
      ifxyo=.true.
!       call outpost (u_avg,v_avg,w_avg,pr_avg,t,   'av_')
!       call outpost (u_rms,v_rms,w_rms,pr_rms,t,   'uu_')
      ifpo=.false. 
!       call outpost (uv_rms,uw_rms,vw_rms,pr,t,    'uv_')
      ifxyo=.false.
      ifpo=.true. 

      call exitt      
      return
      end
      
c---------------------------------------------------------------------------------------

