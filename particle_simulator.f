c---------------------------------------------------------------------c
c             Particle Tracking Routine for Nek5000                   
c   Add "call particle_generator" to your .usr file        
c   in the userchk. To initialize particle location update            
c     particle_init subroutine                                        
c---------------------------------------------------------------------c
      subroutine set_part_pointers
      implicit none  
      include 'SIZE'
      include 'PARTICLES' 
c     Minimal value of ni = 5
c     Minimal value of nr = 14*ndim + 1
      integer n,nr,ni
      common  /iparti/ n,nr,ni
!       common /ptpointers/ jrc,jpt,je0,jps,jai,jgn,nai,jr,jd,jx,jy,jz,jx1
!      $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
!      $         ,ju0,ju1,ju2,ju3,jgu,jgv,jgw,jwd,jwn,jpd,jrh,jdt,jar,nar
      jrc = 1 ! Pointer to fgslib_findpts return code
      jpt = 2 ! Pointer to fgslib_findpts return processor id
      je0 = 3 ! Pointer to fgslib_findpts return element id
      jps = 4 ! Pointer to proc id for data swap
      jai = 5 ! Pointer to auxiliary integers
      jgn = 6 ! Pointer to group number
      jcl = 7 ! Pointer to if the particle should be recycled 
      jcp = 8 ! Pointer to number of times ID has been copied

      nai = ni - (jcp-1)  ! Number of auxiliary integers
      if (nai.le.0) call exitti('Error in nai:$',ni)

      jr  = 1         ! Pointer to fgslib_findpts return rst variables
      jd  = jr + ndim ! Pointer to fgslib_findpts return distance
      jx  = jd + 1    ! Pointer to fgslib_findpts input x value
      jy  = jx + 1    ! Pointer to fgslib_findpts input y value
      jz  = jy + 1    ! Pointer to fgslib_findpts input z value

      jx1 = jx + ndim ! Pointer to xyz at t^{n-1}
      jx2 = jx1+ ndim ! Pointer to xyz at t^{n-2}
      jx3 = jx2+ ndim ! Pointer to xyz at t^{n-3}
      
      ja0 = jx3+ ndim ! Pointer to current particle acceleration
      ja1 = ja0+ ndim ! Pointer to particle acceleration at t^{n-1}
      ja2 = ja1+ ndim ! Pointer to particle acceleration at t^{n-2}
      ja3 = ja2+ ndim ! Pointer to particle acceleration at t^{n-3}
      
      jv0 = ja3+ ndim ! Pointer to current particle velocity
      jv1 = jv0+ ndim ! Pointer to particle velocity at t^{n-1}
      jv2 = jv1+ ndim ! Pointer to particle velocity at t^{n-2}
      jv3 = jv2+ ndim ! Pointer to particle velocity at t^{n-3}

      ju0 = jv3+ ndim ! Pointer to current fluid velocity
      ju1 = ju0+ ndim ! Pointer to fluid velocity at t^{n-1}
      ju2 = ju1+ ndim ! Pointer to fluid velocity at t^{n-2}
      ju3 = ju2+ ndim ! Pointer to fluid velocity at t^{n-3}


      jgu = ju3+ ndim ! Pointer to u grad at current timestep
      jgv = jgu+ ndim ! Pointer to v grad at current timestep
      jgw = jgv+ ndim ! Pointer to w grad at current timestep
      
      jwd = jgw+ndim  ! Pointer to distance to nearest wall
      jwn = jwd+ndim  ! Pointer to nearest wall coordinates
      jpd = jwn+ndim  ! Pointer to particle diameter
      jrh = jpd+1     ! Pointer to particle density 
      jdt = jrh+1     ! Pointer to du/dt at particle location   
      jar = jdt+ndim  ! Particle to stokes number
      jds = jar+2*ndim ! Pointer to stress tensor
      jss = jds+1     ! Pointer to SijSji
      nar = nr - (jss-1)  ! Number of auxiliary reals

      if (nar.le.0) call exitti('Error in nar:$',nr)
      return
      end
c----------------------------------------------------------------------
      subroutine particle_generator ! Particle injection
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      
!       include 'CTIMER'
      integer lr,li
      parameter (lr=27*ldim+1,li=8+1)
      real rpart
      common /rparts/ rpart(lr,lhis)
      integer ipart
      common /iparts/ ipart(li,lhis)
      integer n,nr,ni
      common /iparti/ n,nr,ni
      real  xdrange(2,3) 
      common /domainrange/ xdrange

      integer icalld,ipstep
      save    icalld
      data    icalld  /0/
      
      if(icalld.eq.0) then
        nr=lr
        ni=li
        icalld=1
        
        call rzero(rpart,lr*lhis)
        call izero(ipart,li*lhis)
        call domain_size(xdrange(1,1),xdrange(2,1),xdrange(1,2)
     $                  ,xdrange(2,2),xdrange(1,3),xdrange(2,3))
        call particle_start (rpart,nr,ipart,ni,n)
      endif 
        
        call particle_advect_up  (rpart,nr,ipart,ni,n)
      
      
      if(istep.ge.1000) then
        if(mod(istep,iostep/1000).eq.0)then     
          call volume_stats(rpart,nr,ipart,ni,n)
        endif
      endif
      
      ipstep = iostep
      call output_particles   (rpart,nr,ipart,ni,n,ipstep)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine particle_start(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
      
      call set_part_pointers
      
      if(intpart) then    
        call particle_init       (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      else
        call particle_restart    (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)

      endif
      call particle_interp_all   (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)

      return
      end
c-----------------------------------------------------------------------
      subroutine particle_init (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
c     This version distributes particles throughout the entire domain
c     for specified volume fraction 
c
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      include 'RECYCLE'

      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
      
      real partdiam(pgrp),partdens(pgrp)
      real  xdrange(2,3) 
      common /domainrange/ xdrange
      save partdiam,partdens
      data partdiam /1e-3/
      data partdens /2/

      integer lcount,icount
      save    lcount,icount
      data    lcount,icount /0,0/
      integer llhis,i,j,k,l,nw,npart
      real xstart,ystart,zstart
      real xlen,ylen,zlen
      real ran2,dum,dumx,dumy,dumz


      llhis = lhis
      dum=ran2(-1) !start random number generator

      k  = icount       ! icount = total # particles emitted
      l  = lcount       ! Index into local particle array
      !Will add a routine to track the particles that are upstream and downstream
      !one hundred total particle tracks
      nw = 1 !Number of particles per group  

      do i=1,pgrp
        
        do j=1,nw
          dumx = ran2(2)
          dumy = ran2(2)
          dumz = ran2(2)
          if (j.lt.nw/2) then
            xstart = -upstream_pos
            ystart = -0.49
            zstart = -0.49
            xlen   = upstream_pos+0.48
            ylen   = 0.95
            zlen   = 0.95
          else
            xstart = -0.49
            ystart = -0.50
            zstart = -0.49
            ylen   = -(upstream_pos-0.6)
            xlen   = 0.95
            zlen   = 0.95
          endif

           k = k+1       ! Global particle count
              
           if (mod(k,np).eq.nid) then ! Inject particle _only_ on this processor
              l=l+1  ! local count
              if (l.gt.llhis)then 
               write(6,*)'Not enough space to store more particles'
               call exitt
              endif
              real_parts(jx,l) = xstart + dumx*(xlen -partdiam(i))  
              real_parts(jy,l) = ystart + dumy*(ylen -partdiam(i))   ! Particle Initial position 
              real_parts(jz,l) = zstart + dumz*(zlen -partdiam(i))   !
              real_parts(jx,l) =-6.0
              real_parts(jy,l) =-0.5+partdiam(i)
!                if(j.eq.1)   real_parts(jx,l) = -5.5
!                if(j.eq.2)   real_parts(jy,l) = -1.5
              real_parts(jz,l) = 0.0
              real_parts(jrh,l)=partdens(i)   !Particle density
              real_parts(jpd,l)=partdiam(i)!Particle diameter
              real_parts(jar,l)= 0.0  
              integer_parts(jai,l)=k   !Particle ID
              integer_parts(jgn,l)=i   !Particle Group ID 
              if (abs(real_parts(jx,l)).ge.recycle_pos) then 
                integer_parts(jcl,l)=1  !Particle gets recycled
                !Particles that will be tracked in the domain
                if (k.le.100) integer_parts(jcl,l)=-1 
              else
                integer_parts(jcl,l)=2  !Particle not in recycled portion 
                !Particles that will be tracked in the domain
                if (k.le.100) integer_parts(jcl,l)=-2    
              endif
                integer_parts(jcp,l)=1  !Counter for particle copies in he domain
           endif
         
        enddo
      enddo  
      
      if(nid.eq.0) write(6,*) k ,'total particle count'
        icount = k
        lcount = l
        npart=0
        npart  = max(npart,lcount)
        num_total= npart
        write(6,'(I8,A,I8)') num_total, ' particles on proc: ', nid
        call particle_interp_all (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)

        do j=1,num_total
          do i=0,ndim-1                   !Update the initial particle
            real_parts(jv0+i,j)=real_parts(ju0+i,j) !velocity to fluid velocity                              
          enddo
        enddo
       call Force (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
       call particle_collect (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)

       

      return
      end
c---------------------------------------------------------------------------
      !update this routine to include the case where the number processor
      !change for the next run
      subroutine particle_restart(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      integer i,k
      character*128 filename
      integer icalld
      save    icalld
      data    icalld /0/
      
      write(filename,1) nid+1, 2
 1      format('PART/',I4.4,'_rstpart',i5.5,'.3D')
      open(unit=98,file=filename)
!  900  write(6,*) 'no data on this processor ', nid+1       
!             if (nid.eq.0) then
!            write(6,*) 'Loaded file ', filename(icalld)
!             endif
      read(98,'(I8)') num_total
      read(98,4)((real_parts(k,i),k=1,num_reals)
     $     ,(integer_parts(k,i),k=1,num_ints),i=1,num_total)
 4      format(1p75e17.9,9I9)
      close(98)
          
!       if(nid.eq.0) write(6,*) k ,'total particle count'
      write(6,'(I4,A,I4)') num_total, ' particles on proc: ', nid
  
      return
      end
! c---------------------------------------------------------------------------
      subroutine particle_advect_up(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none 
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'PARTICLES' 
!       real    pt_timers(12), scrt_timers(12)
!       common /trackingtime/ pt_timers, scrt_timers
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
!       common /myparts/ times(0:3),alpha(0:3),beta(0:3)
!       real x_ex,u_ex,v_ex,Stokes_ex,time_ex,g_ex
!       common /exact_sol/ x_ex(3,lhis),u_ex(3,lhis),
!      $                   v_ex(3,lhis),Stokes_ex(lhis),
!      $                    time_ex,g_ex(3),err_ex(3,lhis)

!        integer excalld
!        save excalld
!        data excalld /0/ 
!       real g(3),tau,c(3),rha
       integer i,j,jx0
       real c(3),x_tmp,v_tmp,x_o(3),v_o(3) ,u_o(3)  
       save x_o,v_o,u_o
       
!       data g /0.0,9.80665,0.0/
      
!       tau=100
      call particle_AB3_coef(c)
!       call get_bdf_ext_coefs(beta,alpha,times)
      jx0=jx
       !!!!!!!!REDO THIS ROUTINE!!!!!!!!!!!!!! 
c      write(6,*) rpart(jv0+1,1), 'velocity'
c      write(6,*) rpart(jx0+1,1), 'y-position'
c           Solve for velocity at time t^n

      !!!Don't really need the old fluid velocity data 
      !!!or the particle data, and I only need to
      !!!store up t^-2 for AB3
      call Force (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
!!!!Used for the Exact solution !!!!!!!!!!!!!
!       if (excalld.eq.0) then 
!         call rzero(x_ex,lhis)
!         call rzero(u_ex,lhis)
!         call rzero(v_ex,lhis)
!         call rzero(x_o,ndim)
!         call rzero(u_o,ndim)
!         call rzero(v_o,ndim)
! !         call rzero(Stokes_ex,lhis)
!         time_ex=0
! !         jx0=jx
!         do i=1,n
!           do ii=1,ndim 
!             x_o(ii)=rpart(jx0+ii-1,i)
!             v_o(ii)=rpart(jv0+ii-1,i) 
!             u_o(ii)=rpart(ju0+ii-1,i) 
!           enddo 
!         enddo 
          
!          excalld=1
!       endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!
      !!!!!   I need to rewrite all this code 
      !!!!!   so that it is just one big loop 
      !!!!!   that performs all the operations
      !!!!!   before I make the call to the fgslib_findpts
      !!!!!   and fgslib_crystal router routine
      !!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=0,ndim-1
      do i=1,num_total
         real_parts(ja3+j,i)=real_parts(ja2+j,i)
         real_parts(ja2+j,i)=real_parts(ja1+j,i)
         real_parts(ja1+j,i)=real_parts(ja0+j,i)
         real_parts(ja0+j,i)=0.0
         real_parts(ju3+j,i)=real_parts(ju2+j,i)
         real_parts(ju2+j,i)=real_parts(ju1+j,i)
         real_parts(ju1+j,i)=real_parts(ju0+j,i)
         real_parts(jv3+j,i)=real_parts(jv2+j,i)
         real_parts(jv2+j,i)=real_parts(jv1+j,i)
         real_parts(jv1+j,i)=real_parts(jv0+j,i)
         real_parts(jx3+j,i)=real_parts(jx2+j,i)
         real_parts(jx2+j,i)=real_parts(jx1+j,i)
         real_parts(jx1+j,i)=real_parts(jx0+j,i)
      enddo
      enddo

      
!      
      do i=1,num_total
        do j=0,ldim-1
            
!   c       real_parts(ja0+j,i)=(real_parts(ju1+j,i)-real_parts(jv1+j,i))/tau-g(j+1)
!   c       if(i.eq.1) 
!   c       write(6,*)real_parts(ja0+j,i),j+1,'acceleration',nid,'proc',i,'pt'
          
          !rha=(-real_parts(jv1+j,i))/tau-g(j+1)
!   c        if(istep.eq.0) rha=0.0
          real_parts(jv0+j,i) = real_parts(jv1+j,i)+
     $               dt*(c(1)*real_parts(ja1+j,i)
     $               +c(2)*real_parts(ja2+j,i) 
     $               + c(3)*real_parts(ja3+j,i))
          real_parts(jx0+j,i) =real_parts(jx1+j,i)+
     $              dt*(c(1)*real_parts(jv1+j,i)
     $          + c(2)*real_parts(jv2+j,i)
     $        + c(3)*real_parts(jv3+j,i))
        enddo 
      enddo
!        write(6,*) alpha(1),alpha(2),alpha(3)    
!        write(6,*) beta(1),beta(2),beta(3)    
!        do i=1,n
!         do j=0,ndim-1
!           rhs =     alpha(1)*real_parts(ja1+j,i)
!      $            + alpha(2)*real_parts(ja2+j,i)
!      $            + alpha(3)*real_parts(ja3+j,i)
!      $        +     beta (1)*real_parts(jv1+j,i)
!      $        +     beta (2)*real_parts(jv2+j,i)
!      $        +     beta (3)*real_parts(jv3+j,i)
!           real_parts(jv0+j,i) = rhs / beta(0) ! Implicit solve for v
!           rhx = beta (1)*real_parts(jx1+j,i)
!      $        + beta (2)*real_parts(jx2+j,i)
!      $        + beta (3)*real_parts(jx3+j,i) + real_parts(jv0+j,i)
!           real_parts(jx0+j,i) = rhx / beta(0)     ! Implicit solve for x
!         enddo
!       enddo
      call particle_collect (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      call particle_interp_all (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
    
!      call update_particle_location (real_parts,nr,n)

!!!!!!!!!!!!!!! This was to test the Advection shceme !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Only work when Re_p < 1 !!!!!!!!!!!!!!!!!!!!!!
! !       time_ex=time_ex+dt
!       do i=1,n
!         do j=0,ndim-1
!         !!!Exact solution
!       v_tmp=v_o(j+1)+(exp(time_ex/Stokes_ex(i))-1)
!      $        *(Stokes_ex(i)*g_ex(j+1)+u_o(j+1))
!         v_ex(j+1,i)=exp(-time_ex/Stokes_ex(i))*(v_tmp)
! !         v_ex(j+1,i)=v_ex(j+1,i)*(exp(time_ex/Stokes_ex(i)))
! !      $     -(exp(time_ex/Stokes_ex(i))-1)
! !      $       *(-Stokes_ex(i)*g_ex(j+1)+u_ex(j+1,i))
        
!         x_tmp=exp(-(time_ex)/Stokes_ex(i))
!         x_ex(j+1,i)= x_o(j+1)+v_o(j+1)*Stokes_ex(i)*(x_tmp 
!      $   + 1)-(u_o(j+1)+ Stokes_ex(i)*g_ex(j+1))*(time_ex
!      $    -Stokes_ex(i)*(1-x_tmp))
        
!         err_ex(j+1,i)=abs(v_ex(j+1,i)-real_parts(jv1+j,i))
       
!         enddo
!       enddo 
      
!       write(6,*) v_ex(1,1),err_ex(1,1),time_ex,'check solution_x'
!       write(6,*) v_ex(2,1),err_ex(2,1),time_ex,'check solution_y'
!       write(6,*) v_ex(3,1),err_ex(3,1),time_ex,'check solution_z'
! !       write(6,*) c(1),c(2),c(3),'ccc' 
! !       call exitt
!       time_ex=time_ex+dt
      
      return 
      end       
c-----------------------------------------------------------------------
      subroutine output_particles (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total,particle_io)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      
      integer num_reals,num_ints,num_total,particle_io
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      integer icallp,rstpart
      save    icallp,rstpart
      data    icallp,rstpart /0,0/
      integer i,k
      character*128 fname,FRMT
      
        write(FRMT,4) num_reals,num_ints
 4      format('(1p',I2,'e17.9,',I1,'I9)')
      if (icallp.eq.0.or.mod(istep,particle_io).eq.0
     $   .or.lastep.eq.1) then 
        icallp = icallp+1 
        write(fname,1) nid+1,icallp
 1      format('PART/',I4.4,'_part',i5.5,'.3D')
        open(unit=72,file=fname)!
        write(72,'(I8)') num_total
        write(72,2)((real_parts(k,i),k=jx,jz) ! xyz particles    
     $     ,(real_parts(k,i),k=jv0,jv0+2)     ! velocity particles
!      $     ,(real_parts(k,i),k=ju0,ju0+2)
!      $     ,(real_parts(k,i),k=ja0,ja0+2)
     $     ,(integer_parts(k,i),k=jai,jgn),i=1,num_total) ! ID and Group number
  2      format(1p6e17.9,2I9)
        close(72)
       if(nid.eq.0) write(6,*)'particle output', icallp  
      endif
      !Restart Data
      if (mod(istep,iostep).eq.2.or.lastep.eq.1) then
        rstpart=rstpart+1
        write(fname,3) nid+1,rstpart        
 3      format('PART/',I4.4,'_rstpart',i5.5,'.3D')
        open(unit=73,file=fname)
        write(73,'(I8)') num_total  
        write(73,FMT=FRMT)((real_parts(k,i),k=1,num_reals)
     $     ,(integer_parts(k,i),k=1,num_ints),i=1,num_total)
!  4      format(1p75e17.9,9I9)
        close(73)
        if(nid.eq.0) write(6,*)'restartp output', rstpart 
        if(rstpart.eq.2) rstpart=0
      endif
      
!       call exitt
      return
      end

c-----------------------------------------------------------------------
      subroutine particle_interp_all(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      integer mid,mp,nekcomm,nekgroup,nekreal 
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
      real  wrk(lx1*ly1*lz1*lelt,ldim)  
      common /outtmp/ wrk
      integer nmax
      parameter(nmax=lhis)
      integer lrf,lif 
      parameter (lrf=27*ldim+1,lif=8+1)
      real               loc_rpts(lrf,lhis),rem_rpts(lrf,lhis)
      common /fptspartr/ loc_rpts,rem_rpts
      integer            loc_ipts(lif,lhis),rem_ipts(lif,lhis)
      integer            loc_ptsmap(lhis),rem_ptsmap(lhis)
      integer            n_loc_pts,n_rem_pts 
      common /fptsparti/ loc_ipts,rem_ipts,n_loc_pts,n_rem_pts 
      common /mapss/     loc_ptsmap,rem_ptsmap
      integer ntot,nxyz,num_log
      real tolin
      integer icalld1,wrk_n,i
      save    icalld1
      data    icalld1 /0/

            
      logical partl         ! This is a dummy placeholder, used in cr()
      num_log = 0                ! No logicals exchanged
c      write(6,*) size(loc_ptsmap),'size of local maps' 
      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt
      call opcopy(wrk(1,1),wrk(1,2),wrk(1,3),vx,vy,vz)
      
      if (icalld1.eq.0) then    ! interpolation setup !? intpts_done(ih_intp_v)?
        icalld1 = icalld1+1
        tolin  = 1.e-12
        if (wdsize.eq.4) tolin = 1.e-6
        call fgslib_crystal_setup (cr_handle,nekcomm,np)
        call intpts_setup_part(tolin,ihandle_local)
c        call intpts_setup(tolin,ihandle_remote)
      endif

        if(nid.eq.0) write(6,*) 'call fgslib_findpts'
         call fgslib_findpts(ihandle_local,integer_parts(jrc,1),lif,
     &               integer_parts(jpt,1),lif,
     &               integer_parts(je0,1),lif,
     &               real_parts(jr, 1),lrf,
     &               real_parts(jd, 1),lrf,
     &               real_parts(jx, 1),lrf,
     &               real_parts(jy, 1),lrf,
     &               real_parts(jz, 1),lrf, num_total)
       
        jps = jai-1     ! Pointer to temporary proc id for swapping
         do i=1,num_total        ! Can't use jpt because it messes up particle info
            integer_parts(jps,i) = integer_parts(jpt,i)
         enddo
!       
         call fgslib_crystal_tuple_transfer(cr_handle,num_total,lhis
     $              , integer_parts,num_ints,partl
     $              ,num_log,real_parts,num_reals,jps)
          

c        Sort by element number - for improved local-eval performance
         
!          call fgslib_crystal_tuple_sort    (cr_handle,n 
!      $              , integer_parts,ni,partl,nl,real_parts,nr,je0,1)


c             Interpolate the particles locally          
         call particle_interp_local (real_parts,num_reals,integer_parts
     $       ,num_ints,num_total,wrk,ju0,ndim)
     
      return
      end
c------------------------------------------------------------------------- 
      subroutine intpts_setup_part(tolin,ih)
      implicit none     
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c     
      include 'SIZE'
      include 'GEOM'
      real tolin,tol
      integer ih
      integer nidd,npp,nekcomm,nekgroup,nekreal 
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
      integer n,npt_max,nxf,nyf,nzf,bb_t
      tol = tolin
      if (tolin.lt.0) tol = 1e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by
c
      if(nidd.eq.0) write(6,*) 'initializing intpts(), tol=', tol
      call fgslib_findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
      return
      end        
c------------------------------------------------------------------------

      subroutine particle_interp_local  (loc_r_data,nrf,loc_i_data
     &           ,nif,num_pts,wrk,pt,noflds)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      real loc_r_data(nrf,*)
      integer loc_i_data(nif,*),pt,nrf,nif,num_pts,noflds
      real wrk(1)
      integer ltot,nflds,ifld,iin,iout,is_out
      
      
!       common /iv_intp/ ihandle_local,ihandle_remote,cr_handle      
      ltot = lelt*lx1*ly1*lz1
      nflds  = noflds
c      write(6,*) pt, 'value of pointer'
c      write(6,*)nrf,'nrf',nif,'nif'
c      write(6,*) num_pts,'is this the right number of particles?'
        do ifld = 1,nflds
          iin    = (ifld-1)*ltot + 1
          iout   = (ifld-1)*num_pts + 1
          is_out = 1
            iout   = ifld
            is_out = nflds
          call fgslib_findpts_eval_local(ihandle_local,
     &                     loc_r_data(pt+iout-1,1),nrf,
     &                     loc_i_data(je0,1),nif,
     &                     loc_r_data(jr ,1),nrf,num_pts,
     &                     wrk(iin))
        enddo
        
c      do i=1,num_pts
c      write(6,'(A,1p3e15.7,I3)')'U',(loc_r_data(ju0+k,i),k=0,ndim-1)
c     & ,loc_i_data(je0,i)
c      write(6,'(A,1p3e15.7,I3)')'xyz',(loc_r_data(jx+k,i),k=0,ndim-1)
c     & ,loc_i_data(jai,i)
c      enddo
      
       if (nid.eq.0)
     $ write(6,*) 'local particle interpolation complete'
      
      return
      end
c---------------------------------------------------------------------
      subroutine particle_collect(real_parts,num_reals,
     $   integer_parts,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      include 'RECYCLE'
      include 'COLLISIONS'  
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer e,f,eg,ifld !nek5000 syntax
      integer i,nl !index for particles
      integer icalld
      save    icalld
      data    icalld /0/
      
      integer ntotal
      !!!questions about runtal and cl_up
      logical log,col_temp_l
      logical partl         ! This is a dummy placeholder, used in cr()
      nl = 0                ! No logicals exchanged
      
     
     
      
      
      if (icalld.eq.0) then

        ifld=1
        
!         call rzero(wall_area,nelt)
        call rzero(col_temp_r,rc_index*(2*lhis))  
!         call rzero(col_data_r,rc_index*(2*lhis))  
        call izero(col_temp_i,ic_index*(2*lhis))  
!         call izero(col_data_i,ic_index*(2*lhis))  
        ntotal=nx1*ny1*nz1*nelv
        call cheap_dist(m1,ifld,bound)
        !call distf(d,ifld,b,dmin,emin,xn,yn,zn) !wall
        !call opcopy(wrk(1,1),wrk(1,2),wrk(1,3),xn,yn,zn)
        call gradm1(wrk_n(1,1),wrk_n(1,2),wrk_n(1,3),m1)
        
        !call copy(wrk2,d,ntotal) !wall distance
        call copy(wrk2,m1,ntotal) !wall distance using psuedo

        call fgslib_crystal_setup (col_handle,nekcomm,np) 
        call fgslib_crystal_setup (exit_handle,nekcomm,np)       

        icalld=1
      endif

      !Interpolate the gradient at each partitcle location
      call particle_interp_local (real_parts,num_reals,integer_parts,
     $   num_ints,num_total,wrk_n,jwn,ndim)

      !Interpolate the distance to the wall for each particle
      call particle_interp_local (real_parts,num_reals,integer_parts,
     $   num_ints,num_total,wrk2,jwd,1)

      !Start the collision collection process  
      if(nid.eq.0) write(6,*) 'Checking for collisions...'
      i=0 
      do while (i.lt.num_total.and.num_total.gt.0)
        i=1+i  
        call collisions(real_parts,num_reals,
     $   integer_parts,num_ints,num_total,i) 
        call particle_position (real_parts,num_reals,
     $   integer_parts,num_ints,num_total,i)      
        
      enddo
      call nekgsync
      !call collision_detection
      return
      end
c-----------------------------------------------------------------------
      subroutine collisions(real_parts,num_reals,
     $   integer_parts,num_ints,num_total,p_index)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      include 'RECYCLE'
      include 'COLLISIONS' 
      integer num_reals,num_ints,num_total,p_index
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
      !All variables that belong in this routine
      real opp_dis,rslt_lngth,col_theta,wallpt  
      integer e,f,eg
      integer cl_c_up,colin,mymax2
      real normc,normc_check,dist3d,diff
      real prad
      integer e_col,eg_col,jx0
      integer cl_c
      save    cl_c
      data    cl_c /0/
      character*12 col_position
      col_position='collisns'
      !Start the collision collection process  
      
      jx0 = jx 
     

        
        prad=real_parts(jpd,p_index)*0.5 !Particle radius 
        !Check the wall distance to update for collisions 
        if(real_parts(jwd,p_index)+1e-12.le.prad) then
!            write(6,*) 'Particle collision', integer_parts(jai,p_index)
          cl_c=cl_c + 1 !Count the number of collisions  
          !Coefficient of resitution
          res=e_max*exp(-(e_beta/real_parts(jar,p_index)))

          res=1
          !Determine which wall the paricle hit (1,2,3)
          colin=mymax2(real_parts(jwn:jwn+ndim-1,p_index),ndim)
!           write(6,*) 'Face', colin
          write(6,*) res,colin,'resitution and face'
          !Determine the face numbers
!           ee=integer_parts(je0,p_index)+1 !This needs to be tested, maybe
                                          !just use the coordinate data

          normc=real_parts(jx0+colin-1,p_index)!The actual collision position
!           write(6,*) 'position',(real_parts(jx0+i-1,p_index),i=1,ldim)
!           wallpt=real_parts(jwn+colin-1,p_index)!Not sure if I need this
!           write(6,*)'wall position',normc
          !Assign the wall value to the pt.
          !This needs some reworking
          diff=prad-real_parts(jwd,p_index)
          !!make this dependent on particle size if pd=1e-3 then use 1e4
          normc_check=ceiling(abs(normc)*100)/100
          write(6,*)normc_check,'normc_check'
!           if (diff.lt.0) then 
!           wallpt=abs(normc_check)-abs(diff)
!           else
          wallpt=normc_check-prad
!           endif
          wallpt=sign(wallpt,normc)
!           if(normc.lt.0) wallpt=-wallpt
!            write(6,*) real_parts(jx0+colin-1,p_index),'actualpt'
!            write(6,*) real_parts(jwd,p_index),'distance'
!            write(6,*) abs(prad-real_parts(jwd,p_index)),'diff'

          real_parts(jx0+colin-1,p_index)=wallpt
          write(6,*) wallpt,'wallpt' 
!           write(6,*) 'position2',(real_parts(jx0+i-1,p_index),i=1,ldim)
!           write(6,*)'actual wall position',wallpt
            !!!! Having some issues here !!!!!!
            !Store collision data 
          opp_dis=abs(wallpt-real_parts(jx1+colin-1,p_index))
          rslt_lngth=dist3d(real_parts(jwn,p_index)
     $      ,real_parts(jwn+1,p_index),real_parts(jwn+2,p_index)
     $      ,real_parts(jx1,p_index),real_parts(jx1+1,p_index)
     $      ,real_parts(jx1+2,p_index))
            col_theta=opp_dis/rslt_lngth


            !Save collision data
            e_col=integer_parts(je0,p_index)+1 
            eg_col =lglel(e_col)!Particle local el number
!            call copy(col_temp_r(1:3,cl_c),
!      $       real_parts(jwn:jwn+ndim-1,p_index),ndim)
           call copy(col_temp_r(1:3,cl_c),
     $       real_parts(jx:jz,p_index),ndim)
           call copy(col_temp_r(4:6,cl_c),
     $       real_parts(jv0:jv0+ndim-1,p_index),ndim)
            col_temp_r(7,cl_c)=col_theta !Collision angle
            col_temp_r(8,cl_c)=0.0
            col_temp_r(9,cl_c)=time
            
            col_temp_i(1,cl_c)=eg_col
            col_temp_i(2,cl_c)=integer_parts(jgn,p_index) !Particle GRP
            col_temp_i(3,cl_c)=integer_parts(jai,p_index) !Particle ID
            col_temp_i(4,cl_c)=colin !Now the coord direction ! Collision ID update after transfer
            col_temp_i(5,cl_c)=nid ! real procs number
!             col_temp_i(6,cl_c)=new_nid ! collection procs number
          write(6,*)'old vel',real_parts(jv0+colin-1,p_index)
            !!!Update the velocity
            real_parts(jv0+colin-1,p_index)=
     $     -res*real_parts(jv0+colin-1,p_index)
            real_parts(jv1+colin-1,p_index)=
     $     -real_parts(jv1+colin-1,p_index)
            real_parts(jv2+colin-1,p_index)=
     $     -real_parts(jv2+colin-1,p_index)
            real_parts(jv3+colin-1,p_index)=
     $     -real_parts(jv3+colin-1,p_index)
             !Update the acceleration terms
            real_parts(ja1+colin-1,p_index)=
     $     -real_parts(ja1+colin-1,p_index)
            real_parts(ja2+colin-1,p_index)=
     $     -real_parts(ja2+colin-1,p_index)
            real_parts(ja3+colin-1,p_index)=
     $     -real_parts(ja3+colin-1,p_index)

           write(6,*)'new vel',real_parts(jv0+colin-1,p_index)
           endif
      
      if (cl_c.eq.2*lhis.or.lastep.eq.1) then    
          call data_out (col_temp_r,rc_index,col_temp_i
     $   ,ic_index,cl_c,col_position)
           cl_c=0
       endif 


      return 
      end
c---------------------------------------------------------------------------------- 
      subroutine particle_position (real_parts,num_reals,
     $   integer_parts,num_ints,num_total,p_index)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      include 'RECYCLE'
      integer num_reals,num_ints,num_total,p_index
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      integer ee
      real eul_rpart,lag_rpart
      integer eul_ipart,lag_ipart
      common /euler_lagrang_stats_r/  eul_rpart(27*ldim+1,lhis),
     $                                lag_rpart(27*ldim+1,lhis)    
      common /eulerian__lagrang_stats_i/  eul_ipart(9,lhis),
     $                                    lag_ipart(9,lhis)
      integer eul_out,lag_out 
      save eul_out,lag_out 
      data eul_out ,lag_out /0,0/
      real pst_dwn        !Position to collect data from downstream
      save pst_dwn
      data pst_dwn /-1.5/
      integer icalld
      save    icalld
      data    icalld /0/
      character*12 eul_position (2) 
      save eul_position
      data eul_position /'eulerian','lagrangn' /
      integer jx0,jy0,jz0
      real prad
      jx0=jx
      jy0=jx0+1
      jz0=jy0+1
      prad=real_parts(jpd,p_index)*0.5 !Particle radius 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!             Recycle portion             !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (icalld.eq.0) then 
        call rzero(eul_rpart,num_reals*(lhis))  
        call izero(eul_ipart,num_ints*(lhis)) 
         icalld=icalld+1
      endif
     
      if(real_parts(jx0,p_index)+prad.gt.-recycle_pos
     $     .and.abs(integer_parts(jcl,p_index)).eq.1) then
        !! Create clone
        num_total=num_total+1
        eul_out=eul_out+1
!         write(6,*)'upstream'
!         write(6,*) 'made it here',eul_out,p_index,nid
        !! Index how many clones there are
        integer_parts(jcp,p_index)=integer_parts(jcp,p_index)+1
        
        !! Copy info to euler_stats
        call copy(eul_rpart(1,eul_out),
     $     real_parts(1,p_index),24*ldim+3)
        call icopy(eul_ipart(1,eul_out),
     $     integer_parts(1,p_index),9)
        
        !! Add new particle
        call copy(real_parts(1,num_total),
     $     real_parts(1,p_index),24*ldim+3)
        call icopy(integer_parts(1,num_total),
     $     integer_parts(1,p_index),9)
        write(6,*) 'added new particle',num_total
         
        integer_parts(jcl,num_total)=3


        real_parts(jx0,p_index)=-upstream_pos
     $   +real_parts(jpd,p_index)
        real_parts(jx1,p_index)=-upstream_pos
     $   -abs(real_parts(jx1,p_index)) 
        real_parts(jx2,p_index)=-upstream_pos
     $   -abs(real_parts(jx2,p_index)) 
        real_parts(jx3,p_index)=-upstream_pos
     $   -abs(real_parts(jx3,p_index))
!         Check this to make sure the particles are sent to the right element
!         do ee=1,yzels
!             if (integer_parts(je0,p_index)+1.eq.ee) then
!               integer_parts(jps,p_index)=recycl_map(2,ee)
!               integer_parts(je0,p_index)=recycl_map(4,ee)-1
!             endif  

!         enddo

      endif 
          

 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(real_parts(jy0,p_index)-prad.lt.-outlet_pos) then

        write(6,*)'particle ',integer_parts(jai,p_index),
     $    'has left the domain',nid
        write(6,*)'deleted paritcle', p_index,'new number of particles'
     $     ,num_total-1  
        if (p_index+1.le.num_total) then 
          call copy(real_parts(1,p_index),
     $       real_parts(1,num_total),num_reals)
          call icopy(integer_parts(1,p_index),
     $       integer_parts(1,num_total),num_ints)
        endif
        
        p_index=p_index-1
        num_total=num_total-1
        !write(6,*)real_parts(jy0,i)+prad,integer_parts(jai,i),integer_parts(jgn,i),'ID,group'
      endif
      
      !!Data collection location: Downstream
      !!redo this, need before and after
      !! Maybe move this up before the exit condition
      if(real_parts(jy0,p_index)-prad.le.pst_dwn
     $   .and.real_parts(jx1+1,p_index)-prad.ge.pst_dwn
     $   .and.istep.gt.3) then 
         !! Copy info to euler_stats
!         write(6,*) 'made it inside',eul_out,p_index,nid
!         write(6,*) 'position', real_parts(jy0,p_index)
!         write(6,*) 'prior', real_parts(jx1+1,p_index) 
!         write(6,*)'downstream' 
        eul_out=eul_out+1
        call copy(eul_rpart(1,eul_out),
     $     real_parts(1,p_index),num_reals)
        call icopy(eul_ipart(1,eul_out),
     $     integer_parts(1,p_index),num_ints)  
        
      endif
      
       
       
      if(integer_parts(jcl,p_index).lt.0) then 
        lag_out=lag_out+1
!         if(lag_out.eq.1) write(6,*)'cont.'
        call copy(lag_rpart(1,lag_out),
     $     real_parts(1,p_index),num_reals)
        call icopy(lag_ipart(1,lag_out),
     $     integer_parts(1,p_index),num_ints)  
      endif

      !!Check to see if data needs to be dumped
       if (eul_out.eq.lhis.or.lastep.eq.1) then  
          call data_out (eul_rpart,num_reals,eul_ipart
     $   ,num_ints,eul_out,eul_position(1))
           eul_out=0
       endif 
       
       if (lag_out.eq.lhis.or.lastep.eq.1) then    
          call data_out (lag_rpart,num_reals,lag_ipart
     $   ,num_ints,lag_out,eul_position(2))
           lag_out=0
       endif 


      return
      end
c------------------------------------------------------------------------------------------------
      subroutine data_out (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total,eul_name)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      integer icallp,rstpart
      save    icallp,rstpart
      data    icallp,rstpart /0,0/
      character*12 eul_name
      character*128 fname
      logical exist
      character*128 FRMT
      integer i,k
!       FRMT='(1p75e17.9,9I9)'
      if(num_total.gt.0) then  
        write(fname,1) eul_name, nid+1
 1      format('TIME_SERIES/',A8,'_',I4.4,'.log')
!         write(6,*) eul_name
!         write(6,*) fname
        write(FRMT,2) num_reals,num_ints
 2      format('(1p',I2,'e17.9,',I1,'I9)')
!         write(6,*) FRMT
!  3      format(I9,'I9')
        inquire(file=fname, exist=exist)
        if (exist) then 
          open(unit=74,file=fname,status='old', position='append') 
        else
          open(unit=74,file=fname,status='new')
        endif   
        !write(74,'(I8)') num_total
        write(74,FMT=FRMT)((real_parts(k,i),k=1,num_reals) ! xyz particles    
!      $     ,(real_parts(k,i),k=jv0,jv0+2)     ! velocity particles
!      $     ,(real_parts(k,i),k=ju0,ju0+2)
!      $     ,(real_parts(k,i),k=ja0,ja0+2)
     $     ,(integer_parts(k,i),k=1,num_ints),i=1,num_total) ! ID and Group number
!   4      format('thereals',2I9)
        close(74)
        write(6,*) eul_name  
      endif 
      return
      end

c------------------------------------------------------------------------------------------------
!       subroutine collision_detection
!       include 'SIZE'
!       include 'TOTAL'
!       include 'PARTICLES'
!       include 'RECYCLE'
!       include 'COLLISIONS' 
      
!       integer runtal !I need to intailize this somewhere
!       integer cl_up,cl_out         
      
!                !!!!!!Collision Detection!!!!!!!
!       num_col=cl_c
!       call igop(num_col,wrk_c,'+  ',1)
!       if (cl_c.gt.0.and.cl_c_up.ne.cl_c) then 
!         write(6,*) num_col,cl_c,nid,'number of collisions' 
!         cl_c_up=cl_c
!       endif    
!       if (num_col.ge.lhis.or.lastep.eq.1.or.istep.eq.nsteps) then 
!         ! call exitt 
!         !do i=1,cl_c
!         !write(6,*) col_temp_i(6,i),'ohoh',nid,cl_c
!         !enddo 
!         jps = 6     ! Pointer to temporary proc id for swapping
!         do i=1,n        ! Can't use jpt because it messes up particle info
!           col_temp_i(jps,i) = col_temp_i(6,i)
!         enddo    
               
!         call fgslib_crystal_tuple_transfer(col_handle,cl_c,2*lhis
!      $       ,col_temp_i,ic_index,col_temp_l,nl,col_temp_r,rc_index,jps)
!         call copy(col_data_r,col_temp_r,rc_index*cl_c)
!         call icopy(col_data_i,col_temp_i,ic_index*cl_c)
!          !call nekgsync
!          !   do i=1,cl_c
!          !   write(6,*) col_temp_i(6,i),'ohoh+after',nid,cl_c
!          !  enddo     
!         call rzero(col_temp_r,rc_index*(2*lhis)) 
!         call izero(col_temp_i,ic_index*(2*lhis)) 

!         cl_c=0
!         cl_c_up=cl_c
!         if (nid.eq.new_nid) then 
!           cl_out=num_col 
!           do i=1,cl_out
!             col_data_i(4,i)=runtal+i
!             col_data_i(6,i)=nid
!           !    write(6,*) col_data_i(6,i), nid,'index'
!           enddo
!         endif 
!         runtal=runtal+num_col
!         !  call exitt
!         ! write(6,*)nsteps,'nsteps'
!         ! write(6,*)lastep,'lastep'
!         if (new_nid+1.eq.np.or.lastep.eq.1.or.istep.eq.nsteps) then  
!           call collision_out(col_data_r,rc_index,
!      $            col_data_i,ic_index,cl_out)
!           !do i=1,cl_out
!           !write(6,*) col_data_i(4,i),'yupyup'
!           !enddo  
!           !write(6,*) runtal,'runtal_yo',nid
!           !write(6,*) cl_out,'cl_out_yo',nid
!           !write(6,*) new_nid,'new_nid_yo',nid


!           new_nid=-1
!           runtal =0 
!           call rzero(col_data_r,rc_index*(2*lhis)) 
!           call izero(col_data_i,ic_index*(2*lhis))
!         endif  
!         new_nid=new_nid+1

!       endif
!        !call igop(runtal,wrk_tal,'+  ',1)
!       return
!       end 
! c-----------------------------------------------------------------------
!       subroutine collision_out (rxyz,rdex,ixyz,idex,npart)
!       include 'SIZE'
!       include 'TOTAL'

!       integer ixyz(6,2*lhis),rdex,idex,npart
!       real rxyz(9,2*lhis)
!       integer cllsn_id(2*lhis)

!       common /r_col/ x_tmp(9,2*lhis),work_r(9,2*lhis)
!       common /i_col/ id_tmp(6,2*lhis),work_i(6,2*lhis)
!       character*128 fname

!       integer icalld
!       save    icalld
!       data    icalld  /0/

!       icalld = icalld+1
!       if (nid.eq.0) then
!         write(fname,1) icalld
!  1      format('collision',i5.5,'.3D')
!         open(unit=72,file=fname)
!         write(6,*) 'collision out: ', icalld 
!       endif
!       do i=1,npart
!       cllsn_id(i)=ixyz(4,i)
! !        write(6,*)ixyz(6,i),nid, 'new_nid'
!       enddo 
!       min_points = iglmin(cllsn_id, npart)
!       max_points = iglmax(cllsn_id,npart)
!       n_active   = max_points - min_points + 1
!       npass = n_active / lhis
!       if (n_active.gt.npass*lhis) npass = npass+1
!       ilast=min_points-1
!       if (nid.eq.0) write(72,'(I8)') max_points
!       do ipass=1,npass

!         mpart = min(lhis,max_points-ilast)
!         i0    = ilast
!         i1    = i0 + mpart
!         ilast = i1

!         call rzero(x_tmp,rdex*(2*lhis))
!         call izero(id_tmp,idex*(2*lhis))
!         do ii=1,npart
!           if (i0.lt.ixyz(4,ii).and.ixyz(4,ii).le.i1) then
!             i = ixyz(4,ii)-i0
!             call copy(x_tmp(1,i),rxyz(1,ii),rdex)  ! Real parts
!             call icopy(id_tmp(1,i),ixyz(1,ii),idex) ! Integer parts
                    
!           endif
!         enddo

!         call gop(x_tmp,work_r,'+  ',rdex*(2*lhis))
!         call igop(id_tmp,work_i,'+  ',idex*(2*lhis))
! !       if (nid.eq.0) write(72,'(I8)') max_points 
!         if (nid.eq.0) write(72,2)((x_tmp(k,i),k=1,rdex),
!      $     (id_tmp(k,i),k=1,idex),i=1,mpart)
!  2      format(1p9e17.9,6I8)

!       enddo
!       if (nid.eq.0) close(72)

!       return
!       end

c-------------------------------------------------------------------------
!       subroutine get_bdf_ext_coefs(beta,alpha,times)

!       include 'SIZE'
!       include 'TOTAL'

!       real beta(0:3),alpha(0:3),times(0:3)
!       real c(0:8)

!       integer ilast,ncoef
!       save    ilast,ncoef
!       data    ilast,ncoef / -9 , 0 /

!       do i=3,1,-1
!          times(i)=times(i-1)
!       enddo
!       times(0) = time

!       call rzero(beta ,4)
!       call rzero(alpha,4)
!       if (istep.ne.ilast) then
!          ilast = istep
!          ncoef = ncoef + 1
!          ncoef = min(ncoef,3) ! Maximum 3rd order in time
!       endif
!       ncoefm1 = ncoef - 1

!       call fd_weights_full(times(0),times(1),ncoefm1,0,alpha(1))
!       call fd_weights_full(times(0),times(0),ncoef,1,c)
!       do j=0,ncoef
!          beta(j) = c(ncoef+1+j)
!       enddo
!       do j=1,ncoef
!          beta(j) = -beta(j)  ! Change sign, for convenience
!       enddo
!       return
!       end
c----------------------------------------------------------------------
      subroutine particle_AB3_coef(c)
      !Adam Bashforth Coefficients 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      integer icalld
      save    icalld
      data    icalld  /0/
      real c(3)
      if (.not.intpart) icalld=2  ! For restart
! !        if(nid.eq.0) write(6,*) icalld, 'AB3 Coefficients' 
      if (icalld.eq.0) then      ! AB1(Euler Method)
         c(1) = 1.
         c(2) = 0.
         c(3) = 0.
         icalld=icalld+1
      elseif (icalld.eq.1) then  ! AB2
         c(1) = 3
         c(2) = -1.
         c(3) = 0
         c(1) = c(1)/2.
         c(2) = c(2)/2.
         icalld=icalld+1
      else                      ! AB3
         c(1) = 23.
         c(2) = -16.
         c(3) = 5.
         c(1) = c(1)/12.
         c(2) = c(2)/12.
         c(3) = c(3)/12
      endif
      
      return
      end
! c----------------------------------------------------------------------
      integer function mymax(VEC,VEC2,n)
      implicit none  
      REAL VEC(1), VEC2(1)
      INTEGER VA,n,i
      REAL TMAX,TMAX2,RR
      TMAX =-99.0E20
      VA=0
      do i=1,n
         
         RR=abs(VEC(i)-VEC2(i))
         TMAX2=TMAX         
         TMAX = MAX(TMAX,RR)
         if(TMAX.gt.TMAX2) VA=i
      enddo
      mymax = VA
      return
      END

c----------------------------------------------------------------------
      integer function mymax2(VEC,n)
      implicit none  
      REAL VEC(1)
      INTEGER VA,n,i
      REAL TMAX,TMAX2,RR
      TMAX =-99.0E20
      VA=0
      do i=1,n
         
         RR=abs(VEC(i))

         TMAX2=TMAX         
         TMAX = MAX(TMAX,RR)
         if(TMAX.gt.TMAX2) VA=i
!          write(6,*),VEC(i),i, VA,TMAX,'values'
      enddo
      mymax2 = VA
      return
      END
c----------------------------------------------------------------------- 
      subroutine Force (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 

      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
      integer nnp
      real visc,rho_f,p_Reynolds,Accel,rho_p,mat_div
      real cm,g(ldim),cd_r,Stokes,Stokes_mean,FG,FD,FS,FL,FP,Term_v
      integer icalld
      real vrel(ldim),crsscrl
      save icalld
      data icalld /0/
      save g,rho_f,cm,vrel
      data g /0.0,-9.81,0.0/
!       data g /0.0,0.0,0.0/
      data cm /0.5/
!       integer dij_cons(3,3),h
!       save dij_cons 
!       data dij_cons / 1,2,5
!      $               ,4,2,6
!      $               ,5,6,4/
      real rel_vel,diju
      real d_ratio,dvdt,g_mag,pv,u_mag,v_mag,vu_mag,coeff_D
      integer i,j,k,ii

      visc=param(02)
!       KK=2.594
      pi    = 4.*atan(1.0)
      rho_f=param(01)
      
      call vel_gradient(real_parts,num_reals,
     &   integer_parts,num_ints,num_total)
      call dudt(real_parts,num_reals,
     &   integer_parts,num_ints,num_total)
      call lift_var(real_parts,num_reals,
     &   integer_parts,num_ints,num_total)
      
      do i=1,num_total
        k=0
        crsscrl=0
        FG=0.0  
        FD=FG
        FS=FG
        FL=FS
        FP=FL
        v_mag=sqrt((real_parts(jv0,i))**2+(real_parts(jv0+1,i)**2)
     &   +(real_parts(jv0+2,i)**2))
        u_mag=sqrt((real_parts(ju0,i)**2)+(real_parts(ju0+1,i)**2)
     &   +(real_parts(ju0+2,i)**2))
        
        vu_mag=0
        do ii=0,ndim-1
          vrel(ii+1)=(real_parts(jv0+ii,i)-real_parts(ju0+ii,i))
          vu_mag=vu_mag+(vrel(ii+1)**2) 
        enddo
        vu_mag=sqrt(vu_mag)

        p_Reynolds=abs(vu_mag)*(real_parts(jpd,i))
        p_Reynolds=p_Reynolds/visc
             
        pv=(4./3.)*pi*((real_parts(jpd,i)/2.)**3)
          
        rho_p=real_parts(jrh,i)  
          
        d_ratio=(rho_p/rho_f)  
         
        cd_r=coeff_D(p_Reynolds)
        
        g_mag=sqrt(g(1)**2+g(2)**2+g(3)**2) 
        Term_v=((4/3)*abs(rho_p-rho_f)*real_parts(jpd,i)*g_mag)
        Term_v=sqrt(Term_v/(rho_f*cd_r))
   

        real_parts(jar,i)=((rho_p+cm*rho_f)*v_mag*real_parts(jpd,i))
     $                    /(9*visc)
        Stokes=((rho_p+cm*rho_f)*real_parts(jpd,i)**2)/(18*visc)

        Stokes_mean =real_parts(jar,i)
!         Term_v_ex=Stokes*g_mag*(1-(rho_f/rho_p))
!       !! Exact solution for particle in uniform flow
!           if (icalld.eq.0) then 
!             Stokes_ex(i)=Stokes
!           endif    
!         
         nnp=integer_parts(jai,i)
        if(nnp.eq.1)write(6,*)'Fluid viscosity',visc 
        if(nnp.eq.1)write(6,*)'Particle diameter',real_parts(jpd,i)
        if(nnp.eq.1)write(6,*)'Velocity relative:',vu_mag
        if(nnp.eq.1)write(6,*)'velocity P:',v_mag
        if(nnp.eq.1)write(6,*)'velocity px :',real_parts(jv0,i)
        if(nnp.eq.1)write(6,*)'velocity py :',real_parts(jv0+1,i)
        if(nnp.eq.1)write(6,*)'velocity pz :',real_parts(jv0+2,i)
        if(nnp.eq.1)write(6,*)'velocity F :',u_mag
        if(nnp.eq.1)write(6,*)'velocity fx :',real_parts(ju0,i)
        if(nnp.eq.1)write(6,*)'velocity fy :',real_parts(ju0+1,i)
        if(nnp.eq.1)write(6,*)'velocity fz :',real_parts(ju0+2,i)

        if(nnp.eq.1)write(6,*)'Particle Reynolds is :',p_Reynolds  
        if(nnp.eq.1)write(6,*)'Particle volume :',pv   
        if(nnp.eq.1)write(6,*)'Particle density is :',rho_p 
        if(nnp.eq.1)write(6,*)'Particle density ratio is :',d_ratio
        if(nnp.eq.1)write(6,*)'CD :',cd_r 
        if(nnp.eq.1)write(6,*) 'Stokes_mean number: ',Stokes_mean
        if(nnp.eq.1)write(6,*) 'Stokes number: ',Stokes
        if(nnp.eq.1)write(6,*) 'Term_velocity',Term_v
!         if(nnp.eq.1)write(6,*) 'Term_velocity_ex',Term_v_ex
        do j=0,ndim-1 
          FG=pv*g(j+1)*(rho_p - rho_f)
          if(nnp.eq.1)write(6,*)'Force acting on particle FG', FG, j+1 
          FD=-(real_parts(jpd,i)**2)*rho_f*cd_r*(1./8.0)
          FD= pi*FD*(abs(real_parts(jv0+j,i)-real_parts(ju0+j,i))
     &           *(real_parts(jv0+j,i)-real_parts(ju0+j,i)))
!           if(j.eq.1) real_parts(jgu,i)=FD
!           if(j.eq.1) real_parts(jgv,i)=FG

          if(nnp.eq.1)write(6,*)'Force acting on particle FD', FD, j+1
          dvdt=real_parts(jdt+j,i)
          mat_div=dvdt+real_parts(ju0,i)*real_parts(jgu+k,i)
     &         + real_parts(ju0+1,i)*real_parts(jgu+1+k,i)
     &         + real_parts(ju0+2,i)*real_parts(jgu+2+k,i)
          FS=rho_f*pv*(cm*dvdt+mat_div)
!           if(nnp.eq.1)write(6,*)'dudt',dvdt , j+1
!           if(nnp.eq.1)write(6,*)'dudx',real_parts(jgu+k,i) , j+1
!           if(nnp.eq.1)write(6,*)'dvdx',real_parts(jgu+1+k,i) , j+1
!           if(nnp.eq.1)write(6,*)'dwdx',real_parts(jgu+2+k,i) , j+1
!           if(nnp.eq.1)write(6,*)'mat_div',mat_div , j+1
          if(nnp.eq.1)write(6,*)'Force acting on particle FS', FS, j+1
          
          crsscrl=vrel(1)*real_parts(jgu+k,i)
     &                    +vrel(2)*real_parts(jgv+k,i)
     &                    + vrel(3)*real_parts(jgw+k,i)

          crsscrl=crsscrl-(vrel(1)*real_parts(jgu+k,i)+
     &    vrel(2)*real_parts(jgu+1+k,i)+vrel(3)*real_parts(jgu+2+k,i))
          

          FL=1.6*real_parts(jpd,i)**2*((visc*rho_f)**(0.5))*crsscrl
          FL=FL/sqrt(abs(real_parts(jds+j,i))) 

          if(nnp.eq.1)write(6,*)'Force acting on particle FL', FL, j+1
          
          FP=FG+FD+FS+FL
!           FP=FG
!           if(j.eq.1) real_parts(jgw,i)=FP
          if(nnp.eq.1)write(6,*)'Total force FP', FP, j+1

          Accel=FP/((rho_p+rho_f*cm)*pv)
!           !!!Accel=FP/((rho_p*pv)*(1+cm/d_ratio))

          if(nnp.eq.1)write(6,*)'Particle Acceleration',Accel, j+1
          real_parts(ja0+j,i)=Accel
          k=k+ndim
        enddo
      enddo
!       if(icalld.eq.0) then
!         call copy (g_ex,g,ndim)
! !         call copy (Stokes_ex,Stokes,n)
!         icalld=1
!       endif  
      
      return
      end
       
c----------------------------------------------------------------------- 
      real function coeff_D(p_Reynolds)
      implicit none
      real p_Reynolds,tmp_cd
      
      tmp_cd =0
      if(p_Reynolds.lt.1) then
         tmp_cd=1+(3./16.)*p_Reynolds
          
      elseif(p_Reynolds.gt.1.and.p_Reynolds.lt.285) then
!         tmp_cd=1
        tmp_cd=1+(0.1935)*p_Reynolds**(0.6385)
      elseif(p_Reynolds.gt.285.and.p_Reynolds.lt.2000) then
        tmp_cd=1+(0.015)*p_Reynolds+(0.02283)*p_Reynolds**(0.427)
      elseif(p_Reynolds.gt.2000) then
        tmp_cd=(0.44)*(p_Reynolds/24.)
      endif       
          if(p_Reynolds.eq.0) then 
          coeff_D=0
         else      
          coeff_D=(24/p_Reynolds)*tmp_cd
!           coeff_D=(24/p_Reynolds)

         endif
      return
      end
c----------------------------------------------------------------------- 
      subroutine vel_gradient(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'   
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
!       real  wrk(lx1*ly1*lz1*lelt,ldim)

!       call rzero(wrk,nx1*ny1*nz1*nelv*ndim)
      !grad U
!       call gradm1(wrk(1,1),wrk(1,2),wrk(1,3),vx)
      call particle_interp_local (real_parts,num_reals,integer_parts,
     $   num_ints,num_total,u_grad,jgu,ndim)
!         do j=1,n
!           do i=0,ndim-1
!           write(6,*) real_parts(jgu+i,j),'position dudx',i+1 
!           enddo
!         enddo
      !grad V
!       call rzero(wrk,nx1*ny1*nz1*nelt*ndim)
!       call gradm1(wrk(1,1),wrk(1,2),wrk(1,3),vy)
      call particle_interp_local (real_parts,num_reals,integer_parts,
     $   num_ints,num_total,v_grad,jgv,ndim)
!         do j=1,n
!           do i=0,ndim-1
!           write(6,*) real_parts(jgv+i,j),'position dvdx',i+1 
!           enddo
!         enddo
      !grad W
!       call rzero(wrk,nx1*ny1*nz1*nelt*ndim)
!       call gradm1(wrk(1,1),wrk(1,2),wrk(1,3),vz)
      call particle_interp_local (real_parts,num_reals,integer_parts,
     $   num_ints,num_total,w_grad,jgw,ndim)
!         do j=1,n
!           do i=0,ndim-1
!           write(6,*) real_parts(jgw+i,j),'position dwdx',i+1 
!           enddo
!         enddo       
     
      return
      end
! c-----------------------------------------------------------------------
      subroutine dudt(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      real  wrk1(lx1*ly1*lz1*lelt,ldim)
      real u1tmp(ldim,num_total),u2tmp(ldim,num_total)
      integer i,j
      real v,v1,v2  
       
       !t^n-1
       call opcopy(wrk1(1,1),wrk1(1,2),wrk1(1,3),VXLAG(1,1,1,1,1)
     &                   , VYLAG(1,1,1,1,1)
     &                   , VZLAG(1,1,1,1,1))

       call particle_interp_local (real_parts,num_reals,integer_parts
     &   ,num_ints,num_total,wrk1,jdt,ndim)
        
       do i=1,num_total 
            do j=1,ndim
             u1tmp(j,i)=real_parts(jdt+j-1,i)
            enddo
       enddo
       
       !t^n-2
       call opcopy(wrk1(1,1),wrk1(1,2),wrk1(1,3),VXLAG(1,1,1,1,2)
     &                   , VYLAG(1,1,1,1,2)
     &                   , VZLAG(1,1,1,1,2))

       call particle_interp_local (real_parts,num_reals,integer_parts
     &   ,num_ints,num_total,wrk1,jdt,ndim)
       
       do i=1,num_total 
            do j=1,ndim
             u2tmp(j,i)=real_parts(jdt+j-1,i)
            enddo
       enddo

       do i=1,num_total 
            do j=0,ndim-1
              v=real_parts(ju0+j,i)
              v1=u1tmp(j+1,i) 
              v2=u2tmp(j+1,i)
!                if (j.eq.0) write(6,*) v,v1,v2,dt,'v time'
              real_parts(jdt+j,i)=((3*v)-(4*v1)+v2)/(2*dt)

            enddo
       enddo

       return
       end 

c-------------------------------------------------------------------------
      subroutine lift_var (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)

      !curl of the velocity field 
      call particle_interp_local (real_parts,num_reals,integer_parts
     &   ,num_ints,num_total,dcrl,jds,ndim)
      
!       !Interpolate SijSji 
!       call particle_interp_local (real_parts,num_reals,integer_parts
!      &   ,num_ints,num_total,dsij,jss,1)

      return
      end  
c-------------------------------------------------------------------------
      subroutine volume_stats (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      integer lxyze     
      parameter (lxyze=lx1*ly1*lz1*lelv) 
      real alpha(pgrp,lelv),beta(pgrp,lelv),gamma(pgrp,lelv),summ
      integer e,p_grp,p_el,nxyz,nxyze,gp_el,i,j,k,ii,jj,ifielt
      real pv,pd,total_con,rr
      save pv
      pi    = 4.*atan(1.0)
          
       !!!!This assumes tht all particles are local no global 
       !!!!opertations perfomed   

      nxyze=nx1*ny1*nz1*nelv
      nxyz=nx1*ny1*nz1
      !Need to zero the entire array
      call rzero(concen,lxyze*pgrp)
      call rzero(vpdavg,lxyze*pgrp*ndim)
      call rzero(updavg,lxyze*pgrp*ndim)
      
      if(nio.eq.0) write(6,*) 'Computing particle_volume statistics ...' 
        call rzero(alpha,pgrp*nelv)
        call rzero(beta,pgrp*nelv)
        call rzero(gamma,pgrp*nelv)
          do i = 1,num_total

             pv=(4./3.)*pi*((real_parts(jpd,i)/2.)**3)
             p_grp=integer_parts(jgn,i)
             p_el=integer_parts(je0,i)+1
             gp_el=lglel(p_el)

            gamma(p_grp,p_el)=gamma(p_grp,p_el)+1
            beta(p_grp,p_el)=1/gamma(p_grp,p_el)
            alpha(p_grp,p_el)=1-beta(p_grp,p_el) 

!             call cadd(concen(1,1,1,p_el,p_grp),(pv/volel(p_el)),nxyz)
            call cadd(concen(1,1,1,p_el,p_grp),1.0,nxyz) ! just count the number of paritcles
!             call cadd(concen(1,1,1,p_el,p_grp),1,nxyz)   
                do ii=0,ndim-1
                    jj=ii+1
                
                !!!!!Move this outside the loop !!!!!!!!!!!!!! 
                do k=1,nxyz! Maybe split up

                  
                  vpdavg(k,1,1,p_el,p_grp,jj)=
     &             vpdavg(k,1,1,p_el,p_grp,jj)*alpha(p_grp,p_el)
     &             +beta(p_grp,p_el)*real_parts(jv0+ii,i)

                  updavg(k,1,1,p_el,p_grp,jj)=
     &             updavg(k,1,1,p_el,p_grp,jj)*alpha(p_grp,p_el)
     &             +beta(p_grp,p_el)*real_parts(ju0+ii,i)
                  enddo
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                enddo 
             
          enddo   
              !!!!!Move this outside the loop !!!!!!!!!!!!!! 
!                 do k=2,nxyz! Maybe split up

                  
!                   vpdavg(k,1,1,p_el,p_grp,jj)=vpdavg(k,1,1,p_el,p_grp,jj)
    
!                   updavg(k,1,1,p_el,p_grp,jj)=
!      &             updavg(k,1,1,p_el,p_grp,jj)*alpha(p_grp,p_el)
!      &             +beta(p_grp,p_el)*real_parts(ju0+ii,i)
!                   enddo
!                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!                 enddo 

      
        ifielt = ifield
        ifield = 1
        
        
        !Note maybe I can use the dssum to add the contributions 
        !of particles if they are in one than one element
        do i=1,pgrp
         call col2(concen(1,1,1,1,i),bm1,nxyze)
         call dssum(concen(1,1,1,1,i),nx1,ny1,nz1)
         call col2(concen(1,1,1,1,i),binvm1,nxyze) 
          do j=1,ndim

         call col2(vpdavg(1,1,1,1,i,j),bm1,nxyze)
         call dssum(vpdavg(1,1,1,1,i,j),nx1,ny1,nz1)  
         call col2(vpdavg(1,1,1,1,i,j),binvm1,nxyze)

         call col2(updavg(1,1,1,1,i,j),bm1,nxyze)
         call dssum(updavg(1,1,1,1,i,j),nx1,ny1,nz1)  
         call col2(updavg(1,1,1,1,i,j),binvm1,nxyze)  
          enddo
        enddo  
         ifield = ifielt
       call avg_all_particle
      return
      end  
c----------------------------------------------------------------------
c-----------------------------------------------------------------------
      FUNCTION ran2(idum)
      implicit none  
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
c between 0.0 and 1.0 (exclusive of the endpoint values). 
c Call with idum a negative integer to initialize; thereafter, do not alter 
c idum between successive deviates in a sequence. RNMX should approximate the 
c largest floating value that is less than 1.
      INTEGER idum1,idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then 
         idum1=max(-idum,1) 
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1 
            if (idum1.lt.0) idum1=idum1+IM1 
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1) 
      endif
      k=idum1/IQ1 
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1 
      k=idum2/IQ2 
      idum2=IA2*(idum2-k*IQ2)-k*IR2 
      if (idum2.lt.0) idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1 
      if(iy.lt.1)iy=iy+IMM1 
      ran2=min(AM*iy,RNMX)
      return
      END
c-----------------------------------------------------------------------



