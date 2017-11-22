c---------------------------------------------------------------------c
c             Particle Tracking Routine for Nek5000                   
c   Add "call particle_generator" to your .usr file        
c   in the userchk. To initialize particle location update            
c     particle_init subroutine                                        
c---------------------------------------------------------------------c
      subroutine set_part_pointers
      include 'SIZE'
      include 'PARTICLES' 
c     Minimal value of ni = 5
c     Minimal value of nr = 14*ndim + 1
      common  /iparti/ n,nr,ni
!       common /ptpointers/ jrc,jpt,je0,jps,jai,jgn,nai,jr,jd,jx,jy,jz,jx1
!      $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
!      $         ,ju0,ju1,ju2,ju3,jgu,jgv,jgw,jwd,jwn,jpd,jrh,jdt,jar,nar
      jrc = 1 ! Pointer to findpts return code
      jpt = 2 ! Pointer to findpts return processor id
      je0 = 3 ! Pointer to findpts return element id
      jps = 4 ! Pointer to proc id for data swap
      jai = 5 ! Pointer to auxiliary integers
      jgn = 6 ! Pointer to group number
      jcl = 7 ! Pointer to if the particle should be recycled 
      jcp = 8 ! Pointer to number of times ID has been copied

      nai = ni - (jcp-1)  ! Number of auxiliary integers
      if (nai.le.0) call exitti('Error in nai:$',ni)

      jr  = 1         ! Pointer to findpts return rst variables
      jd  = jr + ndim ! Pointer to findpts return distance
      jx  = jd + 1    ! Pointer to findpts input x value
      jy  = jx + 1    ! Pointer to findpts input y value
      jz  = jy + 1    ! Pointer to findpts input z value

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
      jwn=  jwd+ndim  ! Pointer to nearest wall coordinates
      jpd=  jwn+ndim  ! Pointer to particle diameter
      jrh=  jpd+1     ! Pointer to particle density 
      jdt=  jrh+1     ! Pointer to du/dt at particle location   
      jar = jdt+ndim     ! Pointer to auxiliary reals
                         !Right now used as used to used to store
                         !Particle ,stokes number 
      nar = nr - (jar-1)  ! Number of auxiliary reals

      if (nar.le.0) call exitti('Error in nar:$',nr)
      return
      end
c----------------------------------------------------------------------
      subroutine particle_generator ! Particle injection
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      include 'COLLISIONS'
!       include 'CTIMER'
      parameter (lr=24*ldim+3,li=9)
      common /rparts/ rpart(lr,lpart)
      common /iparts/ ipart(li,lpart)
      common /iparti/ n,nr,ni
!       real  xdrange(2,3) 
!       common /domainrange/ xdrange

      integer icalld,ipstep
      save    icalld
      data    icalld  /0/
      
      if(icalld.eq.0) then
        nr=lr
        ni=li
        icalld=1
        
        call rzero(rpart,lr*lpart)
        call izero(ipart,li*lpart)

        call set_part_pointers
      endif 
        
       call particle_post (rpart,nr,ipart,ni,n)
       
       call exitt
      
!       if(istep.ge.1000) then
!         if(mod(istep,iostep/1000).eq.0)then     
!           call volume_stats(rpart,nr,ipart,ni,n)
!         endif
!       endif
      
!       ipstep = iostep
!       call output_particles   (rpart,nr,ipart,ni,n,ipstep)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine particle_post(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES'
      include 'COLLISIONS'
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
     
      character*128 fname
      integer icalld,check_close
      save    icalld,check_close
      data    icalld,check_close /0,1/
      character*128 eul_name
      logical exist
      integer icount,actual
      real work(lx1,ly1,lz1,lelv)
      common /save_w/ work
      icount=lpart
      eul_name='collisns' 
        write(fname,1) eul_name, nid+1
 1      format('TIME_SERIES/',A8,'_',I4.4,'.log')
      inquire(file=fname, exist=exist)
      if(exist) then
        write(6,*) 'reading file: ',fname
        open(unit=98,file=fname)
      else
      write(6,*)'File does not exist'      
      endif
      
      if (exist) then 
!       do ipass=1,1
!       read(98,4,END=20)((col_temp_r(k,ii),k=1,rc_index)
!      $     ,(col_temp_i(k,ii),k=1,ic_index),ii=1,icount)
!  4      format(1p9e17.9,6I9)
!       endif
!         do ipass=1,2
        if (check_close.ne.0) then     
        do ii=1,icount
            read(98,4,END=20)(col_temp_r(k,ii),k=1,rc_index)
     $     ,(col_temp_i(k,ii),k=1,ic_index)
 4      format(1p9e17.9,6I9)
        enddo 
!    20   close(98,IOSTAT=check_close)
   20    write(6,*) check_close,nid+1,ii-1

!       if (nid.eq.0) write(6,*)'made it here again', ipass
          
!       if(nid.eq.0) write(6,*) k ,'total particle count'
      write(6,'(I12,A,I12)') ii-1, ' particles on proc: ', nid+1
      actual=ii-1
      
      do i=1,actual
        call copy(real_parts(jx:jz,i),
     $       col_temp_r(1:3,i),ndim)
        call copy(real_parts(jv0:jv0+ndim-1,i),
     $       col_temp_r(4:6,i),ndim)
        call copy(real_parts(jx1:jx1+ndim-1,i),
     $       col_temp_r(7:9,i),ndim)
        integer_parts(jgn,i)=col_temp_i(2,i)
        integer_parts(jai,i)=col_temp_i(3,i)
      enddo 


      call particle_interp_all (real_parts,num_reals,integer_parts
     $   ,num_ints,actual)
                
      call collision_count (real_parts,num_reals,integer_parts
     $   ,num_ints,actual)
      endif
!    20   close(98,IOSTAT=check_close)
!       write(6,*) check_close,nid+1
!       enddo
      endif 
    
       

      call nekgsync
!        close(98)

       ifvo=.false.
          ifto=.true.
          ifpo=.false.
          call outpost(vx,vy,vz,pr,work,'yup')
      return
      end
! c---------------------------------------------------------------------------
      subroutine particle_interp_all(real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)
      
!       real  wrk(lx1*ly1*lz1*lelt,ldim)  
!       common /outtmp/ wrk
      parameter(integer nmax=lpart) 
!       common /iv_intp/ ihandle_local,ihandle_remote,cr_handle
      parameter (lrf=24*ldim+3,lif=8+1)
      real               loc_rpts(lrf,lpart),rem_rpts(lrf,lpart)
      common /fptspartr/ loc_rpts,rem_rpts
      integer            loc_ipts(lif,lpart),rem_ipts(lif,lpart)
      integer            loc_ptsmap(lpart),rem_ptsmap(lpart)
      common /fptsparti/ loc_ipts,rem_ipts,n_loc_pts,n_rem_pts 
      common /mapss/     loc_ptsmap,rem_ptsmap
      integer icalld1,wrk_n
      save    icalld1
      data    icalld1 /0/
            
      logical partl         ! This is a dummy placeholder, used in cr()
      num_log = 0                ! No logicals exchanged
c      write(6,*) size(loc_ptsmap),'size of local maps' 
      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt
!       call opcopy(wrk(1,1),wrk(1,2),wrk(1,3),vx,vy,vz)
      
      if (icalld1.eq.0) then    ! interpolation setup !? intpts_done(ih_intp_v)?
        icalld1 = icalld1+1
        tolin  = 1.e-12
        if (wdsize.eq.4) tolin = 1.e-6
        call crystal_setup (cr_handle,nekcomm,np)
        call intpts_setup_part(tolin,ihandle_local)
c        call intpts_setup(tolin,ihandle_remote)
      endif

        if(nid.eq.0) write(6,*) 'call findpts'
         call findpts(ihandle_local,integer_parts(jrc,1),lif,
     &               integer_parts(jpt,1),lif,
     &               integer_parts(je0,1),lif,
     &               real_parts(jr, 1),lrf,
     &               real_parts(jd, 1),lrf,
     &               real_parts(jx, 1),lrf,
     &               real_parts(jy, 1),lrf,
     &               real_parts(jz, 1),lrf, num_total)
       
!         jps = jai-1     ! Pointer to temporary proc id for swapping
!          do i=1,num_total        ! Can't use jpt because it messes up particle info
!             integer_parts(jps,i) = integer_parts(jpt,i)
!          enddo
!          write(6,*) num_total
!          call crystal_tuple_transfer(cr_handle,num_total,lpart
!      $              , integer_parts,num_ints,partl
!      $              ,num_log,real_parts,num_reals,jps)
          

c        Sort by element number - for improved local-eval performance
         
!          call crystal_tuple_sort    (cr_handle,n 
!      $              , integer_parts,ni,partl,nl,real_parts,nr,je0,1)


c             Interpolate the particles locally          
!          call particle_interp_local (real_parts,num_reals,integer_parts
!      $       ,num_ints,num_total,wrk,ju0,ndim)
     
      return
      end
c------------------------------------------------------------------------- 
      subroutine intpts_setup_part(tolin,ih)
c
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

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
      call findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
      return
      end        
c------------------------------------------------------------------------
      subroutine particle_interp_local  (loc_r_data,nrf,loc_i_data
     &           ,nif,num_pts,wrk,pt,noflds)
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      real loc_r_data(nrf,*)
      integer loc_i_data(nif,*),pt
      real wrk(1)
      
      
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
          call findpts_eval_local(ihandle_local,
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
      subroutine collision_count (real_parts,num_reals,integer_parts
     $   ,num_ints,num_total)
      include 'SIZE'
      include 'TOTAL'
      include 'PARTICLES' 
      integer num_reals,num_ints,num_total
      real    real_parts(num_reals,num_total)
      integer integer_parts(num_ints,num_total)     
      parameter (lxyze=lx1*ly1*lz1*lelv) 
      real alpha(pgrp,lelv),beta(pgrp,lelv),gamma(pgrp,lelv)
      integer e,p_grp,p_el,nxyz,gp_el,f
      real pv,pd,total_con,rr
      save pv
      logical IFNORX,IFNORY,IFNORZ,test
      real work(lx1,ly1,lz1,lelv)
      common /save_w/ work
      integer icalld
      save icalld
      data icalld /0/

      pi    = 4.*atan(1.0)
      
       !!!!This assumes tht all particles are local no global 
       !!!!opertations perfomed   

      nxyze=nx1*ny1*nz1*nelv
      nxyz=nx1*ny1*nz1
!       !Need to zero the entire array
      if (icalld.eq.0) then 
            icalld=1
      call rzero(work,nxyze)
!       call rzero(concen,lxyze*pgrp)
      call rzero(vpdavg,lxyze*pgrp*ndim)
!       call rzero(updavg,lxyze*pgrp*ndim)
      endif
      if(nio.eq.0) write(6,*) 'Computing particle_volume statistics ...' 
        call rzero(alpha,pgrp*nelv)
        call rzero(beta,pgrp*nelv)
        call rzero(gamma,pgrp*nelv)
          do i = 1,num_total

! !              pv=(4./3.)*pi*((real_parts(jpd,i)/2.)**3)
             p_grp=integer_parts(jgn,i)
             p_el=integer_parts(je0,i)+1
! !              gp_el=lglel(p_el)

            gamma(p_grp,p_el)=gamma(p_grp,p_el)+1
            beta(p_grp,p_el)=1/gamma(p_grp,p_el)
            alpha(p_grp,p_el)=1-beta(p_grp,p_el) 
               xx=real_parts(jx,i)
               yy=real_parts(jy,i)
               zz=real_parts(jz,i)
               xx=NINT(xx*10)
               yy=NINT(yy*10)
               zz=NINT(zz*10)
!                write(6,*) xx, 'xx'
!                write(6,*) yy, 'yy'
!                write(6,*) zz, 'zz'

! !             call cadd(concen(1,1,1,p_el,p_grp),(pv/volel(p_el)),nxyz)
! !             call cadd(concen(1,1,1,p_el,p_grp),1.0,nxyz) ! just count the number of paritcles
! ! !             call cadd(concen(1,1,1,p_el,p_grp),1,nxyz)   
! ! !                 do ii=0,ndim-1
! ! !                     jj=ii+1
                
               do f=1,2*ndim
                  if (cbc(f,p_el,1).eq.'W  ') then 
               CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,f,p_el)
               if(xx.eq.abs(5)) then
                  test=IFNORX
                  ii=0
               elseif(yy.eq.abs(5)) then 
                  test=IFNORY 
                  ii=1
               elseif(zz.eq.abs(5)) then 
                  test=IFNORZ
                  ii=2
               endif
               if (test) then 
          v_mag=sqrt((real_parts(jv0,i))**2+(real_parts(jv0+1,i)**2)
     &   +(real_parts(jv0+2,i)**2))
                  do k=1,nxyz! Maybe split up
         
!                    vpdavg(k,1,1,p_el,p_grp,1)=
!      &             vpdavg(k,1,1,p_el,p_grp,1)*alpha(p_grp,p_el)
!      &             +beta(p_grp,p_el)*v_mag
                       vpdavg(k,1,1,p_el,p_grp,1)=
     &             vpdavg(k,1,1,p_el,p_grp,1)+1.0
                  enddo

       call facev(work,p_el,f,vpdavg(1,1,1,1,p_grp,1),nx1,ny1,nz1)
               endif


!                   updavg(k,1,1,p_el,p_grp,jj)=
!      &             updavg(k,1,1,p_el,p_grp,jj)*alpha(p_grp,p_el)
!      &             +beta(p_grp,p_el)*real_parts(ju0+ii,i)
!                   enddo
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                  endif 
                enddo 
             
          enddo   
!               !!!!!Move this outside the loop !!!!!!!!!!!!!! 
! !                 do k=2,nxyz! Maybe split up

                  
! !                   vpdavg(k,1,1,p_el,p_grp,jj)=vpdavg(k,1,1,p_el,p_grp,jj)
    
! !                   updavg(k,1,1,p_el,p_grp,jj)=
! !      &             updavg(k,1,1,p_el,p_grp,jj)*alpha(p_grp,p_el)
! !      &             +beta(p_grp,p_el)*real_parts(ju0+ii,i)
! !                   enddo
! !                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! !                 enddo 

           call nekgsync()
        ifielt = ifield
        ifield = 1
!           ifvo=.false.
!           ifto=.true.
!           ifpo=.false.
!           call outpost(vx,vy,vz,pr,work,'yup')
        
!         !Note maybe I can use the dssum to add the contributions 
!         !of particles if they are in one than one element
!         do i=1,pgrp
!          call col2(concen(1,1,1,1,i),bm1,nxyze)
!          call dssum(concen(1,1,1,1,i),nx1,ny1,nz1)
!          call col2(concen(1,1,1,1,i),binvm1,nxyze) 
!           do j=1,ndim

         call col2(work,bm1,nxyze)
         call dssum(work,nx1,ny1,nz1)  
         call col2(work,binvm1,nxyze)

!          call col2(updavg(1,1,1,1,i,j),bm1,nxyze)
!          call dssum(updavg(1,1,1,1,i,j),nx1,ny1,nz1)  
!          call col2(updavg(1,1,1,1,i,j),binvm1,nxyze)  
!           enddo
!         enddo  
         ifield = ifielt
!        call avg_all_particle
      return
      end
c----------------------------------------------------------------------
c-----------------------------------------------------------------------

