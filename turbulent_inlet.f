      !                                                           !
      !          Routine for recycled boundary condition          !
      !          Just place call set_inflow in user check         ! 
      !          To adjust the position of the recycle boundary   !
      !          change the value of post(downstream) and         !
      !          upstream_pos (upstream)                                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine set_inflow
      include 'SIZE'
      include 'TOTAL'

      common /ibc/ gsh_bc

      integer gsh_bc
      common /rbc/ ubc(lx1,ly1,lz1,lelt)
     $           , vbc(lx1,ly1,lz1,lelt)
     $           , wbc(lx1,ly1,lz1,lelt) 
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /c_is1/ glo_num(lx1*ly1*lz1*lelv)
!       common /inlet_hndl/ il_handle 
      integer*8 glo_num,eg,iloc,locc,e
      integer*8 glo_loc,max_eg,n_in_pipe
      real nslab,elpul,tolin,scale
      integer icalld
      save    icalld
      data    icalld  /0/
c      data max_eg /0/
      
      integer f,ne_per_slab
      tolin=1e-8
      
!       nxyz  = nx1*ny1*nz1Fl 
          ! mesh parameters

!       nface = 2*ndim
       nxyz  = nx1*ny1*nz1
       n     = nxyz*nelv
     
!        ifxyo = .true.
!        if (istep.gt.iostep) ifxyo = .false.
     
! c     Connect slab next to the last with inflow (recycling method)
      if(icalld.eq.0) then
!         call crystal_setup (il_handle,nekcomm,np) 
        call get_points 
        call gs_setup(gsh_bc,glo_num,n,nekcomm,mp)
        icalld=1
      endif

      scale = 1.
      call v_mag(scale)
c  
!       
      call cmult2(ubc,vx,scale,n)
      call cmult2(vbc,vy,scale,n)
      call cmult2(wbc,vz,scale,n)
!        call outpost(ubc,vbc,wbc,pr,t,'ok1   ')

      call col2(ubc,v1mask,n)        ! Zero out except inflow
      call col2(vbc,v2mask,n)
      call col2(wbc,v3mask,n)
!       call outpost(ubc,vbc,wbc,pr,t,'ok2   ')

      call gs_op(gsh_bc,ubc,1,1,0)  ! 1 ==> +  ! uin(inflow) = vx(slab_k)
      call gs_op(gsh_bc,vbc,1,1,0)  ! 1 ==> +
      call gs_op(gsh_bc,wbc,1,1,0)  ! 1 ==> +
!        
!        call outpost(ubc,vbc,wbc,pr,t,'ok3   ')
!      call exitt   
      return
      end
      

c-----------------------------------------------------------------------
	    subroutine v_mag(scale)
      include 'SIZE'
      include 'TOTAL'
      include 'RECYCLE'
      integer*8 maxx,e,eg
      real mag, wt(nx1*ny1*nz1,nelt),scale
      post=0.5
      n= nx1*ny1*nz1*nelt
      nxyz  = nx1*ny1*nz1
        call rzero(wt,n)
      do e=1,nelv
         eg = lglel(e)         
!          if (eg.le.maxx)
          do i=1,nxyz
            
        if(xm1(i,1,1,e).le.-recycle_pos) then
              call copy(wt(i,e),bm1(i,1,1,e),1)         
!             write(6,*) xm1(i,1,1,e),'xm1'
!             write(6,*) ym1(i,1,1,e),'ym1'
!             write(6,*) zm1(i,1,1,e),'zm1'
           endif 
         enddo
      enddo
      vol_pipe = glsum(wt,n)
      
        ubar = glsc2(vx,wt,n)/vol_pipe
        vbar = glsc2(vy,wt,n)/vol_pipe
        wbar = glsc2(vz,wt,n)/vol_pipe
        mag= sqrt(ubar**2+vbar**2+wbar**2)
        if(istep.le.10.or.mod(istep,iostep).eq.0) then
        if (nid.eq.0)write(6,*) vol_pipe,'volume'
        if (nid.eq.0)write(6,*) ubar,'ubar'
        if (nid.eq.0)write(6,*) vbar,'vbar'
        if (nid.eq.0)write(6,*) wbar,'wbar'
!         if(istep.lt.100) mag=1
        if (nid.eq.0)write(6,*) mag,'velocity magnitude'
        if (nid.eq.0)write(6,*) 1/mag,'scale'
        endif
      
       
      scale=1
      scale=1/mag
      
      return 
      end	

c-----------------------------------------------------------------------
      subroutine get_points 
      include 'SIZE'
      include 'TOTAL'
      include 'RECYCLE'
!       include 'NEKUSE'  
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /c_is1/ glo_num(lx1*ly1*lz1*lelv)
      common /rbc_hndl/ up_handle,down_handle,rcylc_hand
!       integer lxyze,yzels,elcount
      !!!note that the space allocated for lxyze is based on the mesh !!!
      !!! This need to be upodated before running the simulation      !!!
!       parameter(yzels=4**2)
!       parameter(lxyze=ly1*lz1*yzels)
      real upstream_r(3,lxyze), downstream_r(3,lxyze)
      integer upstream_i(6,lxyze), downstream_i(6,lxyze)
!       integer recycl_map(4,yzels),wrk_map(4,yzels)
!        common /recycle_condition / recycl_map
!      $      downstream_r,downstream_i 

      
      integer*8 glo_num,glo_loc,before,after
      integer uc,dc,r_index,i_index,temp_e(ly1*lz1),temp_n(ly1*lz1)
      common /rbc_index/ uc,dc,r_index,i_index
      integer f,e,eg,gg,cc
      integer pgln,pln,pgel,plel,px,py,pz,rnid,stpr
      save pgln,pln,pgel,plel,rnid,stpr
      save px,py,pz
      data pgln,pln,pgel,plel,rnid,stpr /1,2,3,4,5,6/
      data px,py,pz /1,2,3/
      logical tempu_l !Dummy place holder for the crystal router 
      tolin=1e-4
      nface = 2*ndim
      nxyz  = nx1*ny1*nz1
      n     = nxyz*nelv
!       write(6,*) 'made it here' 
!       recycle_pos=0.5 ! Downstream
!       upstream_pos=1.5 ! Upstream
      uc=0
      dc=0
      nl=0
      i_index=6
      r_index=3
      
      
      call rzero(upstream_r,r_index*lxyze)
      call izero(upstream_i,i_index*lxyze)
      call rzero(downstream_r,r_index*lxyze)
      call izero(downstream_i,i_index*lxyze)  
      call izero (recycl_map,4*yzels)
!       call rzero(temp_e,ly1*)
!       call izero(temp_n,i_index*lxyze)  
      
      call crystal_setup (up_handle,nekcomm,np) 
      call crystal_setup (down_handle,nekcomm,np)
      call crystal_setup (rcylc_hand,nekcomm,np) 
       
      do e=1,nelv
         eg = lglel(e)
         do i=1,nxyz
            
            locc=i + nxyz*(e-1)
            glo_loc=i + nxyz*(eg-1)
c            gl_gl(glo_loc)=glo_loc 
            glo_num(locc) = glo_loc  !numbering all points with their
                                     !global number
            !! Upstream 
            eval_pt=xm1(i,1,1,e)                           
            if(abs(eval_pt+upstream_pos).lt.tolin) then
            uc=uc+1
            !Get all the upstream data points with their global location
            upstream_i(pgln,uc)=i
            upstream_i(pln,uc)=locc
            upstream_i(pgel,uc)=eg
            upstream_i(plel,uc)=e
            upstream_i(rnid,uc)=nid
            upstream_i(stpr,uc)=0   !setup to send to processor 0
            upstream_r(px,uc)=xm1(i,1,1,e)
            upstream_r(py,uc)=ym1(i,1,1,e)
            upstream_r(pz,uc)=zm1(i,1,1,e)
!              write(6,*)'u ',glo_loc,locc,i,e,eg,nid
            endif
            !! Downstream
!             if (nid.eq.0) write(6,*) recycle_pos, 'position to rec'
            eval_pt2=xm1(nx1,1,1,e)  
            if (abs(eval_pt+recycle_pos).lt.tolin) then
!      $       abs(eval_pt-eval_pt2).lt.tolin) then !right-hand rule
!              write(6,*)eval_pt,eval_pt2,'pts'
                  dc=dc+1  
                  downstream_i(pgln,dc)=i
                  downstream_i(pln,dc)=locc
                  downstream_i(pgel,dc)=eg
                  downstream_i(plel,dc)=e
                  downstream_i(rnid,dc)=nid
                  downstream_i(stpr,dc)=0  !setup to send to processor 0 
                  downstream_r(px,dc)=xm1(i,1,1,e)
                  downstream_r(py,dc)=ym1(i,1,1,e)
                  downstream_r(pz,dc)=zm1(i,1,1,e)
!             write(6,*)'d ',glo_loc,locc,i,e,eg,nid        
             endif   
!             
!             write(6,*)'global info',glo_num(locc),locc,eg,uc,nid
!             write(6,*)'element number', e,eg
          
        enddo
      enddo
      
      call nekgsync 
      if (uc.ne.0)  write(6,*) uc, nid, 'before_i'
      if (dc.ne.0)       write(6,*) dc, nid, 'before_k'
!       jps = plel+1     ! Pointer to temporary proc id for swapping
!          do i=1,uc        ! Can't use jpt because it messes up particle info
!             upstream_i(jps,i) = upstream_i(ll,i)
            

!          enddo
!       if (uc.ne.0) then
!         do i=1,uc
!       write(6,*) upstream_i(pgel,i),nid,i,'up global element number_b'
!       write(6,*) upstream_r(px,i),nid,i,'up global element number_b'
!         enddo
!       endif
!       if (dc.ne.0) then
!         do i=1,dc
!       write(6,*)downstream_i(pgel,i),nid,i,'dw global element number_b'
!         enddo
!       endif

      call crystal_tuple_transfer(up_handle,uc,lxyze
     $       ,upstream_i,i_index,upstream_l,nl,upstream_r,r_index,stpr)
      
      call crystal_tuple_transfer(down_handle,dc,lxyze
     $ ,downstream_i,i_index,downstream_l,nl,downstream_r,r_index,stpr)

!       if (uc.ne.0) then
!         do i=1,uc
!       write(6,*) upstream_i(pgel,i),nid,i,' up global element number'
!         enddo
!       endif
!       if (dc.ne.0) then
!         do i=1,dc
!       write(6,*) downstream_i(pgel,i),nid,i,'dw global element number'
!         enddo
!       endif

      call nekgsync 
      if (uc.ne.0)  write(6,*) uc, nid ,'after_i'
      if (dc.ne.0)  write(6,*) dc, nid ,'after_k'
!       call exitt
         
!       dc=0

      elcount=0
      do i=0,uc-ny1*nz1,ny1*nz1 ! upstream loop
!          write(6,*) i, 'counter'
             
             !upstream          
               
        do k=0,dc-ny1*nz1,ny1*nz1 ! downstream loop
         kk=0 

         do j=i+1,ny1*nz1+i !Upstream elements
          do l=k+1,ny1*nz1+k !Downstream elements
          
            if(abs(upstream_r(py,j)-
     &         downstream_r(py,l)).lt.tolin) then
            if(abs(upstream_r(pz,j)-
     &         downstream_r(pz,l)).lt.tolin) then
               kk=kk+1
!                write(6,*)'made it here',kk
               temp_n(kk)=downstream_i(pgln,l)
               temp_e(kk)=downstream_i(pgel,l)

             if (kk.eq.ny1*nz1) then 
             jj=0
             elcount=elcount+1
               recycl_map(1,elcount)=downstream_i(rnid,l)
               recycl_map(2,elcount)=upstream_i(rnid,j)
               recycl_map(3,elcount)=downstream_i(plel,l)
               recycl_map(4,elcount)=upstream_i(plel,j)
!                write(6,*) recycl_map(1,elcount),'proc',elcount
!                write(6,*) recycl_map(2,elcount),'proc_match',elcount
!                write(6,*) recycl_map(3,elcount),'elem',elcount
!                write(6,*) recycl_map(4,elcount),'elem_match',elcount

              do gg=i+1,ny1*nz1+i
                jj=jj+1
!                 write(6,*) upstream_i(pgel,gg), 'at',temp_e(jj),jj
                upstream_i(pgln,gg)=temp_n(jj)
                upstream_i(pgel,gg)=temp_e(jj)
              enddo
             
              goto 100 
             endif
            endif
            endif   
          enddo
         enddo       
        enddo                
 100   continue
      enddo
!       crystal_ituple_transfer(h, ituple,m,n,max, kp)
      call crystal_ituple_transfer(rcylc_hand,recycl_map,
     $   4,elcount,yzels,2)
      call crystal_tuple_transfer(up_handle,uc,lxyze
     $       ,upstream_i,i_index,upstream_l,nl,upstream_r,r_index,rnid)

      if (uc.ne.0)  write(6,*) uc, nid ,'back_i'
      if (dc.ne.0)  write(6,*) dc,nid ,'back_k'
      cc=0 
      do e=1,nelt
         eg = lglel(e)
         do i=1,nxyz
            locc=i + nxyz*(e-1)
            glo_loc=i + nxyz*(eg-1)
c            gl_gl(glo_loc)=glo_loc 
            
            if(abs(xm1(i,1,1,e)+upstream_pos).lt.tolin) then
              cc=cc+1
              before=glo_num(locc)
              after=upstream_i(pgln,cc)
              eg_new=upstream_i(pgel,cc)
              glo_num(locc) = after+nxyz*(eg_new-1)
!               write(6,*)'global info before',before
!               write(6,*)'global info after',glo_num(locc)

            endif 
         enddo
      enddo  

!        call igop(recycl_map,wrk_map,'+  ',4*yzels)
!       if (nid.eq.1) then 
!        do e=1,yzels
!         write(6,*)recycl_map(3,e),'mapping',nid,e
!        enddo 
!       endif
!        call exitt       
      
      return
      end
c----------------------------------------------------------------------------      
      
      
