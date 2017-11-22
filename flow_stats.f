c-----------------------------------------------------------------------
      subroutine avg_all_flow
c
c     This routine computes running averages E(X),E(X^2),E(X*Y)
c     and outputs to avg*.fld*, rms*.fld*, and rm2*.fld* for all
c     fields.
c
c     E denotes the expected value operator and X,Y two
c     real valued random variables.
c
c     variances and covariances can be computed in a post-processing step:
c
c        var(X)   := E(X^X) - E(X)*E(X) 
c        cov(X,Y) := E(X*Y) - E(X)*E(Y)  
c
c     Note: The E-operator is linear, in the sense that the expected
c           value is given by E(X) = 1/N * sum[ E(X)_i ], where E(X)_i
c           is the expected value of the sub-ensemble i (i=1...N).
c
      include 'SIZE'  
      include 'TOTAL' 
      include 'AVG' !update this routine to indepedent of the this block
      include 'PARTICLES'
      real pa(lx1,ly1,lz1,lelt)
      real pb(lx1,ly1,lz1,lelt)
      real work1(lx1,ly1,lz1,lelt)
      real work2(lx1,ly1,lz1,lelt,ldim)
      logical ifverbose
      integer icalld
      save    icalld
      data    icalld  /0/

      if (ax1.ne.lx1 .or. ay1.ne.ly1 .or. az1.ne.lz1) then
         if(nid.eq.0) write(6,*)
     $     'ABORT: wrong size of ax1,ay1,az1 in avg_all(), check SIZE!'
         call exitt
      endif
      if (ax2.ne.lx2 .or. ay2.ne.ay2 .or. az2.ne.lz2) then
         if(nid.eq.0) write(6,*)
     $     'ABORT: wrong size of ax2,ay2,az2 in avg_all(), check SIZE!'
         call exitt
      endif
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelv
      nto2  = nx2*ny2*nz2*nelv
        
      ! initialization
      if (icalld.eq.0) then
         icalld = icalld + 1
         atime  = 0.
         timel  = time
         call rzero(work1,ntot)
         call rzero(work2,ntot*ldim)

         call rzero(uavg,ntot)
         call rzero(vavg,ntot)
         call rzero(wavg,ntot)
         call rzero(pravg,ntot)
        

         call rzero(urms,ntot)
         call rzero(vrms,ntot)
         call rzero(wrms,ntot)
         call rzero(prrms,ntot)
         
         ! Two point correlations
         call rzero(vwms,ntot)
         call rzero(wums,ntot)
         call rzero(uvms,ntot)
         
         call rzero(pums,ntot)
         call rzero(pvms,ntot)
         call rzero(pwms,ntot)
         ! For the pressure correlations
         call rzero(ppr,ntot)
         call rzero(pa,nxyz)
         call rzero(pb,nxyz)
         call rzero(pdudx,ntot*ldim)
         call rzero(pdvdx,ntot*ldim)
         call rzero(pdwdx,ntot*ldim)

         !Triple correlations
         call rzero(u2ums,ntot)
         call rzero(v2vms,ntot)
         call rzero(w2wms,ntot)
         call rzero(p2pms,ntot)
         
         call rzero(u2vms,ntot)
         call rzero(u2wms,ntot)

         call rzero(v2ums,ntot)
         call rzero(v2wms,ntot)

         call rzero(w2ums,ntot) 
         call rzero(w2vms,ntot)
         call rzero(uvwms,ntot)
         
         !Quadruple correlations
         call rzero(u4m,ntot) 
         call rzero(v4m,ntot) 
         call rzero(w4m,ntot) 
         call rzero(p4m,ntot) 

         !more Stuff

         call rzero(e11,ntot)
         call rzero(e22,ntot)
         call rzero(e33,ntot)

         call rzero(e12,ntot)
         call rzero(e13,ntot)
         call rzero(e23,ntot)

      endif

      dtime = time  - timel
      atime = atime + dtime
      
      ! dump freq
      iastep = param(68)
      if  (iastep.eq.0) iastep=param(15)   ! same as iostep
      if  (iastep.eq.0) iastep=500

      ifverbose=.false.
      if (istep.le.10) ifverbose=.true.
      if  (mod(istep,iastep).eq.0) ifverbose=.true.

      if (atime.ne.0..and.dtime.ne.0.) then
         if(nio.eq.0) write(6,*) 'Computing flow statistics ...'
         beta  = dtime/atime
         alpha = 1.-beta
         !map the pressure to the velocity mesh
         call mappr(ppr,pr,pa,pb)

         ! compute averages E(X)
         call myavg1    (uavg,vx,alpha,beta,ntot ,'um   ',ifverbose)
         call myavg1    (vavg,vy,alpha,beta,ntot ,'vm   ',ifverbose)
         call myavg1    (wavg,vz,alpha,beta,ntot ,'wm   ',ifverbose)
         call myavg1    (pravg,ppr,alpha,beta,ntot ,'prm  ',ifverbose)
        

         ! compute averages E(X^2) 
         call myavg2    (urms,vx,alpha,beta,ntot ,'ums  ',ifverbose)
         call myavg2    (vrms,vy,alpha,beta,ntot ,'vms  ',ifverbose)
         call myavg2    (wrms,vz,alpha,beta,ntot ,'wms  ',ifverbose)
         call myavg2    (prrms,ppr,alpha,beta,ntot ,'prms ',ifverbose)
         

         ! compute averages E(X*Y) (for now just for the velocities)
         call myavg3    (uvms,vx,vy,alpha,beta,ntot,'uvm  ',ifverbose)
         call myavg3    (vwms,vy,vz,alpha,beta,ntot,'vwm  ',ifverbose)
         call myavg3    (wums,vz,vx,alpha,beta,ntot,'wum  ',ifverbose)
         call col3(work1,vx,vy,ntot)
         call myavg3(uvwms,vz,work1,alpha,beta,ntot,'uvwm ',ifverbose)
         ! compute averages E(X*Y) including the pressure terms
         
         !<pdudx>t
         call myavg3    (pdudx(1,1,1,1,1),ppr,u_grad,alpha,beta,ntot 
     $    ,'pdudx',ifverbose)
         call myavg3    (pdudx(1,1,1,1,2),ppr,v_grad,alpha,beta,ntot 
     $    ,'pdudy',ifverbose)
         call myavg3    (pdudx(1,1,1,1,3),ppr,w_grad,alpha,beta,ntot 
     $    ,'pdudz',ifverbose)

         !<pdvdx>t
         call myavg3    (pdvdx(1,1,1,1,1),ppr,u_grad,alpha,beta,ntot 
     $    ,'pdvdx',ifverbose)
         call myavg3    (pdvdx(1,1,1,1,2),ppr,v_grad,alpha,beta,ntot 
     $    ,'pdvdy',ifverbose)
         call myavg3    (pdvdx(1,1,1,1,3),ppr,w_grad,alpha,beta,ntot 
     $    ,'pdvdz',ifverbose)
         
         !<pdwdx>t
         call myavg3    (pdwdx(1,1,1,1,1),ppr,u_grad,alpha,beta,ntot 
     $    ,'pdwdx',ifverbose)
         call myavg3    (pdwdx(1,1,1,1,2),ppr,v_grad,alpha,beta,ntot 
     $    ,'pdwdy',ifverbose)
         call myavg3    (pdwdx(1,1,1,1,3),ppr,w_grad,alpha,beta,ntot 
     $    ,'pdwdz',ifverbose)
         
         call myavg3    (pums,vx,ppr,alpha,beta,ntot,'pum  ',ifverbose)
         call myavg3    (pvms,vy,ppr,alpha,beta,ntot,'pvm  ',ifverbose)
         call myavg3    (pwms,vz,ppr,alpha,beta,ntot,'pwm  ',ifverbose)
         


         ! compute averages E(X^2*Y) for the velocity terms
         call myavg4    (u2ums,vx,vx,alpha,beta,ntot,'u2um ',ifverbose)
         call myavg4    (v2vms,vy,vy,alpha,beta,ntot,'v2vm ',ifverbose)
         call myavg4    (w2wms,vz,vz,alpha,beta,ntot,'w2wm ',ifverbose)
         call myavg4  (p2pms,ppr,ppr,alpha,beta,ntot,'p2pm ',ifverbose)
         call myavg4    (u2vms,vx,vy,alpha,beta,ntot,'u2vm ',ifverbose)
         call myavg4    (u2wms,vx,vz,alpha,beta,ntot,'u2wm ',ifverbose)

         call myavg4    (v2ums,vy,vx,alpha,beta,ntot,'v2um ',ifverbose)
         call myavg4    (v2wms,vy,vz,alpha,beta,ntot,'v2wm ',ifverbose)

         call myavg4    (w2ums,vz,vx,alpha,beta,ntot,'w2um ',ifverbose)
         call myavg4    (w2vms,vz,vy,alpha,beta,ntot,'w2vm ',ifverbose)
         


         ! Quad products
         ! Will update this syntax sucks, but for now this will do
         call col3(work1,vx,vx,ntot)
         call myavg4    (u4m,vx,work1,alpha,beta,ntot,'u4m ',ifverbose)
         call col3(work1,vy,vy,ntot)
         call myavg4    (v4m,vy,work1,alpha,beta,ntot,'v4m ',ifverbose)
         call col3(work1,vz,vz,ntot)
         call myavg4    (w4m,vz,work1,alpha,beta,ntot,'w4m ',ifverbose)
         call col3(work1,ppr,ppr,ntot)
         call myavg4    (p4m,ppr,work1,alpha,beta,ntot,'p4m ',ifverbose)

         !<e11>t : (du/dx)^2 + (du/dy)^2 + (du/dz)^2
         call col3(work2,u_grad,u_grad,ntot*ldim)
         call add4(work1,work2(1,1,1,1,1),work2(1,1,1,1,2),
     $     work2(1,1,1,1,3),ntot)
          
         call myavg1    (e11,work1,alpha,beta,ntot ,'e11  ',ifverbose)
         

         !<e22>t : (dv/dx)^2 + (dv/dy)^2 + (dv/dz)^2
         call col3(work2,v_grad,v_grad,ntot*ldim)
         call add4(work1,work2(1,1,1,1,1),work2(1,1,1,1,2),
     $     work2(1,1,1,1,3),ntot)
          
         call myavg1    (e22,work1,alpha,beta,ntot ,'e22  ',ifverbose)

         !<e33>t : (dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2
         call col3(work2,w_grad,w_grad,ntot*ldim)
         call add4(work1,work2(1,1,1,1,1),work2(1,1,1,1,2),
     $     work2(1,1,1,1,3),ntot)
          
         call myavg1    (e33,work1,alpha,beta,ntot ,'e33  ',ifverbose)
         
         !<e12>t: (du/dx)*(dv/dx) + (du/dy)*(dv/dy) + (du/dz)*(dv/dz)

         call col3(work2,u_grad,v_grad,ntot*ldim)
         call add4(work1,work2(1,1,1,1,1),work2(1,1,1,1,2),
     $     work2(1,1,1,1,3),ntot)
          
         call myavg1    (e12,work1,alpha,beta,ntot ,'e12  ',ifverbose)
         

         !<e13>t: (du/dx)*(dw/dx) + (du/dy)*(dw/dy) + (du/dz)*(dw/dz)

         call col3(work2,u_grad,w_grad,ntot*ldim)
         call add4(work1,work2(1,1,1,1,1),work2(1,1,1,1,2),
     $     work2(1,1,1,1,3),ntot)
          
         call myavg1    (e13,work1,alpha,beta,ntot ,'e13  ',ifverbose)

        !<e23>t: (dv/dx)*(dw/dx) + (dv/dy)*(dv/dy) + (dv/dz)*(dw/dz)
         call col3(work2,v_grad,w_grad,ntot*ldim)
         call add4(work1,work2(1,1,1,1,1),work2(1,1,1,1,2),
     $     work2(1,1,1,1,3),ntot)
          
         call myavg1    (e23,work1,alpha,beta,ntot ,'e23  ',ifverbose)  

      endif
        
c-----------------------------------------------------------------------
      if ( (mod(istep,iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then
        

         time_temp = time
         time      = atime   ! Output the duration of this avg
         ifpo=.false. 
         ifto=.true.
         call outpost(uavg,vavg,wavg,pr,pravg,     'f01')
         call outpost(urms,vrms,wrms,pr,prrms,     'f02')
         call outpost (uvms,vwms,wums,pr,pums,      'f03')
         call outpost (pvms,pwms,pdudx(1,1,1,1,1),
     $     prms,pdudx(1,1,1,1,2),      'f04')
         call outpost (pdudx(1,1,1,1,3),pdvdx(1,1,1,1,1),
     $     pdvdx(1,1,1,1,2),pr,pdvdx(1,1,1,1,3),     'f05')
         call outpost (pdwdx(1,1,1,1,1),pdwdx(1,1,1,1,2),
     $     pdwdx(1,1,1,1,3),pr,u2ums,     'f06')
         call outpost (v2vms,w2wms,u2vms,pr,u2wms,     'f07')
         
         call outpost (v2ums,v2wms,w2ums,pr,w2vms,    'f08') 
         call outpost (p2pms,uvwms,u4m,pr,v4m,    'f09')
         call outpost (w4m,p4m,e11,pr,e22,    'f10')
         call outpost (e33,e12,e13,pr,e23,    'f11')
         ifto=.false.
         ifpo=.true.
         atime = 0.
         time  = time_temp  ! Restore clock

      endif
c
      timel = time
c
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine avg_all_particle
c
c     This routine computes running averages E(X),E(X^2),E(X*Y)
c     and outputs to avg*.fld*, rms*.fld*, and rm2*.fld* for all
c     fields.
c
c     E denotes the expected value operator and X,Y two
c     real valued random variables.
c
c     variances and covariances can be computed in a post-processing step:
c
c        var(X)   := E(X^X) - E(X)*E(X) 
c        cov(X,Y) := E(X*Y) - E(X)*E(Y)  
c
c     Note: The E-operator is linear, in the sense that the expected
c           value is given by E(X) = 1/N * sum[ E(X)_i ], where E(X)_i
c           is the expected value of the sub-ensemble i (i=1...N).
c
      include 'SIZE'  
      include 'TOTAL' 
      include 'PARTICLES'
!       real concen(lx1,ly1,lz1,lelt,pgrp)
!       real vpdavg(lx1,ly1,lz1,lelt,pgrp,ldim)
!       real updavg(lx1,ly1,lz1,lelt,pgrp,ldim)
!       common /concon/ concen,updavg,vpdavg
      logical ifverbose
      integer icalld2,ntot
      save    icalld2
      data    icalld2  /0/
      character*3 prefix (6)

      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelv
      nto2  = nx2*ny2*nz2*nelv
!         write(6,*)'write made it here!!!!'
      ! initialization
      if (icalld2.eq.0) then
         icalld2 = icalld2 + 1
         aptime  = 0.
         ptimel  = time

         call rzero(upavg,ntot*pgrp)
         call rzero(vpavg,ntot*pgrp)
         call rzero(wpavg,ntot*pgrp)
         
        

         call rzero(uprms,ntot*pgrp)
         call rzero(vprms,ntot*pgrp)
         call rzero(wprms,ntot*pgrp)
         
         
         ! Two point correlations
         call rzero(vpwms,ntot*pgrp)
         call rzero(wpums,ntot*pgrp)
         call rzero(upvms,ntot*pgrp)
         
         call rzero(ufavg,ntot*pgrp)
         call rzero(vfavg,ntot*pgrp)
         call rzero(wfavg,ntot*pgrp)
         
        

         call rzero(ufrms,ntot*pgrp)
         call rzero(vfrms,ntot*pgrp)
         call rzero(wfrms,ntot*pgrp)
         
         
         ! Two point correlations
         call rzero(vfwms,ntot*pgrp)
         call rzero(wfums,ntot*pgrp)
         call rzero(ufvms,ntot*pgrp)
!          Concentration
         call rzero(allcon,ntot*pgrp)
                  
   
      endif

      dptime = time  - ptimel
      aptime = aptime + dptime

      ! dump freq
      iastep = param(68)
      if  (iastep.eq.0) iastep=param(15)   ! same as iostep
      if  (iastep.eq.0) iastep=500

      ifverbose=.false.
      if (icalld2.lt.5) ifverbose=.true.
      if  (mod(istep,iastep).eq.0) ifverbose=.true.

      if (aptime.ne.0..and.dptime.ne.0.) then
         if(nio.eq.0) write(6,*) 'Computing particle statistics ...'
         beta  = dptime/aptime
         alpha = 1.-beta
         ! compute averages E(X)
         do i=1,pgrp
         if(ifverbose) then  
         if(nid.eq.0) write(6,*)'Group: ',i 
         endif
         call myavg1     (allcon(1+(i-1)*ntot),concen(1,1,1,1,i),
     $                   alpha,beta,ntot ,'cncn ',ifverbose) 
!          write(6,*) allcon(1+(i-1)*ntot),' avg concen'
!      $    ,concen(1,1,1,1,i)
         call myavg1    (upavg(1+(i-1)*ntot),vpdavg(1,1,1,1,i,1),
     $                  alpha,beta,ntot ,'upavg',ifverbose)
         call myavg1    (vpavg(1+(i-1)*ntot),vpdavg(1,1,1,1,i,2),
     $                  alpha,beta,ntot ,'vpavg',ifverbose)
         call myavg1    (wpavg(1+(i-1)*ntot),vpdavg(1,1,1,1,i,3),
     $                   alpha,beta,ntot ,'wpavg',ifverbose)


         call myavg1    (ufavg(1+(i-1)*ntot),updavg(1,1,1,1,i,1),
     $                  alpha,beta,ntot ,'ufavg',ifverbose)
         call myavg1    (vfavg(1+(i-1)*ntot),updavg(1,1,1,1,i,2),
     $                  alpha,beta,ntot ,'vfavg',ifverbose)
         call myavg1    (wfavg(1+(i-1)*ntot),updavg(1,1,1,1,i,3),
     $                   alpha,beta,ntot ,'wfavg',ifverbose)
!          all myavg1    (pavg,pr,alpha,beta,nto2 ,'prmp ',ifverbose)
         
         call myavg2    (uprms(1+(i-1)*ntot),vpdavg(1,1,1,1,i,1),
     $                  alpha,beta,ntot ,'uprms',ifverbose)
         call myavg2    (vprms(1+(i-1)*ntot),vpdavg(1,1,1,1,i,2),
     $                  alpha,beta,ntot ,'vprms',ifverbose)
         call myavg2    (wprms(1+(i-1)*ntot),vpdavg(1,1,1,1,i,3),
     $                   alpha,beta,ntot ,'wprms',ifverbose)


         call myavg2    (ufrms(1+(i-1)*ntot),updavg(1,1,1,1,i,1),
     $                  alpha,beta,ntot ,'ufrms',ifverbose)
         call myavg2    (vfrms(1+(i-1)*ntot),updavg(1,1,1,1,i,2),
     $                  alpha,beta,ntot ,'vfrms',ifverbose)
         call myavg2    (wfrms(1+(i-1)*ntot),updavg(1,1,1,1,i,3),
     $                   alpha,beta,ntot ,'wfrms',ifverbose)
         ! compute averages E(X^2) 
         
         call myavg3    (vpwms(1+(i-1)*ntot),vpdavg(1,1,1,1,i,2),
     $                  vpdavg(1,1,1,1,i,3),
     $                  alpha,beta,ntot ,'vpwms',ifverbose)
         call myavg3    (wpums(1+(i-1)*ntot),vpdavg(1,1,1,1,i,3),
     $                  vpdavg(1,1,1,1,i,2),       
     $                  alpha,beta,ntot ,'wpums',ifverbose)
         call myavg3    (upvms(1+(i-1)*ntot),vpdavg(1,1,1,1,i,1),
     $                  vpdavg(1,1,1,1,i,2),       
     $                  alpha,beta,ntot ,'upvms',ifverbose)


         call myavg3    (vfwms(1+(i-1)*ntot),updavg(1,1,1,1,i,2),
     $                  updavg(1,1,1,1,i,3),
     $                  alpha,beta,ntot ,'vfwms',ifverbose)
         call myavg3    (wfums(1+(i-1)*ntot),updavg(1,1,1,1,i,3),
     $                  updavg(1,1,1,1,i,1),       
     $                  alpha,beta,ntot ,'wfums',ifverbose)
         call myavg3    (ufvms(1+(i-1)*ntot),updavg(1,1,1,1,i,1),
     $                  updavg(1,1,1,1,i,2),       
     $                  alpha,beta,ntot ,'ufvms',ifverbose)

        
        enddo
         
      endif
c
c-----------------------------------------------------------------------
      if ((mod(istep,iastep).eq.0.and.istep.gt.1) .or.lastep.eq.1) then

         time_temp = time
         time      = aptime   ! Output the duration of this avg
        do i=1,pgrp
        prefix(1)='p1' ! particle averages and concentration
        prefix(2)='p2' ! fluid averages
        prefix(3)='p3' ! particle rms
        prefix(4)='p4' ! fluid rms
        prefix(5)='p5' ! particle cross-correlations
        prefix(6)='p6' ! fluid cross-correlations
        
         do j=1,6
          write(prefix(j),"(A2,I1)") prefix(j),i
!           write(6,*)prefix(j)
        enddo
        ifto=.true.
        ifpo=.false.
        call outpost(upavg(1+(i-1)*ntot),vpavg(1+(i-1)*ntot)
     $    ,wpavg(1+(i-1)*ntot),pr,allcon(1+(i-1)*ntot), prefix(1))
        ifto=.false.
        call outpost(ufavg(1+(i-1)*ntot),vfavg(1+(i-1)*ntot)
     $    ,wfavg(1+(i-1)*ntot),pr,t, prefix(2))

        call outpost(uprms(1+(i-1)*ntot),vprms(1+(i-1)*ntot)
     $    ,wprms(1+(i-1)*ntot),pr,t, prefix(3))

        call outpost(ufrms(1+(i-1)*ntot),vfrms(1+(i-1)*ntot)
     $    ,wfrms(1+(i-1)*ntot),pr,t, prefix(4))

        call outpost(vpwms(1+(i-1)*ntot),wpums(1+(i-1)*ntot)
     $    ,upvms(1+(i-1)*ntot),pr,t, prefix(5))

        call outpost(vfwms(1+(i-1)*ntot),wfums(1+(i-1)*ntot)
     $    ,ufvms(1+(i-1)*ntot),pr,t, prefix(6))
         
         ifpo=.true.
        enddo 
         aptime = 0.
         time  = time_temp  ! Restore clock

      endif
c
      ptimel = time
      icalld2=icalld2+1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine myavg1(avg,f,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n)
      character*5 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nio.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av1mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine myavg2(avg,f,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n)
      character*5 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*f(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nio.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av2mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine myavg3(avg,f,g,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n),g(n)
      character*5 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*g(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nio.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end

c-----------------------------------------------------------------------
      subroutine myavg4(avg,f,g,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n),g(n)
      character*5 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*f(k)*g(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nio.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av4mnx')
      endif
c
      return
      end
c-------------------------------------------------------------------