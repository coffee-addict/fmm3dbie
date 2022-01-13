      implicit none
      integer i,j,k
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: srcvals_all(:,:),srcpos(:,:)
      real *8, allocatable :: srccoefs_all(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer npatches, npts, npts_all, nobj
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)

      complex *16 vf(3)

      real *8, allocatable :: errs(:)
      real *8 rres
      real *8 thet,phi
	  complex * 16  zpars(5), omega, ep0,mu0,ep1,mu1

      integer numit,niter
      character *200 title,fname,fname1,fname2

      integer ipatch_id
      real *8 uvs_targ(2)

      integer ifinout
      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer count1
      real *8 done
      real *8 eps_quad_fmm, eps_gmres

      real *8 t1, t2

      real *8 pi
      real *8 hkrand

      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

!
!   simulation for .go3 files of DFIE
!

!
!	Electromagnetic parameters
!
  
      omega=1.1d0
      ep0=1.0d0
      mu0=1.0d0
      ep1=1.200d0
      mu1=1.00d0

      zpars(1) = omega
      zpars(2) = ep0
      zpars(3) = mu0
      zpars(4) = ep1
      zpars(5) = mu1
  
      write (*,*) 'omega: ',omega
      write (*,*) 'ep1: ',ep1
      
!      fname = '../../geometries/sphere_192_o03.go3'
      fname = '../../geometries/y_coupler_o04_r00.go3'

      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.37d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 5.1d0

      nobj = 1
      call open_gov3_geometry_mem(fname,npatches,npts)
      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,&
     &iptype,npts,srcvals,srccoefs,wts)

      allocate(sigma(4*npts),rhs(4*npts))
      ifinout = 1

      call srand(0)
      hkrand = rand()
      thet = hkrand*pi
      phi = hkrand*2*pi
      write(*,'(a,e15.6)'), 'hkrand = ', hkrand
      
!
!	Set the orientation of the Electric/Magnetic dipoles
!

      vf(1)=1.0d0
      vf(2)=2.0d0
      vf(3)=3.0d0
!	  
      call get_rhs_em_muller_trans(xyz_in,xyz_out,vf,npts,srcvals,zpars, &
        rhs)
  
      numit = 400
      niter = 0
      allocate(errs(numit+1))

!
!	Specify tolerance
!

      eps_quad_fmm = 1d-4
      eps_gmres=1d-6
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      
      call em_muller_trans_solver(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,eps_quad_fmm,zpars,numit,ifinout,rhs, &
        eps_gmres,niter,errs,rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
!C$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)

!
!       test solution at interior and exterior point
!

      call test_accuracy_em_muller_trans(eps_quad_fmm,sigma,zpars,npts,&
        wts,srcvals,xyz_in,vf,xyz_out)
      
      allocate(srcvals_all(12,npts*nobj))
      do j=1,npts
        do k=1,12
          srcvals_all(k,j) = srcvals(k,j)
        enddo
      enddo
      call output_geometry(npts*nobj,srcvals_all)

      stop
      end

      subroutine setup_geom(igeomtype,norder,npatches,ipars,& 
     &srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri,& 
       &triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),&
     &nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &iptr3,iptr4, norder,'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

      subroutine test_exterior_pt(npatches,norder,npts,srcvals,&
     &srccoefs,wts,xyzout,isout)

!
!  this subroutine tests whether the pt xyzin, is
!  in the exterior of a surface, and also estimates the error
!  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
!  centered at the interior point. Whether a point 
!  is in the interior or not is tested using Gauss' 
!  identity for the flux due to a point charge
!
!
!  input:
!    npatches - integer
!       number of patches
!    norder - integer
!       order of discretization
!    npts - integer
!       total number of discretization points on the surface
!    srccoefs - real *8 (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  real *8 (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2

      zk = 0
      ra = 0

      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end

      subroutine output_geometry(npts_all,srcvals_all)
        integer i,npts_all,unit_geo
        real *8 srcvals_all(12,npts_all)
        character *100 fname_geom
        fname_geom = 'geometry.txt'
        open(unit_geo, file=fname_geom, status='REPLACE')
        do i=1,npts_all
          write(unit_geo, '(e15.6,e15.6,e15.6)'), srcvals_all(1,i), &
            srcvals_all(2,i), srcvals_all(3,i)
        enddo
        close(unit_geo)
      end subroutine
   




