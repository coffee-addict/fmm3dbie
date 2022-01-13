program lap_dir_ext_example_tori
  implicit none
  integer i,j,k
  real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
  real *8, allocatable :: srcvals_all(:,:),srcpos(:,:),srccoefs_all(:,:)
  integer igeomtype
  integer npatches,norder,npols,npts
  integer nobj,npsrc,npatches_all,npts_all
  integer ipars(2)
  real *8, allocatable :: targs(:,:),uvs_src(:,:),uvs_targ(:,:)
  integer, allocatable :: ipatch_id_targ(:)
  integer, allocatable :: norders(:),ixyzs(:),iptype(:)
  integer, allocatable :: norders_all(:),ixyzs_all(:),iptype_all(:)
  real *8, allocatable :: xyz_psrc(:,:),radii(:,:),scales(:,:),trans(:,:),rotate(:,:)
  real *8, allocatable :: sigma(:),rhs(:)
  real *8 rhs_tmp
  integer numit,niter
  real *8, allocatable :: errs_gmres(:)
  real *8 eps_fmm,eps_gmres,rres_gmres,l2errnorm_pot
  real *8 dpars(2)
  complex *16 zpars
  integer ifinout
  integer ntarg,ndtarg
  real *8, allocatable :: potex(:),pot(:),charge(:)
  real *8 rel_err,ier
  real *8 pi
  integer unit_geo/999/

  call prini(6,13)
  pi = atan(1.0d0)*4
! solution = \alpha S[\sigma] + \beta D[\sigma]
  dpars = 1.0d0 ! \alpha = \beta = 1
  ifinout = 1 ! exterior problem
  nobj = 4 ! # of objects
  npsrc = nobj ! # of point sources
  igeomtype = 2 ! geometry type (1: sphere, 2: torus)
  if(igeomtype.eq.1) then
    ipars(1) = 1
    npatches = 12*(4**ipars(1))
  else if(igeomtype.eq.2) then
    ipars(1) = 4
    ipars(2) = ipars(1)*3
    npatches = 2*ipars(1)*ipars(2)
  endif
  npatches_all = npatches*nobj

  ndtarg = 3 ! space dimensions of target points
!  ntarg = 1000 ! # of target popints
  ntarg = 1 ! # of target popints
  norder = 6 ! heighest degree of Koornwinder expansion
  npols = (norder+1)*(norder+2)/2 ! # of Koornwinder polynomials of degree < norder
  npts = npatches*npols ! # of discretization points for each patch
  npts_all = npts*nobj ! total # of discretization points for all patches
  write(*,'(a,i0)'), '# of dicretization points = ', npts_all

  allocate(targs(3,ntarg))
  allocate(uvs_src(2,npts))
  allocate(iptype_all(npatches*nobj))
  allocate(radii(3,nobj),scales(3,nobj))
  allocate(trans(3,nobj),rotate(3,nobj))
  allocate(xyz_psrc(3,npsrc))
  call setup_targets(ntarg,targs)
  call setup_radii(nobj,radii)
  call setup_translations(nobj,trans) ! translations of objects from the origin
  call setup_rotation(nobj,rotate)    ! rotations of objects about x,y,z-axes
  call setup_point_sources(nobj,radii,trans,rotate,npsrc,xyz_psrc)

  ! set up common variables for GMRES&FMM solver
  allocate(srcvals(12,npts),srcvals_all(12,npts_all))
  allocate(srccoefs(9,npts),srccoefs_all(9,npts_all))
  allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
  allocate(norders_all(npatches*nobj),ixyzs_all(npatches*nobj+1))
  norders = norder
  iptype = 1
  iptype_all = 1
  do i=1,npatches
    ixyzs(i) = 1 + (i-1)*npols
  enddo
  ixyzs(npatches+1) = 1 + npols*npatches
  uvs_src = 0
  do i=1,nobj
    do j=1,npatches
      norders_all(j+npatches*(i-1)) = norders(j) 
      ixyzs_all(j+npatches*(i-1)) = ixyzs(j) + npols*npatches*(i-1)
    enddo
  enddo
  ixyzs_all(npatches*nobj+1) = npols*npatches*nobj + 1

  ! set up geometric information & right hand side 
  ! passed to GMRES&FMM solver
!  eps_fmm = 1.0d-6
  eps_fmm = 1.0d-10
  allocate(srcpos(3,npts),charge(nobj),rhs(npts_all))
  charge= 1.0d0
  do i=1,nobj
    call setup_geom(nobj,i,igeomtype,npts,norder,npols,npatches,ipars, &
      srcvals,srccoefs,radii,scales,trans,rotate)
    do j=1,npts
      do k=1,3
        srcpos(k,j) = srcvals(k,j)
      enddo
      do k=1,12
        srcvals_all(k,j+npts*(i-1)) = srcvals(k,j)
      enddo
      do k=1,9
        srccoefs_all(k,j+npts*(i-1)) = srccoefs(k,j)
      enddo
    enddo
    ! construct the right hand side of the linear system
    call lfmm3d_t_c_p(eps_fmm,npsrc,xyz_psrc,charge,npts,srcpos,rhs(npts*(i-1)+1),ier)
  enddo
  rhs = rhs/(pi*4.0d0) ! adjust constant factor of lfmm3d_t_c_p

  call output_geometry(npts*nobj,srcvals_all)

  ! solve the linear system with GMRES&FMM solver
  numit = 200
  eps_gmres = eps_fmm
  allocate(errs_gmres(numit+1),sigma(npts_all))
  call lap_comb_dir_solver(npatches_all,norders_all,ixyzs_all,iptype_all, &
    npts_all,srccoefs_all,srcvals_all,eps_fmm,dpars,numit,ifinout,rhs,eps_gmres, &
    niter,errs_gmres,rres_gmres,sigma)

  ! evaluate layer potentials for the solutions using determined density values
  allocate(ipatch_id_targ(ntarg),uvs_targ(2,ntarg))
  ipatch_id_targ = -1
  uvs_targ = 0.0d0
  allocate(pot(ntarg))
  call lpcomp_lap_comb_dir(npatches_all,norders_all,ixyzs_all,iptype_all, &
    npts_all,srccoefs_all,srcvals_all,ndtarg,ntarg,targs,ipatch_id_targ, &
    uvs_targ,eps_fmm,dpars,sigma,pot)

  ! exact solutions at target points
  allocate(potex(ntarg))
  call lfmm3d_t_c_p(eps_fmm,npsrc,xyz_psrc,charge,ntarg,targs,potex,ier)
  potex = potex/(pi*4.0d0)

  write(*,'(a)'), ''
  write(*,'(a,i0)'), '# of target points = ', ntarg
  write(*,'(a,i0)'), '# of tori = ', nobj
  l2errnorm_pot = sqrt(sum(((pot-potex)/potex)**2))
  write(*,'(a,a,e15.6)'), 'L2 norm of relative errors in potentials ', &
    'at target points = ', l2errnorm_pot

!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,npsrc 
    write(*,'(a,i0,a,e15.6,e15.6,e15.6)'), 'point source of object (', i, &
      ') = ', xyz_psrc(1,i), xyz_psrc(2,i), xyz_psrc(3,i)
  enddo
  do i=1,ntarg
    write(*,'(a,i0,a,e15.6)'), 'pot(', i, ') = ', pot(i)
    write(*,'(a,i0,a,e15.6)'), 'potex(', i, ') = ', potex(i)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!

  stop
end

subroutine setup_geom(nobj,iobj,igeomtype,npts,norder,npols,npatches, &
  ipars,srcvals,srccoefs,radii,scales,trans,rotate)
  implicit none
  integer igeomtype,nobj,iobj,npts,npols,norder,npatches,ipars(*)
  real *8 srcvals(12,*), srccoefs(9,*)
  real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4
  integer, pointer :: iptr1,iptr2,iptr3,iptr4
  real *8, target :: p1(10),p2(10),p3(10),p4(10)
  real *8, allocatable, target :: triaskel(:,:,:)
  integer, allocatable :: isides(:)
  real *8 radii(3,nobj),scales(3,nobj),trans(3,nobj),rotate(3,nobj)
  real *8 pi
  procedure (), pointer :: xtri_geometry
  external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval

  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
  allocate(wts(npols))
  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
  if(igeomtype.eq.1) then
    call setup_sphere(npts,norder,npols,npatches,ipars,srcvals,&
      srccoefs,trans(1,iobj))
  else if(igeomtype.eq.2) then
    call setup_torus(nobj,iobj,npts,norder,npols,npatches,ipars,srcvals,&
      srccoefs,trans(1,iobj),rotate(1,iobj),radii,scales)
  endif
  return  
end subroutine

subroutine setup_sphere(npts,norder,npols,npatches,ipars,srcvals, &
  srccoefs,trans)
  implicit none
  integer npts,norder,npols,npatches,igeomtype,ntri,ipars(*)
  real *8 srcvals(12,*), srccoefs(9,*)
  real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

  real *8, pointer :: ptr1,ptr2,ptr3,ptr4
  real *8, target :: p1(10),p2(10),p3(10),p4(10)
  real *8, allocatable, target :: triaskel(:,:,:)
  integer, allocatable :: isides(:)
  real *8 trans(3)
  procedure (), pointer :: xtri_geometry

  external xtri_sphere_eval
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
  allocate(wts(npols))

  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

  igeomtype = 2
  allocate(triaskel(3,3,npatches))
  allocate(isides(npatches))
  ntri = npatches
  call xtri_platonic(igeomtype,ipars(1),ntri,npatches,triaskel,isides)

  xtri_geometry => xtri_sphere_eval
  ptr1 => triaskel(1,1,1)
  ptr2 => p2(1)
  ptr3 => p3(1)
  ptr4 => p4(1)

  call getgeominfo_translate(npatches,xtri_geometry,ptr1,ptr2, &
    ptr3,ptr4,npols,uvs,umatr,trans,srcvals,srccoefs)

  return  
end subroutine

subroutine setup_torus(nobj,iobj,npts,norder,npols,npatches,ipars,srcvals, &
  srccoefs,trans,rotate,radii,scales)
  implicit none
  integer nobj,iobj,npts,norder,npols,npatches,nover,ipars(*)
  real *8 srcvals(12,*), srccoefs(9,*)
  real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)
  real *8, pointer :: ptr1,ptr2,ptr3,ptr4
  integer, pointer :: iptr1,iptr2,iptr3,iptr4
  real *8, allocatable, target :: triaskel(:,:,:)
  integer, target :: oscillaiton
  real *8 :: radii(3,nobj), scales(3,nobj)
  real *8, target :: radii0(3), scales0(3)
  real *8 trans(3),rotate(3)
  real *8 umin,umax,vmin,vmax
  real *8 pi
  procedure (), pointer :: xtri_geometry

  external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval
  allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
  allocate(wts(npols))

  call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

  radii0(1) = radii(1,iobj) ! minor
  radii0(2) = radii(2,iobj) ! major
  radii0(3) = radii(3,iobj) ! oscillation?
  scales0(1) = 1.0d0 ! x
  scales0(2) = 1.0d0 ! y
  scales0(3) = 1.0d0 ! z
  pi = atan(1.0d0)*4
  umin = 0
  umax = 2*pi
  vmin = 0
  vmax = 2*pi
  allocate(triaskel(3,3,npatches))
  nover = 0
  call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),&
    nover,npatches,npatches,triaskel)

  xtri_geometry => xtri_wtorus_eval
  oscillaiton = 0

  ptr1 => triaskel(1,1,1)
  ptr2 => radii0(1)
  ptr3 => scales0(1)
  iptr4 => oscillaiton

  call getgeominfo_translate_rotate(npatches,xtri_geometry,ptr1,ptr2,&
    ptr3,iptr4,npols,uvs,umatr,trans,rotate,srcvals,srccoefs)

  return  
end subroutine

subroutine setup_targets(ntarg,targs)
  integer ntarg
  real *8 targs(3,ntarg)

  targs(1,1) = 0.0d0
  targs(2,1) = 0.0d0
  targs(3,1) = 0.0d0

!  targs(1,2) = 0.0d0 
!  targs(2,2) = 1.0d0
!  targs(3,2) = 0.0d0

!!  do i=3,ntarg
!  do i=1,ntarg
!    targs(1,i) = -15.0d0 + 30.0d0*real(i,8)/real(ntarg,8)
!    targs(2,i) = 15.0d0
!    targs(3,i) = 15.0d0
!  enddo
end subroutine

subroutine setup_radii(nobj,radii)
  ! set up radii of tori
  integer nobj
  real *8 radii(3,nobj)

  radii(1,1) = 1.0d0
  radii(2,1) = 3.0d0
  radii(3,1) = 0.0d0

  radii(1,2) = 1.0d0 
  radii(2,2) = 4.0d0
  radii(3,2) = 0.0d0

  radii(1,3) = 1.0d0 
  radii(2,3) = 5.0d0
  radii(3,3) = 0.0d0

  radii(1,4) = 1.0d0 
  radii(2,4) = 4.0d0
  radii(3,4) = 0.0d0
end subroutine

subroutine setup_translations(nobj,trans)
  ! set translation parameters for objects from the origin
  implicit real *8 (a-h,o-z)
  integer nobj
  real *8 trans(3,nobj)

  trans(1,1) = 7.0d0
  trans(2,1) = 0.0d0
  trans(3,1) = 0.0d0

  trans(1,2) = -10.0d0
  trans(2,2) = 0.0d0
  trans(3,2) = -2.0d0

  trans(1,3) = 0.0d0
  trans(2,3) = -10.0d0
  trans(3,3) = -5.0d0

  trans(1,4) = -10.0d0
  trans(2,4) = 6.0d0
  trans(3,4) = 8.0d0
end subroutine

subroutine setup_rotation(nobj,rotate)
  ! set rotation parameters for objects about x,y,z-axes
  implicit real *8 (a-h,o-z)
  integer nobj
  real *8 rotate(3,nobj),pi

  pi = atan(1.0d0)*4

  rotate(1,1) = 0.0d0
  rotate(2,1) = pi*0.25d0
  rotate(3,1) = 0.0d0

  rotate(1,2) = 0.0d0
  rotate(2,2) = 0.0d0
  rotate(3,2) = 0.0d0

  rotate(1,3) = pi*0.25d0
  rotate(2,3) = 0.0d0
  rotate(3,3) = 0.0d0

  rotate(1,4) = 0.0d0
  rotate(2,4) = pi*0.5d0
  rotate(3,4) = 0.0d0
end subroutine

subroutine setup_point_sources(nobj,radii,trans,rotate,npsrc,xyz_psrc)
  implicit none
  integer i
  integer nobj,npsrc
  real *8 xyz_psrc(3,npsrc),radii(3,nobj),trans(3,nobj),rotate(3,nobj),tmp(3)

  do i=1,npsrc
    ! put source at (r,0,0) inside the object
    xyz_psrc(1,i) = radii(2,i)

    ! rotations
    ! about x-axis
    tmp(2) = xyz_psrc(2,i)*cos(rotate(1,i)) - xyz_psrc(3,i)*sin(rotate(1,i));
    tmp(3) = xyz_psrc(2,i)*sin(rotate(1,i)) + xyz_psrc(3,i)*cos(rotate(1,i));
    xyz_psrc(2,i) = tmp(2);
    xyz_psrc(3,i) = tmp(3);
    ! about y-axis
    tmp(1) = xyz_psrc(1,i)*cos(rotate(2,i)) + xyz_psrc(3,i)*sin(rotate(2,i));
    tmp(3) = -xyz_psrc(1,i)*sin(rotate(2,i)) + xyz_psrc(3,i)*cos(rotate(2,i));
    xyz_psrc(1,i) = tmp(1);
    xyz_psrc(3,i) = tmp(3);
    ! about z-axis
    tmp(1) = xyz_psrc(1,i)*cos(rotate(3,i)) - xyz_psrc(2,i)*sin(rotate(3,i));
    tmp(2) = xyz_psrc(1,i)*sin(rotate(3,i)) + xyz_psrc(2,i)*cos(rotate(3,i));
    xyz_psrc(1,i) = tmp(1);
    xyz_psrc(2,i) = tmp(2);

    ! transrations
    xyz_psrc(1,i) = xyz_psrc(1,i) + trans(1,i)
    xyz_psrc(2,i) = xyz_psrc(2,i) + trans(2,i)
    xyz_psrc(3,i) = xyz_psrc(3,i) + trans(3,i)
    write(*,'(a,i0,a,e15.6,e15.6,e15.6)'), 'point source of object (', i, &
      ') = ', xyz_psrc(1,i), xyz_psrc(2,i), xyz_psrc(3,i)
  enddo
end subroutine

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

subroutine getgeominfo_translate(npatches, patchpnt, par1, par2, &
    par3, par4, npols, uvs, umatr, trans, srcvals, srccoefs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,npols), srcvals(12,*)
  real *8 :: umatr(npols,npols),srccoefs(9,*)
  external patchpnt

  real *8 :: xyz(3),dxyzduv(3,10),xyznorm(3),trans(3)
  real *8 :: xyztang1(3),xyztang2(3)

  !
  !       This subroutine return all points, normals and tangents from
  !       geometry descriptor
  !
  !       Input parameters:
  !
  !         npatches: integer: the number of patches
  !         patchpnt: external: subroutine evaluating points along
  !               the surface, given patch by patch, must be of the form
  !                     patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
  !         par1,par2,par3,par4: extra parameters
  !         npols: integer: the total number of polynomials for each patch
  !         uvs: real *8(2,npols): local u-discretization points for each patch
  !         umatr: real *8(npols,npols): values to coeffs matrix on standard patch 
  !         trans: real *8(3): translation from the origin in R^3 space
  !
  !       Output parameters:
  !
  !         srcvals: real *8(12,npts): geometry info with first derivatives
  !               srcvals(1:3,:) - xyz
  !               srcvals(4:6,:) - dxyz/du
  !               srcvals(7:9,:) - dxyz/dv
  !               srcvals(10:12,:) - xyznorms
  !
  !         srccoefs: real *8 (9,npts): geometry info as koornwinder expansion
  !                     coefficients
  !                    
  !         npts: integer: the total number of points in discretization
  !

  do ipatch=1,npatches
    do i=1,npols

      u=uvs(1,i)
      v=uvs(2,i)

      ipt = (ipatch-1)*npols + i 

      call patchpnt(ipatch,u,v,srcvals(1,ipt),srcvals(4,ipt),par1, &
             par2,par3,par4)

      srcvals(1,ipt) = srcvals(1,ipt) + trans(1)
      srcvals(2,ipt) = srcvals(2,ipt) + trans(2)
      srcvals(3,ipt) = srcvals(3,ipt) + trans(3)

      call cross_prod3d(srcvals(4,ipt),srcvals(7,ipt),srcvals(10,ipt))

      ds = sqrt(srcvals(10,ipt)**2 + srcvals(11,ipt)**2 + &
              srcvals(12,ipt)**2)
      srcvals(10,ipt) = srcvals(10,ipt)/ds
      srcvals(11,ipt) = srcvals(11,ipt)/ds
      srcvals(12,ipt) = srcvals(12,ipt)/ds

    end do

    do i=1,npols
      ipt = (ipatch-1)*npols + i
      do j=1,9
        srccoefs(j,ipt) = 0
        do l=1,npols
          lpt = (ipatch-1)*npols + l
          srccoefs(j,ipt) = srccoefs(j,ipt) + umatr(i,l)*srcvals(j,lpt)
        end do
      end do
    end do
  end do

  npts = npatches*npols

  return
end subroutine getgeominfo_translate

subroutine getgeominfo_translate_rotate(npatches, patchpnt, par1, par2, &
    par3, par4, npols, uvs, umatr, trans, theta, srcvals, srccoefs)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,npols), srcvals(12,*)
  real *8 :: umatr(npols,npols),srccoefs(9,*)
  external patchpnt

  real *8 :: xyz(3),dxyzduv(3,10),xyznorm(3),trans(3),theta(3),tmp(3)
  real *8 :: xyztang1(3),xyztang2(3)

  !
  !       This subroutine return all points, normals and tangents from
  !       geometry descriptor
  !
  !       Input parameters:
  !
  !         npatches: integer: the number of patches
  !         patchpnt: external: subroutine evaluating points along
  !               the surface, given patch by patch, must be of the form
  !                     patchpnt(ipatch,u,v,xyz,dxyzduv,par1,par2,par3,par4)
  !         par1,par2,par3,par4: extra parameters
  !         npols: integer: the total number of polynomials for each patch
  !         uvs: real *8(2,npols): local u-discretization points for each patch
  !         umatr: real *8(npols,npols): values to coeffs matrix on standard patch 
  !         trans: real *8(3): translation from the origin in R^3 space
  !         theta: real *8(3): counter crockwise rotation angles in radian 
  !                            about x,y,z-axes 
  !
  !       Output parameters:
  !
  !         srcvals: real *8(12,npts): geometry info with first derivatives
  !               srcvals(1:3,:) - xyz
  !               srcvals(4:6,:) - dxyz/du
  !               srcvals(7:9,:) - dxyz/dv
  !               srcvals(10:12,:) - xyznorms
  !
  !         srccoefs: real *8 (9,npts): geometry info as koornwinder expansion
  !                     coefficients
  !                    
  !         npts: integer: the total number of points in discretization
  !

  do ipatch=1,npatches
    do i=1,npols

      u=uvs(1,i)
      v=uvs(2,i)

      ipt = (ipatch-1)*npols + i 

      call patchpnt(ipatch,u,v,srcvals(1,ipt),srcvals(4,ipt),par1, &
             par2,par3,par4)

      ! rotations
      ! about x-axis
      tmp(2) = srcvals(2,ipt)*cos(theta(1)) - srcvals(3,ipt)*sin(theta(1));
      tmp(3) = srcvals(2,ipt)*sin(theta(1)) + srcvals(3,ipt)*cos(theta(1));
      srcvals(2,ipt) = tmp(2);
      srcvals(3,ipt) = tmp(3);
      ! about y-axis
      tmp(1) = srcvals(1,ipt)*cos(theta(2)) + srcvals(3,ipt)*sin(theta(2));
      tmp(3) = -srcvals(1,ipt)*sin(theta(2)) + srcvals(3,ipt)*cos(theta(2));
      srcvals(1,ipt) = tmp(1);
      srcvals(3,ipt) = tmp(3);
      ! about z-axis
      tmp(1) = srcvals(1,ipt)*cos(theta(3)) - srcvals(2,ipt)*sin(theta(3));
      tmp(2) = srcvals(1,ipt)*sin(theta(3)) + srcvals(2,ipt)*cos(theta(3));
      srcvals(1,ipt) = tmp(1);
      srcvals(2,ipt) = tmp(2);

      ! translation
      srcvals(1,ipt) = srcvals(1,ipt) + trans(1)
      srcvals(2,ipt) = srcvals(2,ipt) + trans(2)
      srcvals(3,ipt) = srcvals(3,ipt) + trans(3)

      call cross_prod3d(srcvals(4,ipt),srcvals(7,ipt),srcvals(10,ipt))

      ds = sqrt(srcvals(10,ipt)**2 + srcvals(11,ipt)**2 + &
              srcvals(12,ipt)**2)
      srcvals(10,ipt) = srcvals(10,ipt)/ds
      srcvals(11,ipt) = srcvals(11,ipt)/ds
      srcvals(12,ipt) = srcvals(12,ipt)/ds

    end do

    do i=1,npols
      ipt = (ipatch-1)*npols + i
      do j=1,9
        srccoefs(j,ipt) = 0
        do l=1,npols
          lpt = (ipatch-1)*npols + l
          srccoefs(j,ipt) = srccoefs(j,ipt) + umatr(i,l)*srcvals(j,lpt)
        end do
      end do
    end do
  end do

  npts = npatches*npols

  return
end subroutine getgeominfo_translate_rotate
