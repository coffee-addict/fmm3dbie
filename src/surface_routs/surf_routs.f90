!    This file contains the following user callable subroutines
!      surf_fun_error - estimate the resolution of a function on
!         a surface
!
!      get_centroid_rads - compute the centroid and bounding sphere
!          radii for a collection of patches
!         
!      get_centroid_rads_tri - compute the centroid and bounding
!          sphere radius for a given triangular patch
!
!      oversample_geom - oversample xyz,dxyz/du,dxyz/dv on the surface
!
!      oversample_fun_surf - oversample functions defined on a surface
!
!      oversample_fun_tri - oversample function on a triangle
!
!      surf_vals_to_coefs - given function values defined on a patch
!         convert to its koornwinder expansion coeffs 
!
!      vals_to_coefs_tri - given function values defined on a patch
!         convert to its koornwinder expansion coeffs 
!
!      get_qwts - compute smooth quadrature weights on the surface
!
!      get_near_far_split_pt - split the near field into the near-near
!         field (done via adaptive integration) and the near-far field
!         (done via oversampled quadrature (more oversampling than the
!         far-far field))
!
!      get_patch_id_uvs - for all boundary points, get patch id they are on
!          and u,v location on patch
!
!      get_patch_distortion - estimate shape distortion parameter for a given
!          collection of patches
!
!      get_nximat - estimate number of elements in the union of 
!        interpolation matrices for all patches
!
!      get_ximats - get all interpolation matrices for the surface
!   
!      get_first_fundamental_form - get the first fundamental form at all
!        points
!
!      get_inv_first_fundamental_form - get inverse of first fundamental 
!         form
!
!      get_surf_grad - compute surface gradient of a scalar function
!
!      get_surf_grad_fast - compute surface gradient of a scalar function 
!        (With precomputed inverse of first fundamental form)
!
!      get_surf_uv_grad - compute the uv gradient of a function defined
!         on a surface, i.e. \partial_{u}f, \partial_{v}f
!
!      get_surf_div - compute the surface divergence of a vector field
!         defined on a surface
!
!      col_ind_to_patch_node_ind - convert a list of column
!         indices to a list of patch and node indices
!
!         



!
!
!
!
!

subroutine surf_fun_error(nd,npatches,norders,ixyzs,iptype, &
   npts,dvals,wts,errp,errm)

!
!  This subroutine computes the error in resolving
!  a function of the surface whose patchwise
!  koornwinder expansion coefficients are provided
!
!  nd - number of functions
!  npatches - number of patches
!  norders(npatches) - order of discretization
!  ixyzs(npatches+1) - location where expansion coeffs for patch i start
!  iptype(npatches) - type of patch
!      iptype = 1, triangular patch with RV nodes
!  dvals(nd,npts) - values of the function at discretization
!                   nodes
!  wts(npts) - quadrature weights for integrating smooth functions on
!              the surface
!  
!
!  output:
!  errp(npatches) - error per patch
!  errm - max of errp
!

  implicit real *8 (a-h,o-z)
  integer norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches)
  real *8 dvals(nd,npts),errp(npatches),errm
  real *8 wts(npts)
  real *8, allocatable :: dcoefs(:,:),dtail(:,:),uvs(:,:),pols(:)
  real *8 rl2s,rl2tails

  allocate(dcoefs(nd,npts),dtail(nd,npts))

  do i=1,npts
    do idim=1,nd
      dcoefs(idim,i) = 0
      dtail(idim,i) = 0
    enddo
  enddo

  call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
   dvals,dcoefs)
  
  
  errm = 0
  rl2s = 0
  rl2tails = 0


  nmax = 20
  npmax = (nmax+1)*(nmax+2)

  allocate(uvs(2,npmax),pols(npmax))
  do ip=1,npatches
    norder = norders(ip)
    istart0 = ixyzs(ip)
    npols = ixyzs(ip+1)-ixyzs(ip)

    if(iptype(ip).eq.1) then
      norderhead = norder-1
      i1 = (norderhead+1)*(norderhead+2)/2
      call get_vioreanu_nodes(norder,npols,uvs)
    endif


    if(iptype(ip).eq.1) then
      do i=1,npols
        ipt = ixyzs(ip)+i-1
        call koorn_pols(uvs(1,i),norder,npols,pols)
        do idim=1,nd
          dtail(idim,ipt) = 0
          do j=i1+1,npols
            ic = ixyzs(ip)+j-1
            dtail(idim,ipt) = dtail(idim,ipt) + pols(j)*dcoefs(idim,ic)
          enddo
        enddo
      enddo
    endif

    errp(ip) = 0
    do idim=1,nd
      do i=1,npols
        ipt = ixyzs(ip)+i-1
        errp(ip) = errp(ip) + dtail(idim,ipt)**2*wts(ipt)
        rl2s = rl2s + dvals(idim,ipt)**2*wts(ipt)
        rl2tails = rl2tails + dtail(idim,ipt)**2*wts(ipt)
      enddo
    enddo
    errp(ip) = sqrt(errp(ip))
    if(errp(ip).gt.errm) errm = errp(ip)
  enddo
  
  rl2s = sqrt(rl2s)
  rl2tails = sqrt(rl2tails)

  do ip=1,npatches
    errp(ip) = errp(ip)/errm*rl2tails/rl2s
  enddo

  errm = rl2tails/rl2s


end subroutine surf_fun_error
!
!
!
!
!




subroutine get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
    srccoefs,cms,rads)
!
!   this subroutine computes the centroid of each patch and the radius 
!   of the bounding sphere centered at the centroid
!
!   Input:
!      npatches - integer
!         number of patches
!      norders - integer
!         order of discretization of each patches
!      ixyzs - integer(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!             and srcvals array corresponding to patch i
!   
!      iptype - integer(npatches)
!         type of patch
!          iptype = 1, triangular patch discretized using RV nodes
!
!     npts - integer
!         total number of discretization points
!      srccoefs - real *8 (9,npts)
!         koornwinder expansions of the geometry info
!           x,y,z,,dx/du,dy/du,dz/du,dx/dv,dy/dv,dz/dv
!
!

  implicit none
  integer npatches,norders(npatches),npols,npts,ixyzs(npatches+1)
  integer iptype(npatches)
  real *8 srccoefs(9,npts),cms(3,npatches),rads(npatches)

  integer i,istart


  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    if(iptype(i).eq.1) call get_centroid_rads_tri(norders(i),npols, &
      srccoefs(1,istart),cms(1,i),rads(i))
  enddo
end subroutine get_centroid_rads
!
!
!
!

subroutine get_centroid_rads_tri(norder,npols,srccoefs,cms,rads) 

  integer norder,npols
  real *8 srccoefs(9,npols),cms(3),rads
  real *8 uv(2,3)
  real *8 xyz(3,3),rad
  real *8, allocatable :: pols(:,:)

  integer i,j,l,m,lpt,np
  

  uv(1,1) = 0
  uv(2,1) = 0

  uv(1,2) = 0
  uv(2,2) = 1

  uv(1,3) = 1
  uv(2,3) = 0

  allocate(pols(npols,3))
  do i=1,3
    call koorn_pols(uv(1,i),norder,npols,pols(1,i))
  enddo

  do l=1,3
    cms(l) = 0
  enddo
  rads = 0

!
!     loop over all three vertices
!
  do j=1,3
    do m=1,3
      xyz(m,j) = 0
    enddo
    
    do l=1,npols
      do m=1,3
        xyz(m,j) = xyz(m,j) + srccoefs(m,l)*pols(l,j)
      enddo
    enddo
  enddo

  do m=1,3
    do l=1,3
      cms(m) = cms(m) + xyz(m,l)
    enddo
    cms(m) = cms(m)/3
  enddo
!
!    compute radius of bounding sphere
!
  do m=1,3
    rad = (cms(1)-xyz(1,m))**2 + (cms(2)-xyz(2,m))**2 + &
        (cms(3)-xyz(3,m))**2
    if(rad.ge.rads) rads = rad
  enddo
  rads = sqrt(rads)
end subroutine get_centroid_rads_tri
!
!
!
!
!
subroutine oversample_geom(npatches,norders,ixyzs,iptype, &
    npts,srccoefs,srcvals,nfars,ixyzso,nptso,srcover)
!
!
!  This subroutine oversamples geometry information
!  given expansion coeffs of geometry info on patches
!
!
!  Extreme care must be taken when an oversampled node is identical
!  to a discretization node. The routine ensures that in this event
!  the source information is copied over from the original
!  discretization
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - nfars: integer(npatches)
!        oversampled order of patch i
!    - ixyzso: integer(npatches+1)
!        starting location of oversampled points on patch i
!    - nptso: integer
!        total number of oversampled points
!
!  Output arguments:   
!    
!    - srcover: double precision(12,nptso)
!        Oversampled geometry information
!
!----------------------------------------------
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches),ixyzso(npatches+1)
  integer, intent(in) :: nfars(npatches)
  integer, intent(in) :: npts,nptso
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: srcover(12,nptso)
  integer nfar,norder
  real *8 tmp(3),rr
  

  real *8, allocatable :: pmat(:,:),pols(:)
  real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
  character *1 transa,transb
  real *8 alpha, beta

  integer i,istart,istartover,nfar_pols,j,jpt,npols,l
  integer n1,n2,ipt,i1,i2


  transa = 'n'
  transb = 'n'
  alpha = 1
  beta = 0
  do i=1,npatches
    nfar = nfars(i)
    norder = norders(i)

    if(nfar.ne.norder) then
      if(iptype(i).eq.1) then
        nfar_pols = ixyzso(i+1)-ixyzso(i)
        npols = ixyzs(i+1)-ixyzs(i)
        allocate(uvs(2,nfar_pols))
        call get_vioreanu_nodes(nfar,nfar_pols,uvs)


        allocate(pmat(npols,nfar_pols))
    
        do j=1,nfar_pols
          call koorn_pols(uvs(1,j),norder,npols,pmat(1,j))
        enddo
      endif

      istart = ixyzs(i)
      istartover = ixyzso(i)
      call dgemm_guru(transa,transb,9,nfar_pols,npols,alpha, &
          srccoefs(1,istart),9,pmat,npols,beta,srcover(1,istartover),12)
      do j=1,nfar_pols
        jpt = istartover+j-1
        call cross_prod3d(srcover(4,jpt),srcover(7,jpt),tmp)
        rr = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
        srcover(10,jpt) = tmp(1)/rr
        srcover(11,jpt) = tmp(2)/rr
        srcover(12,jpt) = tmp(3)/rr
      enddo

!
!    Fix clash in coordinates across orders
!
      if(iptype(i).eq.1) then
        n1 = mod(norder,3)
        n2 = mod(nfar,3)

        if(n1.eq.0.and.n2.eq.0) then
          i1 = 1
          i2 = 1
          if(norder.eq.15) i1 = 76
          if(norder.eq.18) i1 = 106

          if(nfar.eq.15) i2 = 76
          if(nfar.eq.18) i2 = 106

          ipt = ixyzs(i)+i1-1
          jpt = ixyzso(i) + i2-1
          do l=1,12
            srcover(l,jpt) = srcvals(l,ipt)
          enddo
        endif
      endif
      deallocate(uvs,pmat)
    else
      npols = ixyzs(i+1)-ixyzs(i)
      istart = ixyzs(i)
      istartover = ixyzso(i)
      do j=1,npols
        do l=1,12
          srcover(l,j+istartover-1) = srcvals(l,j+istart-1)
        enddo
      enddo
    endif
  enddo
end subroutine oversample_geom





subroutine oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,npts, &
   u,nfars,ixyzso,nptso,uover)
!
!  This subroutine oversamples a function defined on a surface
!
!  input
!    nd - integer
!      number of functions
!    npatches - integer
!      number of patches on the surface
!    norders - integer(npatches)
!      order of discretization of patch i
!    ixyzs - integer(npatches+1)
!      starting location of points on patch i
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangle discretized using RV nodes
!    npts - integer
!      total number of points on the surface
!    u - real *8 (nd,npts)
!      functions at the discretization points on the surface
!    nfars - integer(npatches)
!      oversampled order of patch i
!    ixyzso - integer(npatches+1)
!      starting location of oversampled points on patch i
!    nptso - integer
!      total number of oversampled points
!  
!
!  output
!    uover - real *8 (nd,nptso)
!       oversampled function
!    
   

  implicit none
  integer nd,npatches,norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches),npts
  real *8 u(nd,npts)
  integer nfars(npatches),ixyzso(npatches+1),nptso
  real *8 uover(nd,nptso)
  integer i,istart,istarto,npols,npolso

  do i=1,npatches
    
    istart = ixyzs(i)
    istarto = ixyzso(i)
    npols = ixyzs(i+1)-ixyzs(i)
    npolso = ixyzso(i+1) - ixyzso(i)
    if(iptype(i).eq.1) &
      call oversample_fun_tri(nd,norders(i),npols,u(1,istart),&
      nfars(i),npolso,uover(1,istarto))
  enddo

end subroutine oversample_fun_surf



!
!
!
!
subroutine oversample_fun_tri(nd,norder,npols,u,nfar,nfar_pols, &
    uover)
!
!  This subroutine oversample a function on the standard triangle
!
!  input
!    nd - integer
!       number of functions
!    norder - integer
!       order of original discretization
!    npols - integer
!       number of discretization points corresponding to norder
!       npols = (norder+1)*(norder+2)/2
!    u - real *8 (nd,npols)
!       function values tabulated at original grid
!    nfar - integer
!       oversampled order of discretization
!    nfar_pols - integer
!       number of discretization nodes corresponding to oversampled
!         order
!    
!  output
!    uover - real *8 (nd,nfar_pols)
!      oversampled function values
!


  implicit none
  integer norder,npover,nd,npols,nfar
  real *8 u(nd,npols),uover(nd,nfar_pols)
  real *8, allocatable :: ucoefs(:,:)
  real *8, allocatable :: pmat(:,:),pols(:)
  real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
  real *8, allocatable :: umat0(:,:),uvs0(:,:),vmat0(:,:),wts0(:)
  character *1 transa,transb
  real *8 alpha, beta
  integer i,nfar_pols


  allocate(uvs(2,nfar_pols),wts(nfar_pols))
  allocate(umat(nfar_pols,nfar_pols),vmat(nfar_pols,nfar_pols))
  call vioreanu_simplex_quad(nfar,nfar_pols,uvs,umat,vmat,wts)

  allocate(ucoefs(nd,npols))

  allocate(umat0(npols,npols),vmat0(npols,npols),uvs0(2,npols))
  allocate(wts0(npols))


  call vioreanu_simplex_quad(norder,npols,uvs0,umat0,vmat0,wts0)


!
!  compute ucoefs
!

  transa = 'n'
  transb = 't'
  alpha = 1
  beta = 0
  call dgemm_guru(transa,transb,nd,npols,npols,alpha,u,nd,umat0,npols,&
     beta,ucoefs,nd)

  allocate(pmat(npols,nfar_pols))
    
  do i=1,nfar_pols
    call koorn_pols(uvs(1,i),norder,npols,pmat(1,i))
  enddo

  transa = 'n'
  transb = 'n'
  call dgemm_guru(transa,transb,nd,nfar_pols,npols,alpha, &
        ucoefs,nd,pmat,npols,beta,uover,nd)
 
end subroutine oversample_fun_tri
!
!
!
!
!

subroutine surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
   u,ucoefs)

!
!f2py intent(in) nd,npatches,norders,ixyzs,iptype,npts,u
!f2py intent(out) ucoefs
!
!
!  This subroutine oversamples a function defined on a surface
!
!  input
!    nd - integer
!      number of functions
!    npatches - integer
!      number of patches on the surface
!    norders - integer(npatches)
!      order of discretization of patch i
!    ixyzs - integer(npatches+1)
!      starting location of points on patch i
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangle discretized using RV nodes
!    npts - integer
!      total number of points on the surface
!    u - real *8 (nd,npts)
!      functions at the discretization points on the surface
!
!  output
!    ucoefs - real *8 (nd,npts)
!       coeffs of funs in poly basis
!    
   

  implicit none
  integer nd,npatches,norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches),npts
  real *8 u(nd,npts),ucoefs(nd,npts)
  integer i,istart,npols

  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    if(iptype(i).eq.1) &
      call vals_to_coefs_tri(nd,norders(i),npols,u(1,istart),&
      ucoefs(1,istart))
  enddo

end subroutine surf_vals_to_coefs

!
!
subroutine vals_to_coefs_tri(nd,norder,npols,u,ucoefs)
!
!  This subroutine converts function to koornwinder expansion
!  coefs
!
!  input
!    nd - integer
!       number of functions
!    norder - integer
!       order of original discretization
!    npols - integer
!       number of discretization points corresponding to norder
!       npols = (norder+1)*(norder+2)/2
!    u - real *8 (nd,npols)
!       function values tabulated at original grid
!    
!  output
!    ucoefs - real *8 (nd,npols)
!       koornwinder expansion coefs
!


  implicit none
  integer norder,nd,npols
  real *8 u(nd,npols),ucoefs(nd,npols)
  real *8, allocatable :: umat0(:,:),uvs0(:,:),vmat0(:,:),wts0(:)
  character *1 transa,transb
  real *8 alpha, beta
  integer i



  allocate(umat0(npols,npols),vmat0(npols,npols),uvs0(2,npols))
  allocate(wts0(npols))


  call vioreanu_simplex_quad(norder,npols,uvs0,umat0,vmat0,wts0)


!
!  compute ucoefs
!

  transa = 'n'
  transb = 't'
  alpha = 1
  beta = 0
  call dgemm_guru(transa,transb,nd,npols,npols,alpha,u,nd,umat0,npols,&
     beta,ucoefs,nd)
 
end subroutine vals_to_coefs_tri 
!
!
!
!
!
!
subroutine get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,qwts)

!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srcvals
!f2py intent(out) qwts
!
!  this subroutine returns the quadrature weights
!  for integrating smooth functions on the surface
!
!
!  The surface is assumed to be discretizaed using
!  a collection of patches with the same order discretization
!  nodes
!
!  Input:
!    npatches - number of patches
!    norders(npatches) - order of discretization
!    ixyzs(npatches+1) - starting location of points on patch i
!    iptype(npatches) - type of patch
!    npts - total number of points on the surface
!    srcvals(12,npts) - surface info at the discretization nodes 
!
!  Output:
!    qwts - smooth quadrature weights at the discretization nodes
!
  implicit none
  integer norders(npatches),npatches,npts,npols,i,j,ipt
  integer ixyzs(npatches+1),iptype(npatches)
  real *8 srcvals(12,npts),qwts(npts),tmp(3)
  integer istart
  real *8, allocatable :: wts0(:)


  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    allocate(wts0(npols))
    if(iptype(i).eq.1) call get_vioreanu_wts(norders(i),npols,wts0)
    do j=1,npols
      ipt = istart + j-1
      call cross_prod3d(srcvals(4,ipt),srcvals(7,ipt),tmp)
      qwts(ipt) = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts0(j)
    enddo
    deallocate(wts0)
  enddo
end subroutine get_qwts

!
!
!
!
!
subroutine get_near_far_split_pt(nd,n,xyz,r,xyz0,nf,nn,xyz_f,xyz_n, &
    iind_f,iind_n)
!
!
!   This subroutine splits a given set of targets into
!   near points and far points with respect to a reference
!   point.
!
!   A target is labelled near if |t-xyz0|<=r \,, 
!   where xyz0 is the reference point, t is the target
!   
!
!   input
!     nd - integer
!        dimension of target array (must at least be 3, with the
!        first three components being the coordinates)
!     n - integer 
!       number of targets
!     xyz - real *8 (nd,n)
!       location of targets
!     r - real *8 
!       reference radius
!     xyz0 - real *8(3)
!        reference point
!     
!     output
!      nf - integer
!        number of far targets
!      nn - integer
!        number of near targets
!      xyz_f - real *8 (nd,n)
!        Collection of far targets. Note that
!        only the first (nd,nf) block will
!        be filled out
!      xyz_n - real *8 (nd,n)
!        Collection of near targets. Note that
!        only the first (nd,nn) block will
!        be filled out
!      iind_f - integer (n)
!         iind_f(i) denotes the location in xyz array
!         of the ith target in xyz_f.
!         Note that the first nf entries are filled out
!      iind_n - integer (n)
!         iind_n(i) denotes the location in xyz array
!         of the ith target in xyz_n.
!         Note that the first nn entries are filled out
!      
  implicit none
  integer nd
  integer n,nf,nn,iind_f(n),iind_n(n)
  real *8 xyz(nd,n),xyz0(3),r,rtmp,r2
  real *8 xyz_f(nd,n),xyz_n(nd,n)

  integer i,j

  r2 = r**2
  nf = 0
  nn = 0
  do i=1,n
    rtmp = (xyz(1,i)-xyz0(1))**2 + (xyz(2,i)-xyz0(2))**2 + &
       (xyz(3,i)-xyz0(3))**2
    if(rtmp.le.r2) then
      nn = nn +1

      do j=1,nd
        xyz_n(j,nn) = xyz(j,i)
      enddo

      iind_n(nn) = i
    else
      nf = nf +1
      do j=1,nd
        xyz_f(j,nf) = xyz(j,i)
      enddo
      iind_f(nf) = i
    endif
  enddo

end subroutine get_near_far_split_pt




subroutine get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
    ipatch_id,uvs_pts)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts
!f2py intent(out) ipatch_id,uvs_pts
!
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts,ipatch_id(npts)
  real *8 uvs_pts(2,npts)
  
  integer i,ip,j,npols


  do ip=1,npatches
    do j=ixyzs(ip),ixyzs(ip+1)-1
      ipatch_id(j) = ip
    enddo

    npols = ixyzs(ip+1)-ixyzs(ip)
    if(iptype(ip).eq.1) &
      call get_vioreanu_nodes(norders(ip),npols,uvs_pts(1,ixyzs(ip)))
  enddo


end subroutine get_patch_id_uvs




subroutine get_patch_distortion(npatches,norders,ixyzs,iptype,npts,&
    srccoefs,srcvals,qwts,pdis)

!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!f2py intent(in) qwts
!f2py intent(out) pdis
!
!
!   This subroutine estimates the patch distortion for a collection of 
!   patches:
!
!   For the right triangle, the patch distortion is defined as follows:
!
!   Let J denote the matrix whose columns are dxyz/du, and dxyz/dv
!
!   Let s1, s2 denote the largest and smallest singular values of the matrix
!     J 
!
!   Then the patch distortion is given by 
!      sqrt(\int_{T0} (s1/s2)^2 |det J| du dv/ \int_{T0} |det J| du dv)
!
!   
!   input:
!     npatches - number patches
!     norders(npatches) - order of discretization
!     ixyzs(npatches+1) - starting location of nodes on patch i
!     iptype(npatches) - type of patch
!           iptype = 1, triangular patch discretized with RV ndoes
!     npts - total number of discretization points
!     srcvals(12,npts) - xyz, dxyz/du,dxyz/dv, normals at all nodes
!     srccoefs(9,npts) - koornwinder expansion coeffs
!                         of geometry info
!     qwts(npts) - quadrature weights at discretization nodes
!
!     output:
!      pdis(npatches) - l2 distortion of each patch
!    
!
    
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts
  real *8 srccoefs(9,npts),srcvals(12,npts),qwts(npts),pdis(npatches)
  real *8 rtmp,srat
  integer i,istart,j,ipt,npols

  
  
  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    pdis(i) = 0
    rtmp = 0
    do j=1,npols
      ipt = istart+j-1
      call get_local_asprat(srcvals(4,ipt),srcvals(7,ipt),srat)
      rtmp = rtmp + qwts(ipt)
      pdis(i) = pdis(i) + srat*qwts(ipt)
    enddo
    pdis(i) = sqrt(pdis(i)/rtmp)
  enddo
end subroutine get_patch_distortion




subroutine get_local_asprat(du,dv,srat)
!
!
!  This subroutine computes the ratio of the square of the singular
!  values of the 3x2 jacobian matrix whose columsn are given by 
!  du and dv
!
!  J = [du dv]
!  J^{T}J = [|du|^2     <du,dv> ]
!           [<du,dv>      |dv|^2]
!
!  The singular values of J^{T}J are s1^2, s2^2 and are given by
!    (|du|^2 + |dv|^2 +- sqrt((|du|^2-|dv|^2)^2 + 4 <du,dv>^2))/2
!
!
!  input:
!    du - real *8 (3)
!       first column of jacobian matrix
!    dv - real *8 (3)
!       second column of jacobian matrix
!
!   output:
!      ratio: (s1/s2)^2
!


  implicit none
  real *8 du(3),dv(3),srat,a11,a12,a22,t,s1,s2

  a11 = du(1)**2 + du(2)**2 + du(3)**2
  a22 = dv(1)**2 + dv(2)**2 + dv(3)**2
  a12 = du(1)*dv(1) + du(2)*dv(2) + du(3)*dv(3)

  t = sqrt((a11-a22)**2+4*(a12**2))
  srat = (a11+a22+t)/(a11+a22-t)

end subroutine get_local_asprat








subroutine get_nximat(npatches,ixyzs,ixyzso,nximat)
  implicit none
  integer npatches,ixyzs(npatches+1),ixyzso(npatches+1),nximat,i

  nximat = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:nximat) 
  do i=1,npatches
    nximat = nximat+(ixyzso(i+1)-ixyzso(i))*(ixyzs(i+1)-ixyzs(i))
  enddo
!$OMP END PARALLEL DO  

end subroutine get_nximat




subroutine get_ximats(npatches,iptype,norders,ixyzs,novers,ixyzso, &
  nximat,ximats,iximat)
  
  implicit none
  integer npatches,iptype(npatches),norders(npatches),ixyzs(npatches+1)
  integer novers(npatches),ixyzso(npatches+1),iximat(npatches),nximat
  integer, allocatable :: nximats(:)
  integer i,j,k,l,npols,kpols
  
  real *8 ximats(nximat)

  allocate(nximats(npatches))

!$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,npatches
    nximats(i) = (ixyzso(i+1)-ixyzso(i))*(ixyzs(i+1)-ixyzs(i))
  enddo
!$OMP END PARALLEL DO
 
 iximat(1) = 1
 nximats(1) = nximats(1)+1
 call cumsum(npatches-1,nximats,iximat(2))
 


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(kpols,npols)

  do i=1,npatches
    if(iptype(i).eq.1) then
       kpols = ixyzs(i+1)-ixyzs(i)
       npols = ixyzso(i+1)-ixyzso(i)
       call koorn_oversamp_mat(norders(i),kpols,novers(i),npols, &
         ximats(iximat(i)))
    endif
  enddo
!$OMP END PARALLEL DO


end subroutine get_ximats
!
!
!
!
!
subroutine get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)
!
!------------------------
!  This subroutine computes the first fundamental form at
!  the discretization nodes on the surface.
!
!  The first fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{u} \cdot x_{u} & x_{u} \cdot x_{v} \\
!    x_{u} \cdot x_{v} & x_{v} \cdot x_{v} \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - ffform: double precision(2,2,npts)
!        first fundamental form at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: ffform(2,2,npts)
  integer i
  real *8 xuu,xvv,xuv

  do i=1,npts
    xuu = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
    xvv = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
    xuv = srcvals(4,i)*srcvals(7,i) + srcvals(5,i)*srcvals(8,i)+ &
       srcvals(6,i)*srcvals(9,i)
    ffform(1,1,i) = xuu
    ffform(2,1,i) = xuv
    ffform(1,2,i) = xuv
    ffform(2,2,i) = xvv
  enddo

  return
end subroutine get_first_fundamental_form
!
!
!
!
!

subroutine get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)
!------------------------
!  This subroutine computes the inverse of the first fundamental form at
!  the discretization nodes on the surface.
!
!  The first fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{u} \cdot x_{u} & x_{u} \cdot x_{v} \\
!    x_{u} \cdot x_{v} & x_{v} \cdot x_{v} \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - ffforminv: double precision(2,2,npts)
!        inverse of first fundamental form at the discretization nodes
!--------------------------
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches),npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: ffforminv(2,2,npts)
  integer i
  real *8 xuu,xvv,xuv,d

  do i=1,npts
    xuu = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
    xvv = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
    xuv = srcvals(4,i)*srcvals(7,i) + srcvals(5,i)*srcvals(8,i)+ &
       srcvals(6,i)*srcvals(9,i)
    d = xuu*xvv - xuv*xuv
    ffforminv(1,1,i) = xvv/d
    ffforminv(2,1,i) = -xuv/d
    ffforminv(1,2,i) = -xuv/d
    ffforminv(2,2,i) = xuu/d
  enddo

  return
end subroutine get_inv_first_fundamental_form
!
!
!
!
!
subroutine get_surf_grad(nd,npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,f,df)
!
!-----------------------------
!  Compute the surface gradient of scalar function
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - f: double precision (nd,npts)
!        function values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,2,npts)
!        surface gradient of the input functions
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),f(nd,npts)
  real *8, intent(out) :: df(nd,2,npts)
  real *8, allocatable :: ffforminv(:,:,:)

  allocate(ffforminv(2,2,npts))


  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  call get_surf_grad_fast(nd,npatches,norders,ixyzs,iptype,npts, &
    ffforminv,f,df)

  return
end subroutine get_surf_grad
!
!
!
!
!
!
subroutine get_surf_grad_fast(nd,npatches,norders,ixyzs,iptype,npts, &
  ffforminv,f,df)
!
!-----------------------------
!  Faster version of computing the surface gradient of scalar function,
!  with precomputed inverse of first fundamental form
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - ffforminv: double precision(2,2,npts)
!        inverse of first fundamental form at the discretization nodes
!    - f: double precision (nd,npts)
!        function values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,2,npts)
!        surface gradient of the input functions
!-----------------------------
!
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: ffforminv(2,2,npts),f(nd,npts)
  real *8, intent(out) :: df(nd,2,npts)
  real *8, allocatable :: fcoefs(:,:),dfuv(:,:,:)
  integer i,idim

  allocate(dfuv(nd,2,npts))
  call get_surf_uv_grad(nd,npatches,norders,ixyzs,iptype,npts,f,dfuv)


  do i=1,npts
    do idim=1,nd
      df(idim,1,i) = dfuv(idim,1,i)*ffforminv(1,1,i) + &
         dfuv(idim,2,i)*ffforminv(1,2,i)
      df(idim,2,i) = dfuv(idim,1,i)*ffforminv(2,1,i) + &
         dfuv(idim,2,i)*ffforminv(2,2,i)
    enddo
  enddo


  return
end subroutine get_surf_grad_fast
!
!
!
!
!
subroutine get_surf_uv_grad(nd,npatches,norders,ixyzs,iptype,npts,f, &
   dfuv)
!
!-----------------------------
!  Compute the uv gradient of scalar function
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - f: double precision (nd,npts)
!        function values at the discretization points
!
!  Output arguments:
!
!    - dfuv: double precision(nd,2,npts)
!        uv surface gradient of the input functions
!-----------------------------
!
  implicit none
  integer, intent(in) :: nd,npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches),npts
  real *8, intent(in) :: f(nd,npts)
  real *8, intent(out) :: dfuv(nd,2,npts)

  integer i,istart,npols



  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    if(iptype(i).eq.1) &
      call get_surf_uv_grad_tri(nd,norders(i),npols,f(1,istart),&
        dfuv(1,1,istart))
  enddo


  return
end subroutine get_surf_uv_grad


subroutine get_surf_uv_grad_tri(nd,norder,npols,f,dfuv)
!
!-----------------------------
!  Compute the uv gradient of scalar function on a single triangular
!  patch
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - norder: integer
!        order of discretization
!    - npols: integer
!        total number of points on the patch
!    - f: double precision (nd,npols)
!        function values at the discretization points
!
!  Output arguments:
!
!    - dfuv: double precision(nd,2,npols)
!        uv surface gradient of the input functions
!-----------------------------
!
  implicit none
  integer, intent(in) :: nd,norder,npols
  real *8, intent(in) :: f(nd,npols)
  real *8, intent(out) :: dfuv(nd,2,npols)
  real *8 fcoefs(nd,npols),pols(npols),ders(2,npols)
  real *8 uv(2,npols)
  integer i,idim,j

  call vals_to_coefs_tri(nd,norder,npols,f,fcoefs)
  call get_vioreanu_nodes(norder,npols,uv)

  do i=1,npols
    do idim=1,nd
      dfuv(idim,1,i) = 0
      dfuv(idim,2,i) = 0
    enddo
  enddo


  do i=1,npols
    call koorn_ders(uv(1,i),norder,npols,pols,ders)
    do j=1,npols
      do idim=1,nd
        dfuv(idim,1,i) = dfuv(idim,1,i) + ders(1,j)*fcoefs(idim,j)
        dfuv(idim,2,i) = dfuv(idim,2,i) + ders(2,j)*fcoefs(idim,j)
      enddo
    enddo
  enddo

  return
end subroutine get_surf_uv_grad_tri
!
!
!
!
!
subroutine get_surf_div(nd,npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,f,df)
!
!-----------------------------
!  Compute the surface divergence of a vector function
!
!  Input arguments:
!
!    - nd: integer
!        number of vector fields
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - f: double precision (nd,2,npts)
!        vector field values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,npts)
!        surface divergence of the input vector fields
!-----------------------------
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(in) :: f(nd,2,npts)
  real *8, intent(out) :: df(nd,npts)
  real *8, allocatable :: ffform(:,:,:)

  allocate(ffform(2,2,npts))


  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)

  call get_surf_div_fast(nd,npatches,norders,ixyzs,iptype,npts, &
    ffform,f,df)

  return
end subroutine get_surf_div
!
!
!
!
!


subroutine get_surf_div_fast(nd,npatches,norders,ixyzs,iptype,npts, &
  ffform,f,df)
!
!-----------------------------
!  Faster version of computing the surface divergence of a vector function
!  with precmoputed first fundamental form
!
!  Input arguments:
!
!    - nd: integer
!        number of vector fields
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - f: double precision (nd,2,npts)
!        vector field values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,npts)
!        surface divergence of the input vector fields
!-----------------------------
!
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts,nd
  real *8, intent(in) :: ffform(2,2,npts),f(nd,2,npts)
  real *8, intent(out) :: df(nd,npts)
  integer i,istart,npols

  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    if(iptype(i).eq.1) &
      call get_surf_div_tri(nd,norders(i),npols,ffform(1,1,istart), &
        f(1,1,istart),df(1,istart))
  enddo


  return

  return
end subroutine get_surf_div_fast
!
!
!
!
!
subroutine get_surf_div_tri(nd,norder,npols,ff,f,df)
!
!-----------------------------
!  Compute the surface divergence of a vector field on a single triangular
!  patch
!
!  Input arguments:
!
!    - nd: integer
!        number of functions
!    - norder: integer
!        order of discretization
!    - npols: integer
!        total number of points on the patch
!    - ff: double precision(2,2,npols)
!        first fundamental form at discretization nodes
!    - f: double precision (nd,2,npols)
!        vector field values at the discretization points
!
!  Output arguments:
!
!    - df: double precision(nd,npols)
!        surface divergecen of the input vector fields 
!-----------------------------
!
  implicit none
  integer, intent(in) :: nd,norder,npols
  real *8, intent(in) :: f(nd,2,npols),ff(2,2,npols)
  real *8, intent(out) :: df(nd,npols)
  real *8 fcoefs(nd,2,npols),pols(npols),ders(2,npols)
  real *8 ftmp(nd,2,npols)
  real *8 uv(2,npols),gg(npols)
  integer i,idim,j

  do i=1,npols
    gg(i) = sqrt(ff(1,1,i)*ff(2,2,i) - ff(1,2,i)*ff(2,1,i))
    do idim=1,nd
      ftmp(idim,1,i) = f(idim,1,i)*gg(i)
      ftmp(idim,2,i) = f(idim,2,i)*gg(i)
    enddo
  enddo

  call vals_to_coefs_tri(2*nd,norder,npols,ftmp,fcoefs)
  call get_vioreanu_nodes(norder,npols,uv)

  do i=1,npols
    do idim=1,nd
      df(idim,i) = 0
    enddo
  enddo


  do i=1,npols
    call koorn_ders(uv(1,i),norder,npols,pols,ders)
    do j=1,npols
      do idim=1,nd
        df(idim,i) = df(idim,i) + ders(1,j)*fcoefs(idim,1,j) + &
                                  ders(2,j)*fcoefs(idim,2,j)
      enddo
    enddo

    do idim=1,nd
      df(idim,i) = df(idim,i)/gg(i)
    enddo
  enddo

  return
end subroutine get_surf_div_tri
!
!
!
!

subroutine col_ind_to_patch_node_ind(npatches,ixyzs,ncol,col_ind, &
  patch_ind,node_ind)

  implicit none
  integer, intent(in) :: npatches,ixyzs(npatches+1),ncol
  integer, intent(in) :: col_ind(ncol)
  integer, intent(out) :: patch_ind(ncol),node_ind(ncol)

  integer i,j,j0,npols

  npols = ixyzs(2)-ixyzs(1)
  do i=1,ncol
    node_ind(i) = mod(col_ind(i),npols)
    patch_ind(i) = col_ind(i)/npols+1
    if(node_ind(i).eq.0) then
      node_ind(i) = npols
      patch_ind(i) = patch_ind(i) - 1
    endif
  enddo

!  call sorti(ncol,col_ind,iper)
!
!  do i=1,ncol
!    col_ind_sort(i) = col_ind(iper(i))
!    patch_ind(i) = 0
!    node_ind(i) = 0
!  enddo
!
!
!  j0 = 1
!  do i=1,ncol
!    do j=j0,npatches
!      if(ixyzs(j+1).gt.col_ind_sort(i)) then
!        patch_ind(iper(i)) = j
!        node_ind(iper(i)) = col_ind(iper(i))-ixyzs(j)+1
!        j0 = j
!        goto 1111
!      endif
!    enddo
! 1111 continue  
!  enddo
  

end subroutine col_ind_to_patch_node_ind
