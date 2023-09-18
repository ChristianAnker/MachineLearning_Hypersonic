!---------
!  user defined  variables and functions
!---------
module mod_user

  use mod_param
  use mod_geom
  use mod_var
  use mod_thermo
  use mod_trans
  
#if ibc
   use mod_solid
   use mod_ibc
   use mod_ibcpro
#endif

#if (nompi==0)
  use mod_mpi, only : procid
#endif
  use list_eos

#if amr
  use compreal_amr
  use amr_vars, only : nlevs,ibm_amr
  use multifab_ibm_aux
#endif

  use mod_detect 

  implicit none

  ! private variables within module
  real(wp), private  :: rho0,P0,T0,Ma0,mu0,u0
  real(wp), private, allocatable  :: Lbody(:)
  real(wp), private  :: Re0,E0,Twall,p_solid,rho_solid,pstag,xshock, Rey

  
  contains
  !------------------------------------------------
  !  This subroutine defines the problem
  !------------------------------------------------
  subroutine user_init

    implicit none
    ! local vars
    integer i,j,k
    real(wp), dimension(1:6)  :: temp_store
   
    if (ioproc) then
      write(*,*) ">>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<"
      write(*,*) "   (n)-bodies under uniform flow (AMR-ready)      "
      write(*,*) "      Example by SNM                                 "
      write(*,*) ">>>>>>>>>>>>>>                    <<<<<<<<<<<<<<<<<<<"
    end if

    !-- initialise parameters
    ! free stream
    Pr1 = 0.71_wp

    open(11, file = "input_file")
    read(11,*)(temp_store(j),j=1,6)
      Ma0  = temp_store(1)
      Re0  = temp_store(2)
      T0   = temp_store(3)
      ! wall
      Twall= temp_store(4)

      P0   = temp_store(5)
      xshock = temp_store(6)

    close(11) 

    ! Reynolds per unit length
    Rey = Re0

    !---------------------------
    ! boundary layer estimate

    !delta0 = 4.91_wp* Lref /sqrt(Re0)

    ! soem reference values
    rho0 =  eos_rho(P0,T0)   
    U0   =  Ma0*eos_sound(rho0,P0)
    mu0  =  rho0*U0/Re0 
    ! overwrite reference value
    visc_number = mu0

    ! overwrite visc_ref
    visc_ref = mu0
    E0   =  eos_E(rho0,P0)

    !----  values to starting filling the solid (stagnation velocity)
    pstag     = P0 + r12*U0**2
    p_solid   = P0
    rho_solid = eos_rho(p_solid,Twall)   
 
    if (ioproc) then
      write(*,*) " -----------------------------------"
      write(*,*) "         Flow parameters            "
      write(*,*) " rho0 U0 =",rho0,U0
      write(*,*) " Ma0 Re0 =",Ma0,Re0
      write(*,*) " Ref viscosity (298 K) =",visc_ref
      write(*,*) " Imposed viscosity = ",mu0," visc(T)=",visc_f(T0) 
      write(*,*) " -----------------------------------"
      write(*,*) "          Grid parameters           "
      !write(*,*) "  delta0/dy =",delta0*invDely      
      write(*,*) "  Re (dx) =",Rey*Delx 
#if amr    
      write(*,*) "  Re (dx finest) =",Rey*Delx/(2**(max_levs-1))
      write(*,*) "  dt must be less than = ", 0.3_wp * Delx/(2**(max_levs-1)) / U0
#endif
      write(*,*) " ------------------------------------"
    end if   

    ! bodies positions
    
    
    !--- call initial function
    do i=nxlo,nxup
    do j=nylo,nyup
    do k=nzlo,nzup
      call initial_func(x(i),y(j),z(k),q(i,j,k,:))
    end do
    end do
    end do

    !--- assign bc 
    call user_bc
	
    ! call intial solid
    allocate(Lbody(1))
    ibm%level = 1  !(for non-AMR init)
    call initial_solid(ibm,sld,lvl)

    ! assemble body
    call assemble_solid(ibm,sld,lvl)

    ! find fluid point neighbours to solid
    call neigh_solid(ibm,ngh)  !<<<<<<<<***

!=================================
#if amr
    if (ioproc) write(*,*) "   .......   AMR ACTIVATED ..... "
    !  amr parameters
    max_grid_size = 128  ! default 64
    amr_buf_width = 4     ! default 2 points between levels
    cluster_minwidth = 16  ! default 16  
    cluster_blocking_factor = 2
    cluster_min_eff = 0.7d0  !default 0.7
#endif
	
  ! initialise surfdata tiles


#if debug
    write(*,*) " [user_int] end"
#endif

  open(unit=56, file='global.csv')
  write(56,*) "Ma0, Re0, mu0, P0, T0, Twall, rho0, U0, E0, Reh, dTmax"
  write(56,*) Ma0, Re0, mu0, P0, T0, Twall, rho0, U0, E0, Rey*Delx/(2**(max_levs-1)), 0.3_wp * Delx/(2**(max_levs-1)) / U0
  close(56)
  
  end subroutine user_init
  !----------------------------------------------------------
  !  This subroutine initialise an ibm solid
  !----------------------------------------------------------
  subroutine initial_solid(ibm_n,sld_n,lvl_n, dx_in)

    implicit none

    ! arguments
    type(solid), intent(inout) :: ibm_n
    real(wp),   intent(inout) :: sld_n(nxlo:nxup,nylo:nyup,nzlo:nzup)
    real(wp),   intent(inout) :: lvl_n(nxlo:nxup,nylo:nyup,nzlo:nzup)
    
    real(wp), optional :: dx_in

    ! local vars
    integer :: ib, iflag
    type(body) :: body_n
    real(wp) :: xpos_n, ypos_n
    INTEGER :: points

#if debug
    write(99,*) "[initial_solid]"
    write(99,*) " present dx=?",present(dx_in)
#endif   

    if (present(dx_in)) then
      iflag = 1
    else  
      iflag = 0
    end if
    
    xpos_n = r12 * lx   !
    ypos_n = 0.0_wp
    points = 2998

    if (present(dx_in)) then
      call geo_read("shape.geo",points,xpos_n,ypos_n,0.0_wp,body_n,sld_n,lvl_n,.false.,dx_in)
    else
      call geo_read("shape.geo",points,xpos_n,ypos_n,0.0_wp,body_n,sld_n,lvl_n,.false.)
    end if

    ! add body to structure 
    call addbody_solid(body_n,ibm_n,iflag=0) !<------


  end subroutine initial_solid
  !----------------------------------------------------------
  !  This subroutine defines the initial conditions
  !  called from user_init and user_init_amr
  !----------------------------------------------------------
  subroutine initial_func(xx,yy,zz,phi)

    implicit none
    ! arguments
    real(wp), intent(in) :: xx,yy,zz
    real(wp), intent(inout) :: phi(nvar)

    ! solid marker  (set to 0 outside body)
    phi(nvki) = 0.0_wp
  
    ! post-shock state
    if ( xx .gt.  xshock ) then
      phi(nvux) = 0.0_wp
      phi(nvpr) = P0
      phi(nvtp) = 2.0_wp*T0
    else
      phi(nvux) = U0
      phi(nvpr) = p0     
      phi(nvtp) = T0
    end if   

    phi(nvuy) = 0.0_wp    
   
    call satisfy_eos(phi)
   
  end subroutine initial_func
  !----------------------------------------------------------
  !  This defines the boundary conditions (use with bcsimple)
  !  specify a and b for complex mix/match BC
  !----------------------------------------------------------      
  subroutine user_bc

    implicit none
    ! local vars
    integer :: nv,i

    !===================================
    ! inflow            i=1
    !===================================
    qbc_i1(nvpr) = P0
    qbc_i1(nvux) = U0
    qbc_i1(nvuy) = 0.0_wp
    qbc_i1(nvtp) = T0
    call satisfy_eos(qbc_i1)
    Mi1 = Ma0
    
    ! initial set-up  
    do nv = 1,nvar      
      qbc_i1unst_var(0,:,:,nv)    = qbc_i1(nv)   !
      qbc_i1unst_der(0,:,:,nv)    = 0.0_wp       ! dq/dx
    enddo
  
  end subroutine user_bc
  !----------------------------------------------------------
  !  This subroutine is called at the end  of the program
  !----------------------------------------------------------
  subroutine user_end
 
#if amr    
    use amr_vars, only : mla
    use amr_vars, only : phi_new
    use multifab_operations
#endif
    
    implicit none
    ! local vars
    real(wp) :: aux
    integer i,j,k, nv


    
  end subroutine user_end
  !----------------------------------------------------------
  !  This subroutine is called in main to use unsteady inflow
  !----------------------------------------------------------      
  subroutine user_inflow(time)
    implicit none
    ! arguments
    real(wp), intent(in) :: time
    ! local vars
    real(wp) :: dragP(ndim),drag(ndim)
    character(len=10) :: file_id
    character(len=50) :: file_name
    integer :: num
 
   
  end subroutine user_inflow
  !---------------------------------------------------------
  subroutine user_restart
    implicit none
    ! local varstste
    return  
  end subroutine  user_restart
  !---------------------------------------------------------
  subroutine user_source
    implicit none
    ! local vars
    return  
  end subroutine  user_source
  !--------------------------------------------------------
  !  this subroutine sets-up the farfield conditions
  !  (if option selected)
  !--------------------------------------------------------
  subroutine  user_farfield
    implicit none
    return
  end subroutine  user_farfield
  !-------------------------------------------------------
  subroutine user_output
    implicit none
    return
  end subroutine user_output  

  
  !===========================================================================
  !  AMR subroutines  
  !===========================================================================
  !----------------------------------------------------------
  !  This subroutine initilaiose solution for all levels
  !----------------------------------------------------------
  subroutine user_init_amr(f,dxx,dyy,dzz,sld_m,ifab,ilev)

    implicit none

    real(wp), intent(in) :: dxx,dyy,dzz  ! grid spacing
    real(wp) :: f(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    logical, optional, intent(in) :: sld_m(nxlo:nxup,nylo:nyup,nzlo:nzup)
    integer, optional, intent(in) :: ifab,ilev

    ! local vars
    integer  i,j,k
    real(wp) :: xx,yy,zz

    zz =0.d0

    do i=nxlo,nxup
    do j=nylo,nyup
    do k=nzlo,nzup
      ! compute xx,yy,zz
      xx =  (dble(i)+0.5d0) * dxx  ! xx(0) = dx/2   BoxLib notation
      yy =  (dble(j)+0.5d0) * dyy  
      call initial_func(xx,yy,zz,f(i,j,k,:))
    end do
    end do
    end do

  end subroutine user_init_amr
  
  !----------------------------------------------------------
  subroutine amr_pltuser(f,fplt,lev)

    implicit none

    real(wp) ,intent(in)   :: f(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    real(wp) ,intent(inout):: fplt(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar_plt)
    integer, intent(in) :: lev

   end subroutine  amr_pltuser
  !----------------------------------------------------------
  subroutine user_source_amr(f,dts,rhs)

    implicit none

    real(wp), intent(in) :: dts ! time step
    real(wp), intent(in) :: f(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    real(wp), intent(inout) :: rhs(nxlo:nxup,nylo:nyup,nzlo:nzup,nvarsolve)
    
  end subroutine  user_source_amr

  !----------------------------------------------------------
  !  This subroutine selects the criteria for refinement
  !----------------------------------------------------------
  subroutine user_tagbox(tagbox,f,dxx,dyy,dzz,ilev)

    use mod_param
    use mod_geom
    use mod_var
    use mod_thermo
  
#if amr  
    use compreal_amr, only : ist,ien,jst,jen,kst,ken
#endif

    implicit none

    real(wp), intent(in) :: dxx,dyy,dzz  ! grid spacing
    integer, intent(in) :: ilev
    logical :: tagbox(ist:ien,jst:jen,kst:ken)               
    real(wp) :: f(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)

    ! local vars
    integer  i,j,k
    real(wp) :: xx,yy,zz,rhop

    zz = 0.d0
  
    tagbox(:,:,:)  = .false. ! no tag any cell by default 

    select case(ilev)

    case (1)
    ! level 1 of refinement ---------------------------------------------------------------------
      do i=ist,ien
      do j=jst,jen
      do k=kst,ken
        rhop = f(i,j,k,nvrh)
        ! density variation  drho/rho > 20 %   
        tagbox(i,j,k) = detect_jump(i,j,k,f(:,:,:,nvrh),rhop,0.1d0)
      end do
      end do
      end do
    case (2)
    ! level 2 of refinement ---------------------------------------------------------------------
      do i=ist,ien
      do j=jst,jen
      do k=kst,ken
        rhop = f(i,j,k,nvrh)
        ! density variation  drho/rho > 30 %     
        tagbox(i,j,k) = detect_jump(i,j,k,f(:,:,:,nvrh),rhop,0.2d0) 
      end do
      end do
      end do
    case default

    ! level >= 3 of refinement -----------------------------------------------------------------
      do k=kst,ken  
      do j=jst,jen
      do i=ist,ien
        rhop = f(i,j,k,nvrh)
        ! density variation  drho/rho > 40 %     
        tagbox(i,j,k) = detect_jump(i,j,k,f(:,:,:,nvrh),rhop,0.3d0) 
      end do
      end do
      end do

    end select

    do i=ist,ien
      xx =  (dble(i)+0.5d0) * dxx  
      if (xx .lt. 0.95d0*xshock) tagbox(i,:,:) = .false.
    end do  

    do k=kst,ken  
    do j=jst,jen
    do i=ist,ien
      tagbox(i,j,k) = tagbox(i,j,k) .or. detect_closesolid(i,j,k,f(:,:,:,nvki))
    end do
    end do
    end do    

  end subroutine user_tagbox

  subroutine user_bcamr(f,time_in,dxx,dyy,dzz,abc,bbc)

    implicit none

    real(wp), intent(in)    :: time_in         ! time
    real(wp) ,intent(in)    :: dxx,dyy,dzz  ! mesh spacing
    real(wp) ,intent(in)    :: f(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    real(wp) ,intent(inout) :: abc(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)
    real(wp) ,intent(inout) :: bbc(nxlo:nxup,nylo:nyup,nzlo:nzup,nvar)

    return
  end subroutine user_bcamr
  !===========================================================================
  ! end amr subroutines
  !===========================================================================

  !===========================================================================
  !  IBC subroutines  (will only be compiled if IBC active)
  !===========================================================================

#if ibc
  subroutine user_solidbc(type_bc,phi_bc,dphi_bc,xw,yw,zw)

    implicit none

    integer,  intent(inout) :: type_bc(nvar)             ! type of bcconditions
    real(wp), intent(inout):: phi_bc(nvar),dphi_bc(nvar) ! surface values and derivatives
    real(wp), intent(in) ::   xw,yw,zw                   ! surface coordinates

#if debug
    real(wp) :: aux
    aux = xw + yw + zw
#endif

    ! Initialize all BC to 0-Dirichlet (default)
    type_bc(:) = 1
    phi_bc(:)  = 0.0_wp
    dphi_bc(:) = 0.0_wp

    ! BL approximation dP/dy = 0 .................................. 3
    type_bc(nvpr) = 0
    dphi_bc(nvpr) = 0.0_wp

    ! iso-thermal at Twall
    type_bc(nvtp) = 1
    phi_bc(nvtp)  = Twall

  end subroutine user_solidbc

 !------------------------------------------------
 ! wall/solid values (only at t=0)
 !------------------------------------------------

    subroutine user_fixibc(phi,x1,y1,z1)

   implicit none

   real(wp), intent(in) :: x1,y1             ! x,y
   real(wp), optional,intent(in) :: z1         ! z
   real(wp), intent(inout) :: phi(nvar)      ! vars value at x,y

#if debug
   real(wp) :: aux
   aux = x1 + y1 + z1
#endif


   ! this subroutine
   ! u=v=0   P=p_solid, T=Twall
   !--------------
    phi(nvpr) = p_solid
    phi(nvux) = 0.0_wp
    phi(nvuy) = 0.0_wp
    phi(nvtp) = Twall

    ! given P, T update rho and E 
    call satisfy_eos(phi)
    

   end subroutine user_fixibc

!------------------------------------------------
!  This subroutine defines solid values
!------------------------------------------------
  subroutine user_fixsolid(phi,xx,yy,zz)

    implicit none

    real(wp), intent(inout) :: phi(nvar)
    real(wp), intent(in) :: xx,yy,zz

    !--------------
    phi(nvpr) = p_solid
    phi(nvux) = 0.0_wp
    phi(nvuy) = 0.0_wp
    phi(nvtp) = Twall
    ! given P, T update rho and E 
    call satisfy_eos(phi)
    
    ! is a solid point
    phi(nvki) = 1.0_wp
    
  end subroutine user_fixsolid

#endif 


!---------------------------------------------------------------------------------------------------------------------------
  end module mod_user




