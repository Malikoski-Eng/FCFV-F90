module stabilizationNS

use typeDefinitions

implicit none

contains
subroutine computaTauFace(uHatFace,NuEffFace,hdg,tau,n,iBC,typ,nsd)
  implicit none
  type(hdgdef), intent(in) :: hdg
  real(kind = dp), dimension(:), intent(in) :: uHatFace, n
  real(kind = dp), intent(in) :: NuEffFace
  real(kind = dp) :: tauA(nsd,nsd), tauD(nsd,nsd), eyeNsd(nsd,nsd), epsilon, epN
  real(kind = dp), intent(out) :: tau(nsd,nsd)
  integer, intent(in) :: iBC, typ, nsd
  integer ::  i, tauType

  eyeNsd = 0
  forall(i = 1:nsd) eyeNsd(i,i) = 1

  tauType = hdg%tauType
  epsilon = hdg%epsilon

  selectcase(tauType)
    case(1)
      tau = hdg%tauCST*eyeNsd
    case(2)
      tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
      tauA = maxval([epsilon, 2*sum(uHatFace*n)])*eyeNsd
      tau  = tauD + tauA
    case(4)
      tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
      tauA = maxval([epsilon,2*abs(sum(uHatFace*n))])*eyeNsd  
      tau  = tauD + tauA
    case(6)
      epN = epsilon/(abs(sum(uHatFace*n)))
      tauD = hdg%kappa(typ)*NuEffFace*eyeNsd 
      if (abs(sum(uHatFace*n)) < epsilon/2 .or. abs(sum(uHatFace*n)) == epsilon/2) then
        tauA = epsilon*eyeNsd 
      elseif( abs(sum(uHatFace*n)) > epsilon/2 .and. abs(sum(uHatFace*n)) < epsilon) then
        tauA = sign(1.0_dp,sum(uHatFace*n))*(epN*sum(uHatFace*n)*eyeNsd + (2.0-epN)*spread(uHatFace,2,nsd)*spread(n,1,nsd))
      else
        tauA = sign(1.0_dp,sum(uHatFace*n))*(sum(uHatFace*n)*eyeNsd + spread(uHatFace,2,nsd)*spread(n,1,nsd))
      endif
      tau = tauD + tauA
  end select
 end subroutine

 subroutine computaTauFaceGlobal(uHatFace,NuEffFace,hdg,tau,n,iBC,typ,nsd)
  implicit none
  type(hdgdef), intent(in) :: hdg
  real(kind = dp), dimension(:), intent(in) :: uHatFace, n
  real(kind = dp), intent(in) :: NuEffFace
  real(kind = dp) :: tauA(nsd,nsd), tauD(nsd,nsd), eyeNsd(nsd,nsd), epsilon, epN
  real(kind = dp), intent(out) :: tau(nsd,nsd)
  integer, intent(in) :: iBC, typ, nsd
  integer :: i, tauType

  eyeNsd = 0
  forall(i = 1:nsd) eyeNsd(i,i) = 1

  tauType = hdg%tauType
  epsilon = hdg%epsilon

  selectcase(tauType)
    case(1) 
      tau = hdg%tauCST*eyeNsd
    case(2) 
      if(iBC == ctt%iBC_Dirichlet) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = maxval([epsilon, 2*sum(uHatFace*n)])*eyeNsd
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Outlet) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = 0.0_dp*eyeNsd
         tau =  tauD + tauA
      elseif(iBC == ctt%iBC_Symmetry) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = epsilon*eyeNsd
         tau =  tauD + tauA
      elseif(iBC == ctt%iBC_Neumann) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = 0.0_dp*eyeNsd !maxval([epsilon, 2*sum(uHatFace*n)])*eyeNsd !
         tau =  tauD + tauA
      else
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = maxval([epsilon, 2*sum(uHatFace*n)])*eyeNsd
         tau = tauD + tauA
      endif
    case(4) 
      if(iBC == ctt%iBC_Dirichlet) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = maxval([epsilon,2*abs(sum(uHatFace*n))])*eyeNsd  
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Outlet .or. iBC == ctt%iBC_Symmetry) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = 0.0_dp*eyeNsd
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Neumann) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = 0.0_dp*eyeNsd
         tau = tauD + tauA
      else
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = maxval([epsilon,2*abs(sum(uHatFace*n))])*eyeNsd 
         tau = tauD + tauA
      endif
    case(6)
      epN = epsilon/(abs(sum(uHatFace*n)))
      if(iBC == ctt%iBC_Dirichlet) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        if (abs(sum(uHatFace*n)) < epsilon/2 .or. abs(sum(uHatFace*n)) == epsilon/2) then
          tauA = epsilon*eyeNsd 
        elseif( abs(sum(uHatFace*n)) > epsilon/2 .and. abs(sum(uHatFace*n)) < epsilon) then
          tauA = sign(1.0_dp,sum(uHatFace*n))*(epN*sum(uHatFace*n)*eyeNsd + (2.0-epN)*spread(uHatFace,2,nsd)*spread(n,1,nsd))
        else
          tauA = sign(1.0_dp,sum(uHatFace*n))*(sum(uHatFace*n)*eyeNsd + spread(uHatFace,2,nsd)*spread(n,1,nsd))
        endif
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Outlet .or. iBC == ctt%iBC_Symmetry) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = 0.0_dp*eyeNsd
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Neumann) then
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        tauA = 0.0_dp*eyeNsd
         tau = tauD + tauA
      else
        tauD = hdg%kappa(typ)*NuEffFace*eyeNsd
        if (abs(sum(uHatFace*n)) < epsilon/2 .or. abs(sum(uHatFace*n)) == epsilon/2) then
          tauA = epsilon*eyeNsd 
        elseif( abs(sum(uHatFace*n)) > epsilon/2 .and. abs(sum(uHatFace*n)) < epsilon) then
          tauA = sign(1.0_dp,sum(uHatFace*n))*(epN*sum(uHatFace*n)*eyeNsd + (2.0-epN)*spread(uHatFace,2,nsd)*spread(n,1,nsd))
        else
          tauA = sign(1.0_dp,sum(uHatFace*n))*(sum(uHatFace*n)*eyeNsd + spread(uHatFace,2,nsd)*spread(n,1,nsd))
        endif
        tau = tauD + tauA
      endif
  end select
 end subroutine

 subroutine computaTauFaceT(uHatFace,NuHatFace,hdg,tau,n,iBC,typ,SAmodel)
  implicit none
  type(hdgdef), intent(in) :: hdg
  real(kind = dp), dimension(:), intent(in) :: uHatFace, n
  real(kind = dp), intent(in) :: NuHatFace
  real(kind = dp) :: tauA, tauD, fn1, epsilon
  real(kind = dp), intent(out) :: tau
  integer, intent(in) :: iBC, typ, SAmodel
  integer :: tauType

  if(SAmodel == 1) then
    fn1 = 1.0_dp
  else
    fn1 = (turbCst%Cn1 + NuHatFace**3)/(turbCst%Cn1 - NuHatFace**3)
  end if

  tauType = hdg%tauType
  epsilon = hdg%epsilonSA

  selectcase(tauType)
    case(1)
      tau = hdg%tauCST
    case(2)
      tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
      tauA = maxval([epsilon, sum(uHatFace*n)])
      tau  = tauD + tauA
    case(4)
      tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
      tauA = maxval([epsilon,abs(sum(uHatFace*n))])
      tau  = tauD + tauA
    case(6)
      tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
      tauA = maxval([epsilon,abs(sum(uHatFace*n))])
      tau  = tauD + tauA
  end select
 end subroutine

 subroutine computaTauFaceGlobalT(uHatFace,NuHatFace,hdg,tau,n,iBC,typ,SAmodel)
  implicit none
  type(hdgdef), intent(in) :: hdg
  real(kind = dp), dimension(:), intent(in) :: uHatFace, n
  real(kind = dp), intent(in) :: NuHatFace
  real(kind = dp) :: tauA, tauD, fn1, epsilon
  real(kind = dp), intent(out) :: tau
  integer, intent(in) :: iBC, typ,SAmodel
  integer :: tauType

  if(SAmodel == 1) then
    fn1 = 1.0_dp
  else
    fn1 = (turbCst%Cn1 + NuHatFace**3)/(turbCst%Cn1 - NuHatFace**3)
  end if

  tauType = hdg%tauType
  epsilon = hdg%epsilonSA

  selectcase(tauType)
    case(1)
      tau = hdg%tauCST
    case(2)
      if(iBC == ctt%iBC_Dirichlet) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = maxval([epsilon, sum(uHatFace*n)])
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Outlet) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = 0.0_dp 
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Symmetry) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = epsilon
         tau = tauD + tauA         
      elseif(iBC == ctt%iBC_Neumann) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = 0.0_dp ! maxval([epsilon, sum(uHatFace*n)])
         tau = tauD + tauA
      else
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = maxval([epsilon, sum(uHatFace*n)])
         tau = tauD + tauA
      endif
    case(4)
      if(iBC == ctt%iBC_Dirichlet) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = maxval([epsilon,abs(sum(uHatFace*n))])
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Outlet .or. iBC == ctt%iBC_Symmetry) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = 0.0_dp 
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Neumann) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = 0.0_dp 
         tau = tauD + tauA
      else
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = maxval([epsilon,abs(sum(uHatFace*n))])
         tau = tauD + tauA
      endif
    case(6)
      if(iBC == ctt%iBC_Dirichlet) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = maxval([epsilon,abs(sum(uHatFace*n))])
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Outlet .or. iBC == ctt%iBC_Symmetry) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = 0.0_dp
         tau = tauD + tauA
      elseif(iBC == ctt%iBC_Neumann) then
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = 0.0_dp
         tau = tauD + tauA
      else
        tauD = hdg%kappaT(typ)*(1 + NuHatFace*fn1)/(turbCst%Sigma*problemParam%Rey)
        tauA = maxval([epsilon,abs(sum(uHatFace*n))])
         tau = tauD + tauA
      endif    
  end select
 end subroutine

end module
