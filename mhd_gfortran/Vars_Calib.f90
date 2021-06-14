!---------------------------------------------------------------------------------------------
! modulo de declaracao de variaveis do metodo de calibracao
!---------------------------------------------------------------------------------------------
module vars_calib
implicit none
save

integer, parameter:: npar=10 ! numero maximo de parametros a serem otimizados
real*8,allocatable:: a(:),blow(:),bupp(:),x(:,:),xnstd(:),bound(:),cx(:,:),cf(:)
real*8:: pcento
real,allocatable:: csi_ref(:) ! valor de referencia do csi
integer,allocatable:: solo(:),veg(:) ! indica tipo de vegetacao e de solo no arquivo Blocos.hig 
real,allocatable:: qjcal(:,:,:) !armazena as vazoes de jusante das bacias a montante da bacia a ser calibrada
real mu_ref,cs_ref,alpha_ref,tsub_ref,cb_ref
real ths(13),thr(13),psi(13),ksat(13),b(13),vpar(npar)
integer ibcal,pos(npar),nqcal
integer,allocatable:: iqcal(:) ! vetor com as celulas dos exutorios das bacias a montante e as bacia a ser calibrada
integer nopt,maxn,kstop,iseed,ngs,npg,nps,nspl,mings,npt
integer iniflg,iprint
character*10 xname(npar)
data xname /'    d1','    d2','    d3','  dkss',' dtsub','   dmu','dalpha', &
&'  dcsi','   dcs','   dcb'/

end module vars_calib