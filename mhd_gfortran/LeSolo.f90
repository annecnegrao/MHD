!---------------------------------------------------------------------------------------------
!subrotina que leh os parametros de 13 tipos de solo do triangulo textural
!ParSolo.hig
!estes parametros sao multiplicados pelos coeficientes definidos no arquivo ajuste.hig
!---------------------------------------------------------------------------------------------
subroutine lesolo
use vars_main, only: filsolo,filbacia,filajuste,dir_dados,smax,kss,tsub,alpha,csi,ssmax,mu,srmax,neta, &
                   & cs,cb,secrbs,nu,secrjah,dtp,secrjal,bc1,bc2,bc3,qmesp,qesp,qcons,nb,d1,d2,d3
use vars_calib, only: psi,b,ths,thr,ksat,csi_ref,solo,tsub_ref,mu_ref,alpha_ref,cs_ref
implicit none
integer ib,is,iu
real dkss,dtsub,dmu,dalpha,dcsi,dcs,dcb !parametros usados no ajuste
real psicrh,psicrl,bbc
character texto*20

!inicilizando variaveis
smax=0. ; ssmax=0. ; srmax=0.
kss=0. ; tsub=0. ; mu=1.
neta=1. ; cs=0. ; cb=0. ; csi=0.

open(filsolo,file=dir_dados//'entrada/ParSolo.hig',status='old')	
! le os parametros de Brooks Corey de 13 solos do traingulo textural do USDA
! ksat - condutividade hidraulica saturada (mm/h)
! b - parametro da curva de retencao
! psib - pressao de entrada do ar (kpa)
! ths - umidade volumetrica na saturacao
! thr - umidade volumetrica residual
read(filsolo,*)
read(filsolo,*)
read(filsolo,'(a15,5g10.0)')(texto,ksat(is),psi(is),b(is),ths(is),thr(is),is=1,13)
close(filsolo)
allocate (csi_ref(nu))

do iu=1,nu
    if(solo(iu) /= 14) then ! solo 14 = agua
        ! potencial de entrada do ar e b de brooks corey
        bbc=1./b(solo(iu))
        ! define csi = capacidade de campo = 10 kpa
        csi_ref(iu)=thr(solo(iu))+(ths(solo(iu))-thr(solo(iu)))*(min(psi(solo(iu))/10.,1.))**bbc       
        csi_ref(iu)=(csi_ref(iu)-thr(solo(iu)))/(ths(solo(iu))-thr(solo(iu)))
        ! limite de agua disponivel para evaporacao de solo nu
	    secrbs(iu)=csi_ref(iu)
        ! define fatores de estresse do modelo de jarvis
        psicrh=20 ! etp = 5 mm/dia psi_c= -20 kpa
        psicrl=80 ! etp = 1 mm/dia psi_c= -80 kpa
        secrjah(iu)=(min(psi(solo(iu))/psicrh,1.))**bbc ! saturacao efetiva critica de inicio de estresse para etp=5mm/dia
        secrjal(iu)=(psi(solo(iu))/psicrl)**bbc ! saturacao efetiva critica de inicio de estresse para etp=1mm/dia
        ! transforma os valores para calculo linear
        secrjah(iu)=(secrjah(iu)-secrjal(iu))/((5.0-1.0)*dtp/86400)
        secrjal(iu)=secrjal(iu)-secrjah(iu)*1.*dtp/86400
    endif
enddo

open(filbacia,file=dir_dados//'entrada/ParBacia.hig',status='old')
! le parametros de escoamento in-cell de cada sub-bacia
read(filbacia,*)
read(filbacia,*) tsub_ref    !transmissividade maxima de referencia (m2/dia)
read(filbacia,*) mu_ref      !coeficiente de decaimento da transmissividade com a profundidade
read(filbacia,*) alpha_ref   !fator de anisotropia da camada superriro do solo
read(filbacia,*) cs_ref      !coeficiente do reservatorio linear simples do fluxo sup e subsuperficial
read(filbacia,*) bc1,bc2,bc3 !coeficientes para a largura do rio (m)
read(filbacia,*) qmesp       !vazao media especifica (m3/s/km2)
read(filbacia,*)
read(filbacia,'(a10,3g10.0)')(texto,cb(ib),qesp(ib),qcons(ib),ib=1,nb)
close (filbacia)

! le fatores de ponderacao de cada sub-bacia
open(filajuste,file=dir_dados//'entrada/ParAjuste.hig',status='old')
read(filajuste,*) 
do ib=1,nb
    ! d1	:espessura da camada superior do solo
    ! d2	:espessura da segunda camada do solo
    ! d3    :espessura da terceira camada do solo
    ! dkss	:fator de ponderacao do parametro kss 
    ! dtsub	:fator de ponderacao do parametro tsub
    ! dmu	:fator de ponderacao do parametro mu
 	! dalpha:fator de ponderacao do parametro alpha
 	! dcsi  :fator de ponderacao do parametro csi
 	! dcs   :fator de ponderacao do parametro cs
 	! dcb   :fator de ponderacao do parametro cb
 	read(filajuste,'(a10,10g10.0)') texto,d1(ib),d2(ib),d3(ib),dkss,dtsub,dmu,dalpha,dcsi,dcs,dcb
    do iu=1,nu
        if(solo(iu) /= 14) then
            ! calcula os parametros de cada sub-bacia
            ssmax(ib,iu)=(ths(solo(iu))-thr(solo(iu)))*1000*d1(ib) ! armazenamento maximo camada superior do solo
	        srmax(ib,iu)=(ths(solo(iu))-thr(solo(iu)))*1000*d2(ib) ! armazenamento maximo em mm na segunda camada
            smax(ib,iu)=(ths(solo(iu))-thr(solo(iu)))*1000*d3(ib) ! armazenamento maximo em mm na terceira camada
            alpha(ib,iu)=alpha_ref*dalpha ! coeficiente de anisotropia
            kss(ib,iu)=(ksat(solo(iu))*24./1000.)*dkss ! condutividade hidraulica saturada camada superior m/dia
            neta(iu)=2.5+2.*b(solo(iu)) ! coeficiente da  condutivdade hidraulica nao saturada
            tsub(ib,iu)=tsub_ref*dtsub ! transmissividade maxima
            mu(ib,iu)=mu_ref*dmu ! coeficiente de decaimento da transmissividade com a profundidade
	        csi(ib,iu)=min(csi_ref(iu)*dcsi,1.) ! capacidade de campo
        endif
    enddo
	cs(ib)=cs_ref*dcs
    cb(ib)=cb(ib)*dcb*86400. ! cb calculado em dias e transformado em segundos
enddo

close (filajuste)

write(*,*)
write(*,*) ' leu parametros de solo'

return
end