!---------------------------------------------------------------------------------------------
!Subrotina que calcula alguns parâmetros da célula e do rio
!---------------------------------------------------------------------------------------------

subroutine parcel
use vars_main, only:nc,acel,tcon,hcel,lcel,tkb,dtp,cb,ibac,tks,cs,adren,bc1,brio,bc1,bc2,bc3,brio,qref,qmesp
implicit none
integer ic
real xlado !comprimento do lado da célula

!calcula tempo de concentracao na celula usando kirpich
do ic=1,nc
	xlado=acel(ic)**0.5 !estimativa do comprimento
	tcon(ic)=3600.*((0.868*xlado**3.0)/(hcel(ic)-lcel(ic)))**0.385	! kirpich em segundos
enddo

!calcula os coeficientes do rls da celula
do ic=1,nc
	tkb(ic)=exp(-dtp/(cb(ibac(ic))))       !cb definido em segundos em lesolo.f90 
	tks(ic)=exp(-dtp/(cs(ibac(ic))*tcon(ic)))   !cs definido em lesolo.f90	
enddo

!calcula largura do rio e vazao de referencia
do ic=1,nc
	if(adren(ic).gt.1.0)then	 
	    if(bc1 == 0.) then
	        brio(ic)=bc2*adren(ic)**bc3	 !largura do rio em m
	    else
		    brio(ic)=bc1*adren(ic)*adren(ic)+bc2*adren(ic)+bc3	 !largura do rio
		endif
		qref(ic)=qmesp*adren(ic) !vazao de referencia em m3/s
	else !areas fonte, sem rio, nao interessam
		brio(ic)=0.0
    	qref(ic)=0.0
	endif
enddo
!fim do calculo da largura e vazao de referencia
return
end