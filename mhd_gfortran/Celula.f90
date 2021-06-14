!---------------------------------------------------------------------------------------------
! Subrotina que comanda o loop das celulas
!---------------------------------------------------------------------------------------------
subroutine celula
use vars_main, only:jdia,imes,idia,iano,ncdom,icdom,ibac,icell,evq,afp,vfp,pre,e0ag,vfpt,zfp,zfpb, &
&zfpt,afpt,nu,puso,smax,icalib,sssu,srad,ssub,ec,ett,etp,asat,nbvfo,dsup,dss,dsub,acel,dtp,vsub,vssu,&
&vsup,tkb,qsub,vsub,tks,qssu,vssu,qsup,preb,it,kcb,ecb,ettb,etpb,sssub,ssubb,sradb,asatb,dsupb,dssub,dsubb
implicit none
integer ic,ic2,iu,ib,ifp,k
real vbx,vix,vsx,ppu,fna
integer julday
!valores medios na celula de armazenamentos nas camadas, da evaporacao do dossel e total. 
real sssuc,ssubc,sradc,ecc,ettc,etpc,asatc 

!verifica a que mes corresponde o dia juliano
call caldat(jdia,imes,idia,iano)
jdia=jdia-julday(1,1,iano)+1

do ic2=1,ncdom	 	!loop nas celulas do dominio
    ic=icdom(ic2)
	ecc=0.0 !zera variaveis temporarias
	ettc=0.0
	sssuc=0.0
	ssubc=0.0
	sradc=0.0
	etpc=0.
	asatc=0.
	ib=ibac(ic)
    ifp=icell(ic)
	evq(ic)=0.0 !evaporacao direta da superficie liquida
	fna=1. !fracao da celula sem alagamento
	! area alagada e tratada com um uso separado, com parametros iguais a agua
	if(ifp > 0) then ! celula de planicie
        if(afp(ifp) > 0.001) then
            ! area alagada, a precipitacao menos evaporacao altera o volume armazenado, e nao gera escoamento
            vfp(ifp)=vfp(ifp)+afp(ifp)*(pre(ic)-e0ag)/1000.
            ! atualiza a altura e a area de alagamento devido as trocas verticais
            if(vfp(ifp) <= vfpt(ifp,1)) then
                zfp(ifp)=zfpb(ifp)
                afp(ifp)=0.
            else
                zfp(ifp)=zfpt(ifp,50)
                afp(ifp)=1.
                do k=2,50
                    if(vfp(ifp) < vfpt(ifp,k-1)) then
                        afp(ifp)=afpt(ifp,k-1)+(vfp(ifp)-vfpt(ifp,k-1))*(afpt(ifp,k)-afpt(ifp,k-1))/(vfpt(ifp,k)-vfpt(ifp,k-1))
                        zfp(ifp)=zfpt(ifp,k-1)+(vfp(ifp)-vfpt(ifp,k-1))*(zfpt(ifp,k)-zfpt(ifp,k-1))/(zfpt(ifp,k)-zfpt(ifp,k-1))
                        exit
                    endif
                enddo
            endif
            fna=(1.-afp(ifp)/100) ! fracao de area nao alagada
		    ! considera que a area alagavel nao gera escoamento sup, subs ou de base
		    ! logo, dsup=dss=dsub=dd=0.0
        endif
    endif

	do iu=1,nu !loop dos usos
	    ! ajusta proporcao de uso no tempo t corrigido pela area alagada
        ppu=puso(ic,iu)/100.*fna
		if(puso(ic,iu) < 0.001) cycle !nao tem este uso nesta celula, passa para o proximo uso
		call evaporacao(ic,iu)  !calcula la radiacion liquida e interceptacao	
	
		if(smax(ib,iu).gt.0.001) then !faz o balanco hidrico do solo, se for area coberta por agua (bloco agua, smax=0.0)
			call escoamentos(ic,iu,ib)
			call transpiracao(ic,iu,ib)
    		if(icalib == 0) then ! esta simulando
			    sssuc=sssuc+sssu(ic,iu)*ppu !armazenamento medio na camada superior na celula
			    sradc=sradc+srad(ic,iu)*ppu !armazenamento na camada intermediaria na celula
			    ssubc=ssubc+ssub(ic,iu)*ppu !armazenamento na camada inferior na celula
			    ecc=ecc+(ec(ic,iu))*ppu !interceptacao média na celula
			    ettc=ettc+ett(ic,iu)*ppu !evapotranspiracao total média na celula
			    etpc=etpc+etp(ic,iu)*ppu !transpiracao potencial media na celula
			    asatc=asatc+asat(ic,iu)*ppu !area saturada media na celula
            elseif(ib == nbvfo) then
			    ettc=ettc+ett(ic,iu)*ppu !evapotranspiracao total média na celula            
            endif			
		else ! corpo de agua permanente de superficie - cte: reservatorio ou rio de grande porte
			dsup=pre(ic) ! gera escoamento superficial
			dss=0.0
			dsub=0
			evq(ic)=(e0ag*1000.*acel(ic)*ppu)/dtp !evaporacao direta das superficies líquidas em m3/s
			ettc=ettc+e0ag*ppu
		endif

		!corrige as unidades; multiplica pela area da celula (km2); multiplica pela proporcao de uso
		!os valores de dsup, dssu e dsub estao em mm - converte para m3	
		!atualiza volumes (em m3)
		vsub(ic)=vsub(ic)+dsub*acel(ic)*ppu*1000.
		vssu(ic)=vssu(ic)+dss*acel(ic)*ppu*1000.
		vsup(ic)=vsup(ic)+dsup*acel(ic)*ppu*1000.
	enddo !fim do loop dos usos

	!calcula vazoes das celulas usando o rls
	vbx=vsub(ic)*tkb(ic)
	qsub(ic)=(vsub(ic)-vbx)/dtp
	vsub(ic)=vbx
	vix=vssu(ic)*tks(ic)
	qssu(ic)=(vssu(ic)-vix)/dtp
	vssu(ic)=vix
	vsx=vsup(ic)*tks(ic)
	qsup(ic)=(vsup(ic)-vsx)/dtp
	vsup(ic)=vsx
	
    !guarda valores medios para cada bacia
    !ib é o numero da bacia a qual pertence a celula ic
	if(icalib == 0) then ! esta simulando
	    preb(ib,it)=preb(ib,it)+pre(ic)/kcb(ib) !precipitacao
	    ecb(ib,it)=ecb(ib,it)+ecc/kcb(ib) !lamina interceptada
	    ettb(ib,it)=ettb(ib,it)+ettc/kcb(ib) !evapotranspiracao total
	    etpb(ib,it)=etpb(ib,it)+etpc/kcb(ib) !transpiracao potencial
	    sssub(ib,it)=sssub(ib,it)+sssuc/kcb(ib) !lamina armazenada media camada superior
	    ssubb(ib,it)=ssubb(ib,it)+ssubc/kcb(ib) !lamina armazenada media camada inferior
	    sradb(ib,it)=sradb(ib,it)+sradc/kcb(ib) !lamina armazenada media camada radicular
	    asatb(ib,it)=asatb(ib,it)+100.*asatc/kcb(ib) !area saturada media em %
	    dsupb(ib,it)=dsupb(ib,it)+qsup(ic)*(dtp/1000.)/(acel(ic)*kcb(ib)) !escomento superficial da bacia em mm no deltat
	    dssub(ib,it)=dssub(ib,it)+qssu(ic)*(dtp/1000.)/(acel(ic)*kcb(ib)) !escomento intermediario da bacia em mm
	    dsubb(ib,it)=dsubb(ib,it)+qsub(ic)*(dtp/1000.)/(acel(ic)*kcb(ib)) !escomento subterraneo da bacia	em mm  
	endif
enddo !fim do loop das células

return
end
	
