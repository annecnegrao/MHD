!---------------------------------------------------------------------------------------------
!subrotina que analisa a qualidade do ajuste entre os hidrogramas diarios calculados e 
!observados com base em algumas funcoes objetivo.
!---------------------------------------------------------------------------------------------
subroutine fobj
use vars_main, only: qobs,qr,icalib, nbvfo,it,iniaju,nt,preb,ettb,dtp,adren,fob,dir_dados, &
                   & nb,filaju,iqobs,nobs,qcons,nash,lnash,cdatag,ncdom,icdom,kcb,ibac
implicit none
! variaveis auxiliares
real somares,somalres,somax2,somay2,somaxy,somalx2,somaly2,somalxy,xx,yy 
real somaobs !somatorio dos valores de vazao observados
real somacal !somatorio dos valores de vazao calculados
real somalobs !somatorio dos logaritmos de vazao observados
real somalcal !somatorio dos logaritmos de vazao calculados
real somapreb !somatorio da precipitacao da bacia
real somaecal !somatorio evaporacao total da bacia
real somaeobs !somatorio evaporacaco observada 
integer ib, itvalobs !conta o numero de dias em que existe vazao observada para comparar
real xmobs,xmlobs !vazao media observada e vazao media dos logaritmos das observacoes
real xmcal,xmlcal !vazao media calculada e vazao media dos logaritmos das vazoes calculadas
real r2(nb),errv(nb),erre(nb),r2l(nb) !erro nos volumes
real etcal(nb),etobs(nb) !evapotranspiracao calculada e observada
integer ibacaux,ibflag,ic,ic2,itvalsim

nash=0.0
r2=0.0
r2l=0.0
errv=99.0
erre=99.0
etcal=0.0
etobs=0.0

do ib=1,nobs
    if(icalib == 1 .and. ib /= nbvfo) cycle ! se estiver calibrando, a funcao objetiva calculada eh so do posto 
	somaobs=0.0
	somacal=0.0
	somalobs=0.0
	somalcal=0.0
	itvalobs=0
	itvalsim=0
	somax2=0.
	somay2=0.
	somaxy=0.
	somalx2=0.
	somaly2=0.
	somalxy=0.
	somares=0.
	somalres=0.
	somapreb=0.
    somaecal=0.
    	
	if(iniaju <= nt)then
		do it=iniaju,nt
			!calcula vazao media
			if(qobs(ib,it) >= 0.0) then
				somaobs=somaobs+qobs(ib,it)
				somacal=somacal+qr(ib,it)
				xx=qobs(ib,it)-qr(ib,it)
				somares= somares + xx*xx ! soma residuos
				!logaritmos
				xx=log(qobs(ib,it)+0.0000001)
				yy=log(qr(ib,it)+0.0000001)
				somalobs=somalobs+xx !soma um valor muito pequeno para evitar logaritmo de zero
				somalcal=somalcal+yy !soma um valor muito pequeno para evitar logaritmo de zero
				somalres=somalres+(xx-yy)*(xx-yy) ! soma residuos dos logaritmos								
				itvalobs=itvalobs+1
			endif
			
			!calcula prec e ettsim
			ibacaux=ib
			call dominiobacia(ibacaux)
			
			ibflag=0            
            !loop sobre as celulas do dominio da bacia
            do ic2=1, ncdom    
                ic = icdom(ic2)                    
                
                if(ibflag/=ibac(ic))then                    
                    ibflag=ibac(ic) 
                    somapreb=somapreb+preb(ibflag,it)*kcb(ibflag)                                            
                endif                
            enddo              
			somaecal=somaecal+ettb(ib,it) ! soma evaporacao total calculada da bacia
			itvalsim=itvalsim+1
			
		enddo
        if(itvalobs > 25) then
			xmobs=somaobs/itvalobs ! media das vazoes observacoes
			xmlobs=somalobs/itvalobs ! media do logaritmo das observacoes
			xmcal=somacal/itvalobs ! media das vazoes calculadas
			xmlcal=somalcal/itvalobs ! media dos logaritmos das vazoes calculadas
			do it=iniaju,nt
				if(qobs(ib,it) >= 0.0) then
				    yy=qr(ib,it)-xmcal
				    xx=qobs(ib,it)-xmobs
					somaxy=somaxy + xx*yy
				    somax2= somax2 + xx*xx
				    somay2= somay2 + yy*yy
				    yy=log(qr(ib,it)+0.0000001)-xmlcal
				    xx=log(qobs(ib,it)+0.0000001)-xmlobs
					somalx2= somalx2 + xx*xx
				    somaly2= somaly2 + yy*yy
				    somalxy= somalxy + xx*yy
				endif					
			enddo
			somapreb = (somapreb/real(ncdom))/itvalsim
			
			somaeobs=somapreb-((somaobs+qcons(ib)/itvalobs)*(dtp/1000.)/adren(iqobs(ib))/itvalobs) ! soma evaporacao media observada
			
			!coeficiente de nash
			nash(ib)=1.-somares/somax2
			!coeficiente de nash dos logartimos das vazoes
			lnash(ib)=1.-somalres/somalx2
			!erro de volume
			errv(ib)=(somacal-somaobs)/somaobs
			!coeficiente de correlacao: r2
			r2(ib)=somaxy*somaxy/(somax2*somay2)
			!coeficiente de correlacao dos logaritmos: r2log
			r2l(ib)=somalxy*somalxy/(somalx2*somaly2)
			
			!erro na evaporacao, apenas quando esta simulando			
			etcal(ib)=somaecal/itvalsim*(86400/dtp)
			etobs(ib)=somaeobs 
			erre(ib)=(etcal(ib)-etobs(ib))/etobs(ib)			
		endif
	endif
enddo	

if(icalib == 1) then ! esta calibrando, calcula o valor da funcao objetivo  
    fob = 10000*( (1-nash(nbvfo)) + (1-lnash(nbvfo)))   
    return
else
    write(*,*)	
    write(*,*)'estatisticas:'
    open(filaju,file=dir_dados//'saida/estatisticas.hig',status='unknown')
    write(filaju,*)"intervalo de analise: ",cdatag(iniaju)(1:2),'/',cdatag(iniaju)(3:4),'/',cdatag(iniaju)(5:8), &
                  & "-",cdatag(nt)(1:2),'/',cdatag(nt)(3:4),'/',cdatag(nt)(5:8)
    !10 format(a,i2.2)
    write(*,'(a8,8a8)')'bacia','nash','lnash','r2','r2l','errv','erre','etcal','etobs'
    write(filaju,'(a8,8a8)')'bacia','nash','lnash','r2','r2l','errv','erre','etcal','etobs'
    do ib=1,nobs 
        write(*,'(i8,8f8.3)')ib,nash(ib),lnash(ib),r2(ib),r2l(ib),errv(ib),erre(ib),etcal(ib),etobs(ib)        
        write(filaju,'(i8,8f8.3)')ib,nash(ib),lnash(ib),r2(ib),r2l(ib),errv(ib),erre(ib),etcal(ib),etobs(ib)
    enddo
    write(*,*)
    close(filaju)
endif

return
end