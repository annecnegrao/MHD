!---------------------------------------------------------------------------------------------
! Subrotima para calcular a evaporacao do dossel de acordo com Van Dijk, 2001 - Modelo de Gash
!---------------------------------------------------------------------------------------------
subroutine evaporacao(ic,iu)
use vars_main, only: jdia,ta,tdew,patm,vv,pre,ec,zveg,dveg,z0,iaf,cover,scmax,xcel,ycel,gamma, &
                   & delta,imes,ra,nu,pi,dtp,roc,alb,tontem,e0,e0ag,ihora
implicit none
real,parameter:: stefan=4.903e-9 !constante de stefan bolstzmann em mj m-2 k-4 d-1
integer iu,ic	
real clate,es,ed,rhocp,pl,epr,pint,rol,rn,g,sn,lonmc
real declin,ws,lat,rae,roc0,dr,insol,aux,z,vz,w1,w2,w,sc,t,dth
real rnag,raag

!variaveis locais

! epr	: parametro - razao entre evaporacao umida e intensidade media (0.039)
! clate	: calor latente de vaporizacao mj m-2 dia-1, !calor latente de vaporizacao - em mj/k, temperatura em c
! es	: pressao de saturacao kpa !pressao de saturacao do vapor da agua em kpa
! ed	: pressao de vapor kpa !pressao aual de vapor real em kpa
! rhocp : densidade do ar a pressao constante em kg/m3 vezes o calor especifico do ar umido mj/kg/c 
! mspa	: densidade do ar umido
! z0m   : rugosidade s m-1 ! rugosidade para transferencia de momentum, resistencia aerodinamica de maidment e chou
! dveg  : plano de deslocamento zero
! mesp	: massa especifica da agua
! pl	: capacidade de armazenamento actual da vegetacao mm ! chuva necessaria para saturar o dossel
! pint	: interceptacao no dossel mm
! rol   : radiacao liquida de onda longa
! rnag  : radiacao liquida na agua
! raag  : resistência aerodinвmica da agua
! rn    : radiacao liquida (mj/m2/dia)
! g     : fluxo de calor do solo	
! sn    : radiacao liquida de ondas curtas
! dr    : distвncia terra-sol
! lonmc : longitude no centro do fuso horario local (graus a oeste de Greenwhich)

es=0.6108*exp((17.27*ta(ic))/(237.3+ta(ic))) ! pressao de vapor de saturacao
ed=0.6108*exp((17.27*tdew(ic))/(237.3+tdew(ic))) ! pressao de vapor 
clate=(2.501-0.002361*ta(ic)) ! calor latente mj/c
gamma=0.0016286*patm(ic)/clate ! constante psicrometrica
delta=(4098.0*es)/((237.3+ta(ic))**2.) ! declividade da curva de pressao de vapor de saturacao
rhocp=0.00353128*patm(ic)/(clate*(275.+ta(ic))) ! rho	: massa especifica do ar
if(zveg(iu,imes) > 10.) then ! vegetacao mais alta que o anemometro
	! leva a velocidade do vento a altura do dossel assumindo grama com rugosidade = 0.078 m e deslocamento = 2.65 m 
    vz=vv(ic)*log(12.82*zveg(iu,imes)-3.4)/4.85+0.01
    z=zveg(iu,imes)
else
    vz=vv(ic)+0.01
    z=10.
endif
!aux=log((z-dveg(iu,imes))/z0(iu,imes))
!ra=aux*aux/(0.41*0.41*vz)
aux=log((z-dveg(iu,imes))/z0(iu,imes))*log((z-dveg(iu,imes))/(z0(iu,imes))) 
ra=aux/(0.41*0.41*vz) ! resistencia aerodinamica
aux=log((z-dveg(nu+1,imes))/z0(nu+1,imes))*log((z-dveg(nu+1,imes))/(z0(nu+1,imes)))
raag=aux/(0.41*0.41*vz) ! resistencia aerodinamica da agua
! calcula radiacao de onda longa
dr=1+0.033*cos(2*pi*jdia/365) ! dr	: distancia terra sol
declin=0.4093*sin(2*pi*jdia/365-1.405) ! declinacao
lat=ycel(ic)*pi/180. ! latitude em radianes
ws=acos(-tan(lat)*tan(declin)) 	! angulo solar
if(dtp >= 86400.) then	! delta t maior ou igual a 1 dia
    insol=24*ws/pi ! horas de insolacao
    rae=37.586*dr*(ws*sin(lat)*sin(declin)+cos(lat)*cos(declin)*sin(ws)) ! radiacao extraterrestre mj/m2/dia 
    roc0=max(0.75*rae,roc(ic)) ! radiacao de onda curta incidente em dia de ceu claro
    rol=stefan*(ta(ic)+273.16)**4*(0.34-0.14*sqrt(ed))*(1.35*roc(ic)/roc0-0.35) ! radiacao liquida de onda longa
    sn=roc(ic)*(1.-alb(iu,imes)) ! radicao liquida de onda curta
    rn=sn-rol ! radiacao liquida
    sn=roc(ic)*(1.-alb(nu+1,imes)) ! radicao liquida de onda curta da agua
    rnag=sn-rol ! radiacao liquida da agua
    g=0.38*(ta(ic)-tontem(ic)) ! fluxo de calor do solo do handbook of hydrology
	e0=(delta*(rn+g)+86400.*rhocp*(es-ed)/ra)/(delta+gamma)*dtp/(86400.*clate) !mm no deltat   	
    e0ag=(delta*rnag+86400.*rhocp*(es-ed)/raag)/(delta+gamma)*dtp/(86400.*clate) !mm no deltat    
else ! delta t menor a um dia
    aux = 2*pi*(jdia-81)/364.
    sc = 0.1645*sin(2*aux)-0.1255*cos(aux)-0.025*sin(aux) ! correcao sazonal para hora solar (horas)
    aux=xcel(ic)/15.
    lonmc=int(aux)
    if((aux-lonmc) > 0.5 ) lonmc=lonmc-1
    lonmc=lonmc*15. 
    dth=dtp/3600. ! intervalo de tempo em horas
    t=float(ihora)-0.5*dth ! standard clock time at the midpoint of the period (hs)
    if(t < 0.) t = t + 24.
    w=(pi/12.)*(t+0.06667*(lonmc-xcel(ic))+sc-12.) ! solar time angle at midpoint of hourly or shorter period (rad)
    w1=w-pi*dth/48. ! dividi por 48, era 24
    if(w1 < -ws) then
        w1=-ws
    elseif(w1 > ws) then
        w1=ws
    endif
    w2=w+pi*dth/48. ! dividi por 48, era 24
    if(w2 < -ws) then
        w2=-ws
    elseif(w2 > ws) then
        w2=ws
    endif
	rae=18.793*dr*((w2-w1)*sin(lat)*sin(declin)+cos(lat)*cos(declin)*(sin(w2)-sin(w1))) ! radiacao extraterrestre em mj/m2/hora
	if(rae <= 0.) then
	    roc0=0.5 ! razao assumida para horario noturno (fa0 56) para condicoes umidas, 0.8 para condicoes aridas
	else
        roc0=roc(ic)/max(0.75*rae,roc(ic)) ! relative solar radiation
    endif
    rol=(stefan/24.)*(ta(ic)+273.16)**4*(0.34-0.14*sqrt(ed))*(1.35*roc0-0.35) ! radiacao liquida de onda longa  em mj/m2/hora
    sn=roc(ic)*(1.-alb(iu,imes)) ! radiacao liquida de onda curta
    rn=sn-rol ! radiacao liquida
	if(rae <= 0.) then
	    g = -2.1 * exp(-0.5 * iaf(iu,imes)) * rn ! choudhury (1989) densidade de fluxo de calor no solo (g) sob condicoes noturnas
	else
        g = -0.4 * exp(-0.5 * iaf(iu,imes)) * rn ! choudhury (1989) densidade de fluxo de calor no solo (g) sob condicoes diurnas
    endif
    sn=roc(ic)*(1.-alb(nu+1,imes)) ! radicao liquida de onda curta da agua
    rnag=sn-rol ! radiacao liquida da agua
	e0=(delta*(rn-g)+3600.*rhocp*(es-ed)/ra)/(delta+gamma)*dtp/(3600.*clate) !mm no deltat
    e0ag=(delta*rnag+3600.*rhocp*(es-ed)/raag)/(delta+gamma)*dtp/(3600.*clate) !mm no deltat  
endif
if (scmax(iu) > 0.0) then !exclui o caso de agua livre
	epr =0.031 !razao evaporacao umida/intensidade media (parametro)- germer et al 2006
	pl=-scmax(iu)*cover(iu,imes)/epr*log(1-epr) ! capacidade de saturacao do dossel
    ! calcula a interceptacao do dossel [mm/dia]
	if (pre(ic).le.pl) then
		pint = cover(iu,imes)*pre(ic)
	else
		pint = cover(iu,imes)*(pl + epr*(pre(ic)-pl))
	endif	
	! retira evaporacao do reservatorio de interceptacao
	if(pint.gt.e0) pint=e0
	pint=min(pint,pre(ic))
	ec(ic,iu) = pint
endif
return
end