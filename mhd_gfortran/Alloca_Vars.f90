!---------------------------------------------------------------------------------------------
!subrotina de alocacao de memoria das variaveis principais
!---------------------------------------------------------------------------------------------
subroutine alloca_vars(iop)
use vars_main !,only: nc,nb,nu,nt,icodmusk,evq,qrioini,qcontorm, &
!                  & qcontorj,prebanual,kbanual,qref,puso,lambda,lambdam,lambdan,ac,lambdai, &
!                  & tanb,acel,adren,srio,decl,icell,icdom,idren,xcel,ycel,ibac,hcel,lcel, &
!                  & celjus,prec,brio,iexut,neta,d1,d2,d3,smax,kss,tsub,mu,alpha,csi,ssmax, &
!                  & srmax,tcon,cb,cs,ntrmax,ssub,srad,secrbs,secrjal,secrjah,ntrc,dt,cel,tks, &
!                  & tkb,alb,iaf,zveg,rc
implicit none
integer iop
save

allocv_case: select case (iop) !verifica se alloca ou dealloca
case (0) ! alloca
	allocate (icodmusk(nc)) !codigo que indica linear ou nao linear
	allocate (evq(nc)) !evaporacao direta da superficie liquida da celula em m3/s
	allocate (qrioini(nc,ntrmax+1)) !condicao incial da propagacao muskingum cunge
	allocate (qcontorm(nc+1,25),qcontorj(nc+1,25)) !vazao da condicao de contorno de muskingum cunge
	allocate (qref(nc)) !vazao de referencia
 	allocate (puso(nc,nu)) !proporcao de usos na celula
 	allocate (lambda(nc,nu),lambdam(nc,nu),lambdan(nc,nu),ac(nc,50),lambdai(nc,50),tanb(nc)) !indices topograficos e histograma na celula
	allocate (acel(nc),adren(nc),srio(nc),decl(nc),icell(nc),icdom(nc)) !area da celula, area drenada, comprimento e declividade do rio, no. celula, no. celula do dominio
	allocate (idren(nc)) !vetor que ordena os valores de area de drenagem em ordem crescente para a rotina rede
	allocate (xcel(nc),ycel(nc)) !coordenadas do centro da celula
	allocate (ibac(nc),hcel(nc),lcel(nc),celjus(nc)) !bacia,hmax,hmin,celula de jusante
 	allocate (prec(nc)) !chuva media na celula,
	allocate (brio(nc)) !largura do rio
	allocate (iexut(nb)) !indica celula do exutorio da bacia
	allocate (neta(nu),d1(nb),d2(nb),d3(nb)) ! parametro neta de Brooks Corey e profundidades do solo
	allocate (smax(nb,nu),kss(nb,nu),tsub(nb,nu),mu(nb,nu),alpha(nb,nu),csi(nb,nu),ssmax(nb,nu),srmax(nb,nu)) ! parametros do solo
	allocate (tcon(nc),cb(nb),cs(nb)) !tempo de concentracao e parametros da propagacao na celula
	allocate (ssub(nc,nu)) !armazenamento na camada inferior do solo
	allocate (srad(nc,nu)) !armazenamento na zona radicular do solo
	allocate (secrbs(nu)) !parametro de transpiracao de solo nu
	allocate (secrjal(nu),secrjah(nu)) !parametros de estresse de Jarvis
	allocate (ntrc(nc)) !numero de subtrechos de rio para Muskingum Cunge (diferente para cada celula)
    allocate (dt(nc)) !intervalo de tempo Muskingum Cunge (diferente para cada celula)
	allocate (cel(nc),tks(nc),tkb(nc)) !celeridade no rio e coeficientes do rls da celula
	! nos parametros de vegetacao, o uso nu+1 eh reservado a area alagavel
	allocate (alb(nu+1,12),iaf(nu+1,12),zveg(nu+1,12),rc(nu+1)) ! albedo, indice de area foliar, altura da veg e resistencia minima
	allocate (dveg(nu+1,12),z0(nu+1,12)) ! plano de deslocamento zero e rugosidade			
	allocate (cover(nu+1,12),scmax(nu+1))
	allocate (qesp(nb),qcons(nb)) !vazao especifica de base (m3/s/km2) e consumo de agua no exutorio da bacia (m3/seg)
	allocate (qbsub(nt+ntprev),qbssu(nt+ntprev),qbsup(nt+ntprev)) !vazao de acordo com origem
	allocate (dsupb(nb,nt+ntprev),dssub(nb,nt+ntprev),dsubb(nb,nt+ntprev)) !medias bacia
	allocate (ecb(nb,nt+ntprev),etpb(nb,nt+ntprev),ettb(nb,nt+ntprev),sssub(nb,nt+ntprev), &
             & asatb(nb,nt+ntprev),preb(nb,nt+ntprev),ssubb(nb,nt+ntprev),sradb(nb,nt+ntprev)) !medias bacia
	allocate (cdatag(nt+ntprev)) ! vetor com datas em formato character
	allocate (kcb(nb)) !numero de celulas em cada sub-bacia
	allocate (qr(nobs,nt),qobs(nobs,nt)) !vazao calculada e observada nos exutorios das bacias
	allocate (qlido(numsubst,nt)) !vazao que eh lida para substituir calculada
	allocate (qrg(nhidg,nt)) !hidrogramas para gravacao		
	allocate (qm1(nc+1),qj1(nc),qm2(nc+1),qj2(nc)) !vazoes a montante e a jusante em i
	allocate (qcel1(nc),qcel2(nc)) !vazoes originadas na celula nos instantes t e t+1 na celula i
	allocate (pmsub(nc+1),pmssu(nc+1),pmsup(nc+1),pjsub(nc+1),pjssu(nc+1),pjsup(nc+1)) !proporcoes de origem das vazoes no rio
	allocate (vrsub(nc),vrssu(nc),vrsup(nc)) !volumes das proporcoes no trecho
	allocate (ett(nc,nu)) !evapotranspiracao total
	allocate (etp(nc,nu)) !transpiracao potencial
	allocate (ec(nc,nu)) ! evaporacao no dossel
	allocate (asat(nc,nu)) !armazenamento dossel
	allocate (sssu(nc,nu)) !armazenamento na camad superior do solo
	allocate (qsub(nc),qssu(nc),qsup(nc)) !vazao na celula
	allocate (vsub(nc),vssu(nc),vsup(nc)) !volume na celula
	allocate (tontem(nc)) !temperatura no dia anterior
	allocate (pre(nc)) !chuva no intervalo na celula
	allocate (ta(nc),tdew(nc),vv(nc),roc(nc),patm(nc)) !temperatura, umidade, vento, insol., pressao
	allocate (preall(nc,nt),taall(nc,nt),tdewall(nc,nt),vvall(nc,nt),rocall(nc,nt),patmall(nc,nt))! para todo o periodo de simulacao
	allocate (prad(nu+1,12)) !profundidade radicular
	allocate (nash(nb),lnash(nb)) ! coef. Nash e Log-Nash

	
case (1) ! dealloca
	deallocate (icodmusk) !codigo que indica linear ou nao linear
	deallocate (evq) !evaporacao direta da superficie liquida da celula em m3/s
	deallocate (qrioini) !condicao incial da propagacao Muskingum Cunge
	deallocate (qcontorm,qcontorj) !vazao da condicao de contorno de Muskingum Cunge
	deallocate (qref) !vazao de referencia
 	deallocate (puso)!proporcao de usos na celula
 	deallocate (lambda,lambdam,lambdan,ac,lambdai,tanb) !indices topograficos e histograma na celula
 	deallocate (acel,adren,srio,decl,icell,icdom) !area da celula, area drenada, comprimento e declividade do rio, no. celula, no. celula do dominio
 	deallocate (idren) !vetor que ordena os valores de area de drenagem em ordem crescente para a rotina rede
	deallocate (xcel,ycel) !coordenadas do centro da celula
	deallocate (ibac,hcel,lcel,celjus) !bacia,hmax,hmin,celula de jusante
 	deallocate (prec) !chuva media na celula
	deallocate (brio) !largura do rio
	deallocate (iexut) !indica celula do exutorio da bacia
	deallocate (neta,d1,d2,d3) ! parametro neta de Brooks Corey e profundidades do solo
	deallocate (smax,kss,tsub,mu,alpha,csi,ssmax,srmax) !parametros de solo da celula
    deallocate (tcon,cb,cs) !tempo de concentracao e parametros da propagacao na celula
    deallocate (ssub) !armazenamento na camada inferior do solo
	deallocate (srad) !armazenamento na zona radicular do solo
	deallocate (secrbs) !parametro de transpiracao de solo nu
	deallocate (secrjal,secrjah) !parametros de estresse de Jarvis
	deallocate (ntrc,dt)
	deallocate (cel,tks,tkb)
	deallocate (alb,iaf,zveg,rc)
	deallocate (dveg,z0) ! plano de deslocamento zero e rugosidade	
	deallocate (qesp,qcons) !vazao especifica de base (m3/s/km2) e consumo de agua noe xutorio da bacia (m3/seg)
	deallocate (qbsub,qbssu,qbsup) !vazao de acordo com origem
	deallocate (dsupb,dssub,dsubb) !medias bacia
	deallocate (ecb,etpb,ettb,sssub,preb,ssubb,asatb,sradb) !medias bacia
	deallocate (cdatag)
	deallocate (kcb) !numero de celulas em cada sub-bacia
	deallocate (qr,qobs)!vazao calculada e observada nos exutorios das bacias
	deallocate (qrg) !vazao nos exutorios das bacias
	deallocate (qm1,qj1,qm2,qj2) !vazoes a montante e a jusante em i
	deallocate (pmsub,pmssu,pmsup,pjsub,pjssu,pjsup) !proporcoes de origem das vazoes no rio
	deallocate (vrsub,vrssu,vrsup) !volumes das proporcoes no trecho
	deallocate (ett) !evapotraspiracao total
	deallocate (etp) !transpiracao potencial
	deallocate (ec)	!evaporacao no dossel
	deallocate (asat) !armazenamento dossel
	deallocate (sssu) !armazenamento na camada superior do solo
	deallocate (qsub,qssu,qsup) !vazao na celula
	deallocate (vsub,vssu,vsup) !volume na celula
	deallocate (tontem) !temperatura no dia anterior
	deallocate (pre) !chuva no intervalo na celula
	deallocate (qcel1,qcel2) !vazoes originadas na celula nos instantes t e t+1 na celula i
	deallocate (qlido) !vazao de substituiзao da calculada
	deallocate (ta,tdew,vv,roc,patm) !temperatura, ponto de orvalho, vento, insol., pressao
	deallocate (preall,taall,tdewall,vvall,rocall,patmall)
	deallocate (prad) !profundidade radicular
case default
	stop ' erro: iop desconhecido no rotina alloca_calib!!!'
end select allocv_case

   
end
