!---------------------------------------------------------------------------------------------
!subrotina que propaga na planicie
!---------------------------------------------------------------------------------------------
subroutine floodplain
use vars_main, only: nfp,dtp,dtfp,qrfp,icellfp,ntrc,qm1,qm2,qj1,qj2,rugman,brio,decl,hcel,lcel, &
& zvert,zfp,cvl,bvl,cva,bva,srio,celln,afp,vfp,zfpt,sssu,zfpb,fch,cellw,vfpt,afpt,ibac,nu,puso, &
& acel,ssmax
implicit none
integer ic,iu,ib,ifp,k
real qn(nfp),qw(nfp),qs(nfp),qe(nfp) ! vazao de troca entre celulas de planicie nas direcoes cardinais
real zr(nfp) ! cota do rio
real ppu,duso,aux,area,hf
real qmx,qjx,c
integer ndtfp,idtfp

ndtfp=dtp/dtfp ! numero de subintervalos da rotina floodplain em seg

qrfp=0.
do idtfp=1,ndtfp
    qn=0.
    qw=0.
    qs=0.
    qe=0.
    ! determina a altura da água em cada celula do rio na planicie
    do ifp=1,nfp
        ic=icellfp(ifp)
        if(ntrc(ic) == 0)cycle !celula sem rio
        
        !interpola as vazoes de montante e jusante entre t e t+1
        qmx=qm1(ic)+idtfp*(qm2(ic)-qm1(ic))/ndtfp
        qjx=qj1(ic)+idtfp*(qj2(ic)-qj1(ic))/ndtfp        
        
        ! estima a cota media de agua do rio, zr, usando manning, para secao retangular   
        zr(ifp)=0.5*(qjx**0.6+qmx**0.6)*(rugman/brio(ic))**0.6/decl(ic)**0.3
        zr(ifp)=zr(ifp)+0.5*(hcel(ic)+lcel(ic))
        
        if(zr(ifp) <= zvert(ifp)) then ! o nivel do rio esta abaixo da cota de vertimento
            if(zfp(ifp) > zvert(ifp)) then !planicie verte no rio como vertedor livre
                qrfp(ifp)=cvl*bvl*(zfp(ifp)-zvert(ifp))**1.5
            endif
        else ! rio acima da cota de vertimento
            if(zfp(ifp) > zvert(ifp)) then ! vertedor afogado, o sinal de aux define o sentido
                ! aux > 0 planicie -> rio; aux < 0 rio -> planicie
                aux=zfp(ifp)-zr(ifp)
                qrfp(ifp)=cva*bva*aux*abs(aux)**0.5
            else  ! vertedor livre, o rio verte na planicie
                aux =zr(ifp)-zvert(ifp)
                qrfp(ifp)=-cvl*bvl*(aux)**1.5
            endif
            ! amortece a vazao que sai do rio para evitar esvaziamentos bruscos
            if(qrfp(ifp) < 0.) then ! vazao do rio para a planicie
                aux = abs(aux)*brio(ic)*srio(ic) ! volume do rio acima da cota de vertimento ou da cota na planicie
                c=aux/(qrfp(ifp)*dtfp) ! coeficiente de reducao
                if(c < 1. ) then ! verteu mais do que devia
                    qrfp(ifp)=qrfp(ifp)*c
                endif
            endif
        endif               
        
    enddo ! fim das trocas rio planicie
    ! calculo das trocas entre elementos da planicie
    qn=0.
    qw=0.
    qs=0.
    qe=0.
    do ifp=1,nfp
        ! calcula primeiramente ao norte
        if(celln(ifp) > 0) then ! existe celula de planicie ao norte
            hf=max(zfp(celln(ifp)),zfp(ifp))-max(zfpb(celln(ifp)),zfpb(ifp))
            if(hf > 0) then
                aux=zfp(celln(ifp))-zfp(ifp)
                qn(ifp)=fch*hf**1.67*sqrt(abs(aux)) ! vazao na direcao norte
                if(aux < 0) qn(ifp)=-qn(ifp) ! a celula perde agua
                qs(celln(ifp))=-qn(ifp) ! equivale a - a vazao na direcao da celula ao norte
            endif
        endif
        ! calcula ao oeste
        if(cellw(ifp) > 0) then ! existe celula de planicie ao oeste
            hf=max(zfp(cellw(ifp)),zfp(ifp))-max(zfpb(cellw(ifp)),zfpb(ifp))
            if(hf > 0) then
                aux=zfp(cellw(ifp))-zfp(ifp)
                qw(ifp)=fch*hf**1.67*sqrt(abs(aux)) ! vazao na direcao oeste
                if(aux < 0) qw(ifp)=-qw(ifp) ! a celula perde agua
                qe(cellw(ifp))=-qw(ifp) ! equivale a - (menos) a vazao na direcao da celula ao oeste
            endif
        endif   
    enddo
    ! verifica se as celulas de planicie nao esvaziaram alem do que deviam
    do ifp =1,nfp
        aux=qn(ifp)+qw(ifp)+qs(ifp)+qe(ifp)-qrfp(ifp)
        if(aux < 0.) then ! a celula perde agua
            c = -(zfp(ifp)-zfpb(ifp))*afp(ifp)/(aux*dtfp) ! coeficiente linear de esvaziamento
            if(c < 1.) then ! esvaziou alem do volume armazenado na celula
                qn(ifp)=qn(ifp)*c
                qw(ifp)=qw(ifp)*c
                qe(ifp)=qe(ifp)*c
                qs(ifp)=qs(ifp)*c
                qrfp(ifp)=qrfp(ifp)*c
            endif
            ! compatibiliza as vazoes
            if(celln(ifp) > 0) then ! existe celula de planicie ao norte
                if(qn(ifp) < 0. ) then ! se for menor que zero escolhe a maior vazao
                    qn(ifp)=max(qn(ifp),-qs(celln(ifp)))
                else ! se for maior que zero escolhe a vazao minima
                    qn(ifp)=min(qn(ifp),-qs(celln(ifp)))
                endif
                qs(celln(ifp))=-qn(ifp)
            endif
            if(cellw(ifp) > 0) then ! existe celula de planicie ao oeste
                if(qw(ifp) < 0.) then
                    qw(ifp)=max(qw(ifp),-qe(cellw(ifp)))
                else
                    qw(ifp)=min(qw(ifp),-qe(cellw(ifp)))
                endif
                qe(cellw(ifp))=-qw(ifp)
            endif
        endif
    enddo
    ! balanco hidrico de cada elemento de celula
    do ifp=1,nfp
        vfp(ifp)=vfp(ifp)+(qn(ifp)+qw(ifp)+qs(ifp)+qe(ifp)-qrfp(ifp))*dtfp ! vol armazenado no floodplain
        ! estimativa preliminar de area alagada
        if(vfp(ifp) <= vfpt(ifp,1)) then
            area=0.
        else
            area=100.
            do k=2,50
                if(vfp(ifp) < vfpt(ifp,k)) then
                    area=afpt(ifp,k-1)+(vfp(ifp)-vfpt(ifp,k-1))*(afpt(ifp,k)-afpt(ifp,k-1))/(vfpt(ifp,k)-vfpt(ifp,k-1))
                    exit
                endif
            enddo
        endif
        ic=icellfp(ifp)
        ib=ibac(ic)
        ! correcoes na celula em decorrencia da mudanca de area alagada
        do iu=1,nu
            ppu=puso(ic,iu)/100.*(1.-afp(ifp)/100.) ! percentagem de uso iu no tempo t corrigido pela area alagada
            duso=puso(ic,iu)/100.*(afp(ifp)-area)/100. ! expansao/contracao do uso iu no tempo t+1 em funcao da mudanca de area alagada
            if(duso > 0) then ! expansao da area alagada
                ! corrige o volume alagado em funcao da demanda de agua necessaria para saturar a primeira camada das nova areas alagadas
                vfp(ifp)=vfp(ifp)-duso*acel(ic)*(ssmax(ib,iu)-sssu(ic,iu))/1000.
            else ! contracao de area saturada
                ! corrige o armazenamento na camada superior do solo nas areas que agora nao estao alagadas                
                sssu(ic,iu)=(sssu(ic,iu)*ppu+ssmax(ib,iu)*duso)/(ppu+duso)
            endif
        enddo        
        if(vfp(ifp) < 0.) then
            ! erro de balanco eh alocado totalmente no fluxo planicie rio, pois elimina a necessidade de corrigir toda a planicie
            qrfp(ifp)=qrfp(ifp)+vfp(ifp)/86400.
            vfp(ifp)=0.
        endif
        ! area alagada eh reajustada e calcula-se a altura de alagamento 
        if(vfp(ifp) <= vfpt(ifp,1)) then
            afp(ifp)=0.
            zfp(ifp)=zfpb(ifp)
        else
            afp(ifp)=100.
            zfp(ifp)=zfpt(ifp,50)
            do k=2,50
                if(vfp(ifp) < vfpt(ifp,k-1)) then
                    afp(ifp)=afpt(ifp,k-1)+(vfp(ifp)-vfpt(ifp,k-1))*(afpt(ifp,k)-afpt(ifp,k-1))/(vfpt(ifp,k)-vfpt(ifp,k-1))
                    zfp(ifp)=zfpt(ifp,k-1)+(vfp(ifp)-vfpt(ifp,k-1))*(zfpt(ifp,k)-zfpt(ifp,k-1))/(vfpt(ifp,k)-vfpt(ifp,k-1))
                    exit
                endif
            enddo
        endif
    enddo ! fim do loop celulas de planicie
enddo ! fim do loop de subintervalos de tempo
return
end