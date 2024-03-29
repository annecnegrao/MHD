#!/bin/bash
gfortran -c Vars_Calib.f90
gfortran -c Vars_main.f90
gfortran -c MHD.f90
gfortran -c Alloca_Vars.f90
gfortran -c Caldat.f90
gfortran -c Dominio.f90 
gfortran -c InterpMusk.f90
gfortran -c Atualiza.f90
gfortran -c LeCell.f90 
gfortran -c LeFix.f90
gfortran -c Calibra.f90
gfortran -c Celula.f90
gfortran -c Escoamentos.f90
gfortran -c Evaporacao.f90
gfortran -c Floodplain.f90
gfortran -c FObj.f90 
gfortran -c Julday.f90
gfortran -c LeMet.f90 
gfortran -c LeMetPrev.f90 
gfortran -c LeQobs.f90 
gfortran -c LeSolo.f90 
gfortran -c LeSubst.f90 
gfortran -c LeUso.f90 
gfortran -c LeVeg.f90 
gfortran -c Modelo.f90
gfortran -c Musk.f90 
gfortran -c Musk_NL.f90 
gfortran -c Parcel.f90 
gfortran -c Parcunge.f90 
gfortran -c Previsao.f90 
gfortran -c QQMet.f90 
gfortran -c Rede.f90 
gfortran -c RedeIni.f90 
gfortran -c Simula.f90 
gfortran -c Sort.f90 
gfortran -c Transpiracao.f90
gfortran -c BalancoHidrico.f90
gfortran -o mhd.exe *.o
rm -f *.o
rm -f *.mod
