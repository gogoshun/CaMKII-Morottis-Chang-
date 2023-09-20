using Catalyst
using DifferentialEquations
using Plots
using ModelingToolkit
using Sundials
using ODEInterfaceDiffEq


function get_Morotti_equations()
    @variables t
    D = Differential(t)

    # ECC variables
    @variables Na_m(t) Na_h(t) Na_j(t) ICa_HH4(t) ICa_HH5(t) ICa_HH6(t) ICa_HH7(t) Itos_x(t) Itos_y(t) Itof_x(t) Itof_y(t)
    @variables Ikr(t) IKs(t) RyR_R(t) RyR_O(t) RyR_I(t) NaBj(t) NaBsl(t) TnCL(t) TnCHc(t) TnCHm(t) CaM(t)
    @variables Myosin_ca(t) Myosin_mg(t) SRB(t) SLLj(t) SLLsl(t) SLHj(t) SLHsl(t) Csqn(t) Ca_sr(t) Naj(t) Nasl(t)
    @variables Nai(t) Ki(t) Ca_j(t) Ca_sl(t) Cai(t) Vm(t) Itos_r(t) influx_LTCC(t) influx_PMCA(t) influx_NCX(t) influx_ICa(t)
    @variables Na_late_h(t) CNa2(t) CNa1(t) ONa(t) IFNa(t) I1Na(t) CNa3(t) ICNa2(t) ICNa3(t) LONa(t)
    @variables LCNa1(t) LCNa2(t) LCNa3(t) C2_m1j(t) C1_m1j(t) I1Ca_m1j(t) I2Ca_m1j(t) I1Ba_m1j(t) I2Ba_m1j(t)
    @variables C2_m2j(t) C1_m2j(t) I1Ca_m2j(t) I2Ca_m2j(t) I1Ba_m2j(t) I2Ba_m2j(t)
    @variables C2_m1sl(t) C1_m1sl(t) I1Ca_m1sl(t) I2Ca_m1sl(t) I1Ba_m1sl(t) I2Ba_m1sl(t)
    @variables C2_m2sl(t) C1_m2sl(t) I1Ca_m2sl(t) I2Ca_m2sl(t) I1Ba_m2sl(t) I2Ba_m2sl(t)
    @variables IKs_x(t) IKs1_y(t) Iss(t) IKs2_y(t)

    # CaMDyad / Cyt / SL variables
    @variables CaM_dyad(t) Ca2CaM_dyad(t) Ca4CaM_dyad(t) CaMB_dyad(t) Ca2CaMB_dyad(t) Ca4CaMB_dyad(t) Pb2_dyad(t) Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t)
    @variables Pa_dyad(t) Ca4CaN_dyad(t) CaMCa4CaN_dyad(t) Ca2CaMCa4CaN_dyad(t) Ca4CaMCa4CaN_dyad(t)
    @variables CaM_sl(t) Ca2CaM_sl(t) Ca4CaM_sl(t) CaMB_sl(t) Ca2CaMB_sl(t) Ca4CaMB_sl(t) Pb2_sl(t) Pb_sl(t) Pt_sl(t) Pt2_sl(t)
    @variables Pa_sl(t) Ca4CaN_sl(t) CaMCa4CaN_sl(t) Ca2CaMCa4CaN_sl(t) Ca4CaMCa4CaN_sl(t)
    @variables CaM_cyt(t) Ca2CaM_cyt(t) Ca4CaM_cyt(t) CaMB_cyt(t) Ca2CaMB_cyt(t) Ca4CaMB_cyt(t) Pb2_cyt(t) Pb_cyt(t) Pt_cyt(t) Pt2_cyt(t)
    @variables Pa_cyt(t) Ca4CaN_cyt(t) CaMCa4CaN_cyt(t) Ca2CaMCa4CaN_cyt(t) Ca4CaMCa4CaN_cyt(t)

    # CaMKII variables
    @variables LCC_PKAp(t) LCC_CKdyadp(t) RyR2809p(t) RyR2815p(t) PLBT17p(t) LCC_CKslp(t)
    @variables Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t) Pa_dyad(t) I1p_PP1(t) Pb_sl(t) Pt_sl(t) Pt2_sl(t) Pa_sl(t)

    #BAR variables
    @variables LR(t) LRG(t) RG(t) b1AR_S464(t) b1AR_S301(t) GsaGTPtot(t) GsaGDP(t) Gsby(t) AC_GsaGTP(t) PDEp(t)
    @variables cAMPtot(t) RC_I(t) RCcAMP_I(t) RCcAMPcAMP_I(t) RcAMPcAMP_I(t) PKACI(t) PKACI_PKI(t) RC_II(t) RCcAMP_II(t) RCcAMPcAMP_II(t)
    @variables RcAMPcAMP_II(t) PKACII(t) PKACII_PKI(t) I1p_PP1(t) I1ptot(t) PLBp(t) PLMp(t) LCCap(t) LCCbp(t) RyRp(t)
    @variables TnIp(t) KS79(t) KS80(t) KSp(t) CFTRp(t) KURp(t)


## CamDyad/ CaMSL/ CamCyt
    ## Parameters
    Mg = 1      # [mM]
    K = 135     # [mM]
    Btot_dyad = 0
    CaMKIItot_dyad = 120         # [uM]
    CaNtot_dyad = 3e-3/8.293e-4  # [uM]
    PP1tot_dyad = 96.5           # [uM]

    ## Parameters for Cyt and SL
    Btot = 24.2     # [uM]
    CaMKIItot = 120*8.293e-4  # [uM]
    CaNtot = 3e-3             # [uM]
    PP1tot = 0.57             # [uM]

    ## Parameters
    # Ca/CaM parameters
    if Mg <= 1
        Kd02 = 0.0025*(1+K/0.94-Mg/0.012)*(1+K/8.1+Mg/0.022)  # [uM^2]
        Kd24 = 0.128*(1+K/0.64+Mg/0.0014)*(1+K/13.0-Mg/0.153) # [uM^2]
    else
        Kd02 = 0.0025*(1+K/0.94-1/0.012+(Mg-1)/0.060)*(1+K/8.1+1/0.022+(Mg-1)/0.068)   # [uM^2]
        Kd24 = 0.128*(1+K/0.64+1/0.0014+(Mg-1)/0.005)*(1+K/13.0-1/0.153+(Mg-1)/0.150)  # [uM^2]
    end
    k20 = 10               # [s^-1]      
    k02 = k20/Kd02         # [uM^-2 s^-1]
    k42 = 500              # [s^-1]      
    k24 = k42/Kd24         # [uM^-2 s^-1]
    
    # CaM buffering (B) parameters
    k0Boff = 0.0014        # [s^-1] 
    k0Bon = k0Boff/0.2   # [uM^-1 s^-1] kon = koff/Kd
    k2Boff = k0Boff/100    # [s^-1] 
    k2Bon = k0Bon          # [uM^-1 s^-1]
    k4Boff = k2Boff        # [s^-1]
    k4Bon = k0Bon          # [uM^-1 s^-1]

    # using thermodynamic constraints
    k20B = k20/100 # [s^-1] thermo constraint on loop 1
    k02B = k02     # [uM^-2 s^-1] 
    k42B = k42     # [s^-1] thermo constraint on loop 2
    k24B = k24     # [uM^-2 s^-1]
    
    # CaMKII parameters
    # Wi Wa Wt Wp
    kbi = 2.2      # [s^-1] (Ca4CaM dissocation from Wb)
    kib = kbi/33.5e-3 # [uM^-1 s^-1]
    kib2 = kib
    kb2i = kib2*5
    kb24 = k24
    kb42 = k42*33.5e-3/5
    kpp1 = 1.72     # [s^-1] (PP1-dep dephosphorylation rates)
    Kmpp1 = 11.5    # [uM]
    kta = kbi/1000  # [s^-1] (Ca4CaM dissociation from Wt)
    kat = kib       # [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
    kt42 = k42*33.5e-6/5
    kt24 = k24
    kat2 = kib
    kt2a = kib*5
    
    # CaN parameters
    kcanCaoff = 1              # [s^-1] 
    kcanCaon = kcanCaoff/0.5               # [uM^-1 s^-1] 
    kcanCaM4on = 46            # [uM^-1 s^-1]
    kcanCaM4off = 1.3e-3       # [s^-1]
    kcanCaM2on = kcanCaM4on
    kcanCaM2off = 2508*kcanCaM4off
    kcanCaM0on = kcanCaM4on
    kcanCaM0off = 165*kcanCaM2off
    k02can = k02
    k20can = k20/165
    k24can = k24
    k42can = k20/2508
 
    ## Dyad Fluxes
    
    Ca_Dyad = Ca_j * 1e3
    # CaM Reaction fluxes
    B_dyad = Btot_dyad - CaMB_dyad - Ca2CaMB_dyad - Ca4CaMB_dyad
    rcn02_dyad = k02*Ca_Dyad^2*CaM_dyad - k20*Ca2CaM_dyad 
    rcn24_dyad = k24*Ca_Dyad^2*Ca2CaM_dyad - k42*Ca4CaM_dyad 
    # CaM buffer fluxes
    rcn02B_dyad = k02B*Ca_Dyad^2*CaMB_dyad - k20B*Ca2CaMB_dyad 
    rcn24B_dyad = k24B*Ca_Dyad^2*Ca2CaMB_dyad - k42B*Ca4CaMB_dyad 
    rcn0B_dyad = k0Bon*CaM_dyad*B_dyad - k0Boff*CaMB_dyad 
    rcn2B_dyad = k2Bon*Ca2CaM_dyad*B_dyad - k2Boff*Ca2CaMB_dyad 
    rcn4B_dyad = k4Bon*Ca4CaM_dyad*B_dyad - k4Boff*Ca4CaMB_dyad 
    # CaN reaction fluxes 
    Ca2CaN_dyad = CaNtot_dyad - Ca4CaN_dyad - CaMCa4CaN_dyad - Ca2CaMCa4CaN_dyad - Ca4CaMCa4CaN_dyad 
    rcnCa4CaN_dyad = kcanCaon*Ca_Dyad^2*Ca2CaN_dyad - kcanCaoff*Ca4CaN_dyad 
    rcn02CaN_dyad = k02can*Ca_Dyad^2*CaMCa4CaN_dyad - k20can*Ca2CaMCa4CaN_dyad 
    rcn24CaN_dyad = k24can*Ca_Dyad^2*Ca2CaMCa4CaN_dyad - k42can*Ca4CaMCa4CaN_dyad 
    rcn0CaN_dyad = kcanCaM0on*CaM_dyad*Ca4CaN_dyad - kcanCaM0off*CaMCa4CaN_dyad 
    rcn2CaN_dyad = kcanCaM2on*Ca2CaM_dyad*Ca4CaN_dyad - kcanCaM2off*Ca2CaMCa4CaN_dyad 
    rcn4CaN_dyad = kcanCaM4on*Ca4CaM_dyad*Ca4CaN_dyad - kcanCaM4off*Ca4CaMCa4CaN_dyad 
    # CaMKII reaction fluxes
    Pi_dyad = 1 - Pb2_dyad - Pb_dyad - Pt_dyad - Pt2_dyad - Pa_dyad 
    rcnCKib2_dyad = kib2*Ca2CaM_dyad*Pi_dyad - kb2i*Pb2_dyad 
    rcnCKb2b_dyad = kb24*Ca_Dyad^2*Pb2_dyad - kb42*Pb_dyad 
    rcnCKib_dyad = kib*Ca4CaM_dyad*Pi_dyad - kbi*Pb_dyad 
    T_dyad = Pb_dyad + Pt_dyad + Pt2_dyad + Pa_dyad 
    kbt_dyad = 0.055*T_dyad + 0.0074*T_dyad^2 + 0.015*T_dyad^3 
    rcnCKbt_dyad = kbt_dyad*Pb_dyad - kpp1*PP1tot_dyad*Pt_dyad/(Kmpp1+CaMKIItot_dyad*Pt_dyad) 
    rcnCKtt2_dyad = kt42*Pt_dyad - kt24*Ca_Dyad^2*Pt2_dyad 
    rcnCKta_dyad = kta*Pt_dyad - kat*Ca4CaM_dyad*Pa_dyad 
    rcnCKt2a_dyad = kt2a*Pt2_dyad - kat2*Ca2CaM_dyad*Pa_dyad 
    rcnCKt2b2_dyad = kpp1*PP1tot_dyad*Pt2_dyad/(Kmpp1+CaMKIItot_dyad*Pt2_dyad) 
    rcnCKai_dyad = kpp1*PP1tot_dyad*Pa_dyad/(Kmpp1+CaMKIItot_dyad*Pa_dyad)

    ## SL Fluxes
    
    Ca_SL = Ca_sl * 1e3
    # CaM Reaction fluxes
    B_sl = Btot - CaMB_sl - Ca2CaMB_sl - Ca4CaMB_sl
    rcn02_sl = k02*Ca_SL^2*CaM_sl - k20*Ca2CaM_sl
    rcn24_sl = k24*Ca_SL^2*Ca2CaM_sl - k42*Ca4CaM_sl 
    # CaM buffer fluxes
    rcn02B_sl = k02B*Ca_SL^2*CaMB_sl - k20B*Ca2CaMB_sl 
    rcn24B_sl = k24B*Ca_SL^2*Ca2CaMB_sl - k42B*Ca4CaMB_sl 
    rcn0B_sl = k0Bon*CaM_sl*B_sl - k0Boff*CaMB_sl 
    rcn2B_sl = k2Bon*Ca2CaM_sl*B_sl - k2Boff*Ca2CaMB_sl 
    rcn4B_sl = k4Bon*Ca4CaM_sl*B_sl - k4Boff*Ca4CaMB_sl 
    # CaN reaction fluxes 
    Ca2CaN_sl = CaNtot - Ca4CaN_sl - CaMCa4CaN_sl - Ca2CaMCa4CaN_sl - Ca4CaMCa4CaN_sl 
    rcnCa4CaN_sl = kcanCaon*Ca_SL^2*Ca2CaN_sl - kcanCaoff*Ca4CaN_sl 
    rcn02CaN_sl = k02can*Ca_SL^2*CaMCa4CaN_sl - k20can*Ca2CaMCa4CaN_sl 
    rcn24CaN_sl = k24can*Ca_SL^2*Ca2CaMCa4CaN_sl - k42can*Ca4CaMCa4CaN_sl 
    rcn0CaN_sl = kcanCaM0on*CaM_sl*Ca4CaN_sl - kcanCaM0off*CaMCa4CaN_sl 
    rcn2CaN_sl = kcanCaM2on*Ca2CaM_sl*Ca4CaN_sl - kcanCaM2off*Ca2CaMCa4CaN_sl 
    rcn4CaN_sl = kcanCaM4on*Ca4CaM_sl*Ca4CaN_sl - kcanCaM4off*Ca4CaMCa4CaN_sl 
    # CaMKII reaction fluxes
    Pi_sl = 1 - Pb2_sl - Pb_sl - Pt_sl - Pt2_sl - Pa_sl 
    rcnCKib2_sl = kib2*Ca2CaM_sl*Pi_sl - kb2i*Pb2_sl 
    rcnCKb2b_sl = kb24*Ca_SL^2*Pb2_sl - kb42*Pb_sl 
    rcnCKib_sl = kib*Ca4CaM_sl*Pi_sl - kbi*Pb_sl 
    T_sl = Pb_sl + Pt_sl + Pt2_sl + Pa_sl 
    kbt_sl = 0.055*T_sl + 0.0074*T_sl^2 + 0.015*T_sl^3 
    rcnCKbt_sl = kbt_sl*Pb_sl - kpp1*PP1tot*Pt_sl/(Kmpp1+CaMKIItot*Pt_sl) 
    rcnCKtt2_sl = kt42*Pt_sl - kt24*Ca_SL^2*Pt2_sl 
    rcnCKta_sl = kta*Pt_sl - kat*Ca4CaM_sl*Pa_sl 
    rcnCKt2a_sl = kt2a*Pt2_sl - kat2*Ca2CaM_sl*Pa_sl 
    rcnCKt2b2_sl = kpp1*PP1tot*Pt2_sl/(Kmpp1+CaMKIItot*Pt2_sl) 
    rcnCKai_sl = kpp1*PP1tot*Pa_sl/(Kmpp1+CaMKIItot*Pa_sl)

    ## Cyt Fluxes
    
    Ca_Cyt = Cai * 1e3
    # CaM Reaction fluxes
    B_cyt = Btot - CaMB_cyt - Ca2CaMB_cyt - Ca4CaMB_cyt
    rcn02_cyt = k02*Ca_Cyt^2*CaM_cyt - k20*Ca2CaM_cyt 
    rcn24_cyt = k24*Ca_Cyt^2*Ca2CaM_cyt - k42*Ca4CaM_cyt 
    # CaM buffer fluxes
    rcn02B_cyt = k02B*Ca_Cyt^2*CaMB_cyt - k20B*Ca2CaMB_cyt 
    rcn24B_cyt = k24B*Ca_Cyt^2*Ca2CaMB_cyt - k42B*Ca4CaMB_cyt 
    rcn0B_cyt = k0Bon*CaM_cyt*B_cyt - k0Boff*CaMB_cyt 
    rcn2B_cyt = k2Bon*Ca2CaM_cyt*B_cyt - k2Boff*Ca2CaMB_cyt 
    rcn4B_cyt = k4Bon*Ca4CaM_cyt*B_cyt - k4Boff*Ca4CaMB_cyt 
    # CaN reaction fluxes 
    Ca2CaN_cyt = CaNtot - Ca4CaN_cyt - CaMCa4CaN_cyt - Ca2CaMCa4CaN_cyt - Ca4CaMCa4CaN_cyt 
    rcnCa4CaN_cyt = kcanCaon*Ca_Cyt^2*Ca2CaN_cyt - kcanCaoff*Ca4CaN_cyt 
    rcn02CaN_cyt = k02can*Ca_Cyt^2*CaMCa4CaN_cyt - k20can*Ca2CaMCa4CaN_cyt 
    rcn24CaN_cyt = k24can*Ca_Cyt^2*Ca2CaMCa4CaN_cyt - k42can*Ca4CaMCa4CaN_cyt 
    rcn0CaN_cyt = kcanCaM0on*CaM_cyt*Ca4CaN_cyt - kcanCaM0off*CaMCa4CaN_cyt 
    rcn2CaN_cyt = kcanCaM2on*Ca2CaM_cyt*Ca4CaN_cyt - kcanCaM2off*Ca2CaMCa4CaN_cyt 
    rcn4CaN_cyt = kcanCaM4on*Ca4CaM_cyt*Ca4CaN_cyt - kcanCaM4off*Ca4CaMCa4CaN_cyt 
    # CaMKII reaction fluxes
    Pi_cyt = 1 - Pb2_cyt - Pb_cyt - Pt_cyt - Pt2_cyt - Pa_cyt 
    rcnCKib2_cyt = kib2*Ca2CaM_cyt*Pi_cyt - kb2i*Pb2_cyt 
    rcnCKb2b_cyt = kb24*Ca_Cyt^2*Pb2_cyt - kb42*Pb_cyt 
    rcnCKib_cyt = kib*Ca4CaM_cyt*Pi_cyt - kbi*Pb_cyt 
    T_cyt = Pb_cyt + Pt_cyt + Pt2_cyt + Pa_cyt 
    kbt_cyt = 0.055*T_cyt + 0.0074*T_cyt^2 + 0.015*T_cyt^3 
    #rcnCKbt_cyt = Pb_cyt - kpp1*PP1tot*Pt_cyt/(Kmpp1+CaMKIItot*Pt_cyt) 
    rcnCKbt_cyt = kbt_cyt*Pb_cyt - kpp1*PP1tot*Pt_cyt/(Kmpp1+CaMKIItot*Pt_cyt) 
    rcnCKtt2_cyt = kt42*Pt_cyt - kt24*Ca_Cyt^2*Pt2_cyt 
    rcnCKta_cyt = kta*Pt_cyt - kat*Ca4CaM_cyt*Pa_cyt 
    rcnCKt2a_cyt = kt2a*Pt2_cyt - kat2*Ca2CaM_cyt*Pa_cyt 
    rcnCKt2b2_cyt = kpp1*PP1tot*Pt2_cyt/(Kmpp1+CaMKIItot*Pt2_cyt) 
    rcnCKai_cyt = kpp1*PP1tot*Pa_cyt/(Kmpp1+CaMKIItot*Pa_cyt)
    
    
    ## Ordinary Differential Equations
    Vmyo = 2.1454e-11           # [L]
    Vdyad = 1.7790e-014         # [L]
    VSL = 6.6013e-013           # [L]
    kSLmyo = 8.587e-15          # [L/msec]
    CaMKIItotDyad = 120         # [uM]
    BtotDyad = 1.54/8.293e-4    # [uM]
    CaMtotDyad = CaM_dyad+Ca2CaM_dyad+Ca4CaM_dyad+CaMB_dyad+Ca2CaMB_dyad+Ca4CaMB_dyad+CaMKIItotDyad*(Pb2_dyad+Pb_dyad+Pt_dyad+Pt2_dyad)+CaMCa4CaN_dyad+Ca2CaMCa4CaN_dyad+Ca4CaMCa4CaN_dyad
    Bdyad = BtotDyad - CaMtotDyad                                                  # [uM dyad]
    J_cam_dyadSL = 1e-3*(k0Boff*CaM_dyad - k0Bon*Bdyad*CaM_sl)                     # [uM/msec dyad]
    J_ca2cam_dyadSL = 1e-3*(k2Boff*Ca2CaM_dyad - k2Bon*Bdyad*Ca2CaM_sl)            # [uM/msec dyad]
    J_ca4cam_dyadSL = 1e-3*(k2Boff*Ca4CaM_dyad - k4Bon*Bdyad*Ca4CaM_sl)            # [uM/msec dyad]
    J_cam_SLmyo = kSLmyo*(CaM_sl-CaM_cyt)                                          # [umol/msec]
    J_ca2cam_SLmyo = kSLmyo*(Ca2CaM_sl-Ca2CaM_cyt)                                 # [umol/msec]
    J_ca4cam_SLmyo = kSLmyo*(Ca4CaM_sl-Ca4CaM_cyt)                                 # [umol/msec]

    # CaMDyad equations
    CaMdyad_eqs = [
        D(CaM_dyad) ~ (1e-3*(-rcn02_dyad - rcn0B_dyad - rcn0CaN_dyad)-J_cam_dyadSL),                                                                      # du[1]
        D(Ca2CaM_dyad) ~ (1e-3*(rcn02_dyad - rcn24_dyad - rcn2B_dyad - rcn2CaN_dyad + CaMKIItot_dyad.*(-rcnCKib2_dyad + rcnCKt2a_dyad))-J_ca2cam_dyadSL),      # du[2]
        D(Ca4CaM_dyad) ~ (1e-3*(rcn24_dyad - rcn4B_dyad - rcn4CaN_dyad + CaMKIItot_dyad.*(-rcnCKib_dyad+rcnCKta_dyad))-J_ca4cam_dyadSL),                       # du[3]
        D(CaMB_dyad) ~ 1e-3*(rcn0B_dyad-rcn02B_dyad),                       # du[4]
        D(Ca2CaMB_dyad) ~ 1e-3*(rcn02B_dyad + rcn2B_dyad - rcn24B_dyad),    # du[5]
        D(Ca4CaMB_dyad) ~ 1e-3*(rcn24B_dyad + rcn4B_dyad),                  # du[6]
        # CaMKII equations
        D(Pb2_dyad) ~ 1e-3*(rcnCKib2_dyad - rcnCKb2b_dyad + rcnCKt2b2_dyad),     # du[7]
        D(Pb_dyad) ~ 1e-3*(rcnCKib_dyad + rcnCKb2b_dyad - rcnCKbt_dyad),         # du[8]
        D(Pt_dyad) ~ 1e-3*(rcnCKbt_dyad-rcnCKta_dyad-rcnCKtt2_dyad),             # du[9]
        D(Pt2_dyad) ~ 1e-3*(rcnCKtt2_dyad-rcnCKt2a_dyad-rcnCKt2b2_dyad),         # du[10]
        D(Pa_dyad) ~ 1e-3*(rcnCKta_dyad+rcnCKt2a_dyad-rcnCKai_dyad),             # du[11]
        # CaN equations
        D(Ca4CaN_dyad) ~ 1e-3*(rcnCa4CaN_dyad - rcn0CaN_dyad - rcn2CaN_dyad - rcn4CaN_dyad),      # du[12]
        D(CaMCa4CaN_dyad) ~ 1e-3*(rcn0CaN_dyad - rcn02CaN_dyad),                        # du[13]
        D(Ca2CaMCa4CaN_dyad) ~ 1e-3*(rcn2CaN_dyad+rcn02CaN_dyad-rcn24CaN_dyad),              # du[14]
        D(Ca4CaMCa4CaN_dyad) ~ 1e-3*(rcn4CaN_dyad+rcn24CaN_dyad)                       # du[15]
    ]
    # CaMSL equations
    CaMSL_eqs = [
        D(CaM_sl) ~ (1e-3*(-rcn02_sl - rcn0B_sl - rcn0CaN_sl)+ J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL),                                                                     # du[1]
        D(Ca2CaM_sl) ~ (1e-3*(rcn02_sl - rcn24_sl - rcn2B_sl - rcn2CaN_sl + CaMKIItot.*(-rcnCKib2_sl + rcnCKt2a_sl))+ J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL),        # du[2]
        D(Ca4CaM_sl) ~ (1e-3*(rcn24_sl - rcn4B_sl - rcn4CaN_sl + CaMKIItot.*(-rcnCKib_sl+rcnCKta_sl))+ J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL),                       # du[3]
        D(CaMB_sl) ~ 1e-3*(rcn0B_sl-rcn02B_sl),                     # du[4]
        D(Ca2CaMB_sl) ~ 1e-3*(rcn02B_sl + rcn2B_sl - rcn24B_sl),    # du[5]
        D(Ca4CaMB_sl) ~ 1e-3*(rcn24B_sl + rcn4B_sl),                # du[6]
        # CaMKII equations
        D(Pb2_sl) ~ 1e-3*(rcnCKib2_sl - rcnCKb2b_sl + rcnCKt2b2_sl),     # du[7]
        D(Pb_sl) ~ 1e-3*(rcnCKib_sl + rcnCKb2b_sl - rcnCKbt_sl),         # du[8]
        D(Pt_sl) ~ 1e-3*(rcnCKbt_sl-rcnCKta_sl-rcnCKtt2_sl),             # du[9]
        D(Pt2_sl) ~ 1e-3*(rcnCKtt2_sl-rcnCKt2a_sl-rcnCKt2b2_sl),         # du[10]
        D(Pa_sl) ~ 1e-3*(rcnCKta_sl+rcnCKt2a_sl-rcnCKai_sl),             # du[11]
        # CaN equations
        D(Ca4CaN_sl) ~ 1e-3*(rcnCa4CaN_sl - rcn0CaN_sl - rcn2CaN_sl - rcn4CaN_sl),      # du[12]
        D(CaMCa4CaN_sl) ~ 1e-3*(rcn0CaN_sl - rcn02CaN_sl),                        # du[13]
        D(Ca2CaMCa4CaN_sl) ~ 1e-3*(rcn2CaN_sl+rcn02CaN_sl-rcn24CaN_sl),              # du[14]
        D(Ca4CaMCa4CaN_sl) ~ 1e-3*(rcn4CaN_sl+rcn24CaN_sl)                       # du[15]
    ]
    # CaMCyt equations
    CaMcyt_eqs = [
        D(CaM_cyt) ~ (1e-3*(-rcn02_cyt - rcn0B_cyt - rcn0CaN_cyt) + J_cam_SLmyo/Vmyo),                                                                    # du[1]
        D(Ca2CaM_cyt) ~ (1e-3*(rcn02_cyt - rcn24_cyt - rcn2B_cyt - rcn2CaN_cyt + CaMKIItot.*(-rcnCKib2_cyt + rcnCKt2a_cyt)) + J_ca2cam_SLmyo/Vmyo),       # du[2]
        D(Ca4CaM_cyt) ~ (1e-3*(rcn24_cyt - rcn4B_cyt - rcn4CaN_cyt + CaMKIItot.*(-rcnCKib_cyt+rcnCKta_cyt)) + J_ca4cam_SLmyo/Vmyo),                       # du[3]
        D(CaMB_cyt) ~ 1e-3*(rcn0B_cyt-rcn02B_cyt),                      # du[4]
        D(Ca2CaMB_cyt) ~ 1e-3*(rcn02B_cyt + rcn2B_cyt - rcn24B_cyt),    # du[5]
        D(Ca4CaMB_cyt) ~ 1e-3*(rcn24B_cyt + rcn4B_cyt),                 # du[6]
        # CaMKII equations
        D(Pb2_cyt) ~ 1e-3*(rcnCKib2_cyt - rcnCKb2b_cyt + rcnCKt2b2_cyt),     # du[7]
        D(Pb_cyt) ~ 1e-3*(rcnCKib_cyt + rcnCKb2b_cyt - rcnCKbt_cyt),         # du[8]
        D(Pt_cyt) ~ 1e-3*(rcnCKbt_cyt-rcnCKta_cyt-rcnCKtt2_cyt),             # du[9]
        D(Pt2_cyt) ~ 1e-3*(rcnCKtt2_cyt-rcnCKt2a_cyt-rcnCKt2b2_cyt),         # du[10]
        D(Pa_cyt) ~ 1e-3*(rcnCKta_cyt+rcnCKt2a_cyt-rcnCKai_cyt),             # du[11]
        # CaN equations
        D(Ca4CaN_cyt) ~ 1e-3*(rcnCa4CaN_cyt - rcn0CaN_cyt - rcn2CaN_cyt - rcn4CaN_cyt),      # du[12]
        D(CaMCa4CaN_cyt) ~ 1e-3*(rcn0CaN_cyt - rcn02CaN_cyt),                        # du[13]
        D(Ca2CaMCa4CaN_cyt) ~ 1e-3*(rcn2CaN_cyt+rcn02CaN_cyt-rcn24CaN_cyt),              # du[14]
        D(Ca4CaMCa4CaN_cyt) ~ 1e-3*(rcn4CaN_cyt+rcn24CaN_cyt)                       # du[15]
    ]
    ## For adjusting Ca buffering in EC coupling model
    JCaDyad = 1e-3*(2*CaMKIItot_dyad*(rcnCKtt2_dyad-rcnCKb2b_dyad) - 2*(rcn02_dyad+rcn24_dyad+rcn02B_dyad+rcn24B_dyad+rcnCa4CaN_dyad+rcn02CaN_dyad+rcn24CaN_dyad))   # [uM/msec]
    JCaSL = 1e-3*(2*CaMKIItot*(rcnCKtt2_sl-rcnCKb2b_sl) - 2*(rcn02_sl+rcn24_sl+rcn02B_sl+rcn24B_sl+rcnCa4CaN_sl+rcn02CaN_sl+rcn24CaN_sl))   # [uM/msec]
    JCaCyt = 1e-3*(2*CaMKIItot*(rcnCKtt2_cyt-rcnCKb2b_cyt) - 2*(rcn02_cyt+rcn24_cyt+rcn02B_cyt+rcn24B_cyt+rcnCa4CaN_cyt+rcn02CaN_cyt+rcn24CaN_cyt))   # [uM/msec]


## CaMKII
    ## RATE CONSTANTS and KM VALUES
    # L-Type Ca Channel (LTCC) parameters
    k_ckLCC = 0.4                   # [s^-1]
    k_pp1LCC = 0.1103               # [s^-1] 
    k_pkaLCC = 13.5                 # [s^-1] 
    k_pp2aLCC = 10.1                # [s^-1] 
    KmCK_LCC = 12                   # [uM] 
    KmPKA_LCC = 21                  # [uM] 
    KmPP2A_LCC = 47                 # [uM] 
    KmPP1_LCC = 9                   # [uM] 

    # Ryanodine Receptor (RyR) parameters
    k_ckRyR = 0.4                   # [s^-1] 
    k_pkaRyR = 1.35                 # [s^-1] 
    k_pp1RyR = 1.07                 # [s^-1] 
    k_pp2aRyR = 0.481               # [s^-1] 

    # Basal RyR phosphorylation (numbers based on param estimation)
    kb_2809 = 0.51                  # [uM/s] - PKA site
    kb_2815 = 0.35                  # [uM/s] - CaMKII site

    KmCK_RyR = 12                   # [uM] 
    KmPKA_RyR = 21                  # [uM] 
    KmPP1_RyR = 9                   # [uM] 
    KmPP2A_RyR = 47                 # [uM] 

    # Phospholamban (PLB) parameters
    k_ckPLB = 8e-3                  # [s^-1]
    k_pp1PLB = 0.0428               # [s^-1]

    KmCK_PLB = 12
    KmPP1_PLB = 9

    # Okadaic Acid inhibition params (based on Huke/Bers [2008])
    # Want to treat OA as non-competitive inhibitor of PP1 and PP2A
    Ki_OA_PP1 = 0.78                # [uM] - Values from fit
    Ki_OA_PP2A = 0.037              # [uM] - Values from fit

    # Default PKA level
    PKAc = 95.6 * 0.54

    ## Parameters for CaMKII module
    LCCtotDyad = 31.4*0.9      # [uM] - Total Dyadic [LCC] - (umol/l dyad)
    LCCtotSL = 0.0846          # [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
    RyRtot = 382.6             # [uM] - Total RyR (in Dyad)
    PP1_dyad = 95.7            # [uM] - Total dyadic [PP1]
    PP1_SL = 0.57              # [uM] - Total Subsarcolemmal [PP1]
    PP2A_dyad = 95.76          # [uM] - Total dyadic PP2A
    OA = 0                     # [uM] - PP1/PP2A inhibitor Okadaic Acid
    plb_val=106 # MOUSE
    PLBtot = plb_val           # [uM] - Total [PLB] in cytosolic units

    ## OA inhibition term (non-competitive) for PP1 and PP2A
    OA_PP1 = 1/(1 + (OA/Ki_OA_PP1)^3)
    OA_PP2A = 1/(1 + (OA/Ki_OA_PP2A)^3)

    CaMKIItotDyad = 120             # [uM]
    CaMKIItotSL = 120*8.293e-4      # [uM]
    PP1_PLBtot = 0.89               # [uM] - [umol/L cytosol]

    ## ODE EQUATIONS
    # LTCC states (note: PP2A is acting on PKA site and PP1 on CKII site)

    # Variables related to camdyad_ODEfile
    CaMKIIact_Dyad = CaMKIItotDyad .* (Pb_dyad + Pt_dyad + Pt2_dyad + Pa_dyad)
    CaMKIIact_SL = CaMKIItotSL .* (Pb_sl + Pt_sl + Pt2_sl + Pa_sl)
    PP1_PLB_avail = 1 - I1p_PP1/PP1_PLBtot + 0.081698
    # CaMKII phosphorylation of Dyadic LCCs
    LCC_CKdyadn = LCCtotDyad - LCC_CKdyadp 
    LCCDyad_PHOS = (k_ckLCC*CaMKIIact_Dyad*LCC_CKdyadn)/(KmCK_LCC+LCC_CKdyadn) 
    LCCDyad_DEPHOS = (k_pp1LCC*PP1_dyad*LCC_CKdyadp)/(KmPP1_LCC+LCC_CKdyadp)*OA_PP1 
    LCC_CKsln = LCCtotSL - LCC_CKslp 
    LCCSL_PHOS = (k_ckLCC*CaMKIIact_SL*LCC_CKsln)/(KmCK_LCC+LCC_CKsln) 
    LCCSL_DEPHOS = (k_pp1LCC*PP1_SL*LCC_CKslp)/(KmPP1_LCC+LCC_CKslp)*OA_PP1 
    LCC_PKAn = LCCtotDyad - LCC_PKAp 
    RyR2815n = RyRtot - RyR2815p 
    RyR_BASAL = kb_2815*RyR2815n 
    RyR_PHOS = (k_ckRyR*CaMKIIact_Dyad*RyR2815n)/(KmCK_RyR+RyR2815n) 
    RyR_PP1_DEPHOS = (k_pp1RyR*PP1_dyad*RyR2815p)/(KmPP1_RyR+RyR2815p)*OA_PP1 
    RyR_PP2A_DEPHOS = (k_pp2aRyR*PP2A_dyad*RyR2815p)/(KmPP2A_RyR+RyR2815p)*OA_PP2A 
    RyR2809n = RyRtot - RyR2809p 
    PP1_PLB = PP1_dyad*PP1_PLB_avail  # Inhibitor-1 regulation of PP1_dyad included here
    PLBT17n = PLBtot - PLBT17p 
    PLB_PHOS = (k_ckPLB*PLBT17n*CaMKIIact_Dyad)/(KmCK_PLB+PLBT17n) 
    PLB_DEPHOS = (k_pp1PLB*PP1_PLB*PLBT17p)/(KmPP1_PLB+PLBT17p)*OA_PP1 

    CaMKII_eqs = [
        D(LCC_CKdyadp) ~ (LCCDyad_PHOS - LCCDyad_DEPHOS)*1e-3, # du[2]
        # CaMKII phosphorylation of Sub-sarcolemmal LCCs
        D(LCC_CKslp) ~ (LCCSL_PHOS - LCCSL_DEPHOS)*1e-3, # du[6]
        # PKA phosphorylation (currently unused elsewhere)
        D(LCC_PKAp) ~ ((k_pkaLCC*PKAc*LCC_PKAn)/(KmPKA_LCC+LCC_PKAn) - (k_pp2aLCC*PP2A_dyad*LCC_PKAp)/(KmPP2A_LCC+LCC_PKAp)*OA_PP2A)*1e-3, # du[1]
        # RyR states
        D(RyR2815p) ~ (RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS)*1e-3, # du[4]
        # PKA phosphorylation of Ser 2809 on RyR (currently unused elsewhere)
        D(RyR2809p) ~ (kb_2809*RyR2809n + (k_pkaRyR*PKAc*RyR2809n)/(KmPKA_RyR+RyR2809n) - (k_pp1RyR*PP1_dyad*RyR2809p)/(KmPP1_RyR+RyR2809p)*OA_PP1)*1e-3, # du[3]
        # PLB states
        D(PLBT17p) ~ (PLB_PHOS - PLB_DEPHOS)*1e-3 # du[5]
    ]


## BAR
    # Drug concentrations
    Ligtot = 0.0               # [uM] - SET LIGAND CONCENTRATION (0 or 0.1)
    FSK=0
    IBMX=0
    LCCtotBA = 0.025           # [uM] - [umol/L cytosol]
    plb_val=106                # MOUSE
    PP1_PLBtot = 0.89          # [uM] - [umol/L cytosol]
    PLMtotBA = 48              # [uM] - [umol/L cytosol] MOUSE
    PLBtotBA = plb_val                     # [uM] - [umol/L cytosol]
    ISO=Ligtot

    ## b-AR module
    b1ARtot         = 0.00528        # (uM) total b1-AR protein # MOUSE
    #b1ARtot=0.028  # RABBIT
    kf_LR           = 1              # (1/[uM ms]) forward rate for ISO binding to b1AR
    kr_LR           = 0.285          # (1/ms) reverse rate for ISO binding to b1AR
    kf_LRG          = 1              # (1/[uM ms]) forward rate for ISO:b1AR association with Gs
    kr_LRG          = 0.062          # (1/ms) reverse rate for ISO:b1AR association with Gs
    kf_RG           = 1              # (1/[uM ms]) forward rate for b1AR association with Gs
    kr_RG           = 33             # (1/ms) reverse rate for b1AR association with Gs
    Gstot           = 3.83           # (uM) total Gs protein
    k_G_act         = 16e-3          # (1/ms) rate constant for Gs activation
    k_G_hyd         = 0.8e-3         # (1/ms) rate constant for G-protein hydrolysis
    k_G_reassoc     = 1.21           # (1/[uM ms]) rate constant for G-protein reassociation
    kf_bARK         = 1.1e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
    kr_bARK         = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by b1ARK
    kf_PKA          = 3.6e-6         # (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
    kr_PKA          = 2.2e-6         # (1/ms) reverse rate for b1AR phosphorylation by PKA
    b1ARact = b1ARtot - b1AR_S464 - b1AR_S301
    b1AR = b1ARact - LR - LRG - RG
    Gs = Gstot - LRG - RG - Gsby
    bARK_desens = kf_bARK*(LR+LRG)
    bARK_resens = kr_bARK*b1AR_S464
    PKA_desens = kf_PKA*PKACI*b1ARact
    PKA_resens = kr_PKA*b1AR_S301
    G_act = k_G_act*(RG+LRG)
    G_hyd = k_G_hyd*GsaGTPtot
    G_reassoc = k_G_reassoc*GsaGDP*Gsby

    bar_eqs = [
        D(LR) ~ (kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs),    # du[1]
        D(LRG) ~ (kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG),                 # du[2]
        D(RG) ~ (kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG),                    # du[3]
        D(b1AR_S464) ~ (bARK_desens - bARK_resens),                         # du[4]
        D(b1AR_S301) ~ (PKA_desens - PKA_resens),                           # du[5]
        D(GsaGTPtot) ~ (G_act - G_hyd),                                     # du[6]
        D(GsaGDP) ~ (G_hyd - G_reassoc),                                    # du[7]
        D(Gsby) ~ (G_act - G_reassoc)                                       # du[8]
    ]

    ## cAMP module
    ACtot           = 70.57e-3        # (uM) total adenylyl cyclase # MOUSE
    # ACtot=47e-3  # RABBIT
    ATP             = 5e3             # (uM) total ATP
    k_AC_basal      = 0.2e-3          # (1/ms) basal cAMP generation rate by AC
    Km_AC_basal     = 1.03e3          # (uM) basal AC affinity for ATP
    Kd_AC_Gsa       = 0.4             # (uM) Kd for AC association with Gsa
    kf_AC_Gsa       = 1               # (1/[uM ms]) forward rate for AC association with Gsa
    kr_AC_Gsa       = Kd_AC_Gsa       # (1/ms) reverse rate for AC association with Gsa
    k_AC_Gsa        = 8.5e-3          # (1/ms) basal cAMP generation rate by AC:Gsa
    Km_AC_Gsa       = 315.0           # (uM) AC:Gsa affinity for ATP
    Kd_AC_FSK       = 44.0            # (uM) Kd for FSK binding to AC
    k_AC_FSK        = 7.3e-3          # (1/ms) basal cAMP generation rate by AC:FSK
    Km_AC_FSK       = 860.0           # (uM) AC:FSK affinity for ATP
    PDEtot          = 22.85e-3        # (uM) total phosphodiesterase
    k_cAMP_PDE      = 5e-3            # (1/ms) cAMP hydrolysis rate by PDE
    k_cAMP_PDEp     = 2*k_cAMP_PDE    # (1/ms) cAMP hydrolysis rate by phosphorylated PDE
    Km_PDE_cAMP     = 1.3             # (uM) PDE affinity for cAMP
    Kd_PDE_IBMX     = 30.0            # (uM) Kd_R2cAMP_C for IBMX binding to PDE
    k_PKA_PDE       = 7.5e-3          # (1/ms) rate constant for PDE phosphorylation by type 1 PKA
    k_PP_PDE        = 1.5e-3          # (1/ms) rate constant for PDE dephosphorylation by phosphatases
    cAMP = cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II)
    AC = ACtot-AC_GsaGTP
    GsaGTP = GsaGTPtot - AC_GsaGTP
    AC_FSK = FSK*AC/Kd_AC_FSK
    AC_ACT_BASAL = k_AC_basal*AC*ATP/(Km_AC_basal+ATP)
    AC_ACT_GSA = k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP)
    AC_ACT_FSK = k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP)
    PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX
    PDE = PDEtot - PDE_IBMX - PDEp
    PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP)

    cAMP_eqs = [
        D(AC_GsaGTP) ~ (kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP),
        D(PDEp) ~ (k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp),
        D(cAMPtot) ~ (AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT)
    ]


    ## PKA module
    PKItot          = 0.18             # (uM) total PKI
    kf_RC_cAMP      = 1                # (1/[uM ms]) Kd for PKA RC binding to cAMP
    kf_RCcAMP_cAMP  = 1                # (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
    kf_RcAMPcAMP_C  = 4.375            # (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
    kf_PKA_PKI      = 1                # (1/[uM ms]) Ki for PKA inhibition by PKI
    kr_RC_cAMP      = 1.64             # (1/ms) Kd for PKA RC binding to cAMP
    kr_RCcAMP_cAMP  = 9.14             # (1/ms) Kd for PKA RC:cAMP binding to cAMP
    kr_RcAMPcAMP_C  = 1                # (1/ms) Kd for PKA R:cAMPcAMP binding to C
    kr_PKA_PKI      = 2e-4             # (1/ms) Ki for PKA inhibition by PKI
    epsilon         = 10               # (-) AKAP-mediated scaling factor
    PKI = PKItot - PKACI_PKI - PKACII_PKI

    PKA_eqs = [
        D(RC_I) ~ (- kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I),
        D(RCcAMP_I) ~ (- kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I),
        D(RCcAMPcAMP_I) ~ (- kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI),
        D(RcAMPcAMP_I) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I),
        D(PKACI) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI),
        D(PKACI_PKI) ~ (- kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI),
        D(RC_II) ~ (- kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II),
        D(RCcAMP_II) ~ (- kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II),
        D(RCcAMPcAMP_II) ~ (- kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII),
        D(RcAMPcAMP_II) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II),
        D(PKACII) ~ (- kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI),
        D(PKACII_PKI) ~ (- kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI)
    ]


    ## I-1/PP1 module
    I1tot           = 0.3             # (uM) total inhibitor 1
    k_PKA_I1        = 60e-3           # (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
    Km_PKA_I1       = 1.0             # (uM) Km for I-1 phosphorylation by type 1 PKA
    Vmax_PP2A_I1    = 14.0e-3         # (uM/ms) Vmax for I-1 dephosphorylation by PP2A
    Km_PP2A_I1      = 1.0             # (uM) Km for I-1 dephosphorylation by PP2A
    Ki_PP1_I1       = 1.0e-3          # (uM) Ki for PP1 inhibition by I-1
    kf_PP1_I1       = 1               # (uM) Ki for PP1 inhibition by I-1
    PP1tot                      = PP1_PLBtot      # PP1tot = 0.89  # (uM) total phosphatase 1
    kr_PP1_I1                   = Ki_PP1_I1       # (uM) Ki for PP1 inhibition by I-1
    I1 = I1tot - I1ptot
    PP1 = PP1tot - I1p_PP1
    I1p = I1ptot - I1p_PP1
    I1_phosph = k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1)
    I1_dephosph = Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot)

    PP1_eqs = [
        D(I1p_PP1) ~ (kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1),
        D(I1ptot) ~ (I1_phosph - I1_dephosph)
    ]


    ## PLB module
    PLBtot = PLBtotBA   # [uM]
    k_PKA_PLB = 54e-3   # [1/ms]
    Km_PKA_PLB = 21     # [uM]
    k_PP1_PLB = 8.5e-3  # [1/ms]
    Km_PP1_PLB = 7.0    # [uM]

    PLB = PLBtot - PLBp
    PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB)
    PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp)
    
    PLB_eqs = [
        D(PLBp) ~ (PLB_phosph - PLB_dephosph)
    ]


    ## PLM module (included 09/18/12) MOUSE
    PLMtot = PLMtotBA   # [uM]
    k_PKA_PLM = 54e-3   # [1/ms]
    Km_PKA_PLM = 21     # [uM]
    k_PP1_PLM = 8.5e-3  # [1/ms]
    Km_PP1_PLM = 7.0    # [uM]
    PLM = PLMtot - PLMp
    PLM_phosph = k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM)
    PLM_dephosph = k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp)
    
    PLM_eqs = [
        D(PLMp) ~ (PLM_phosph - PLM_dephosph)
    ]

    ## LCC module
    PKAIItot = 0.059        # (uM) total type 2 PKA # MOUSE
    LCCtot = LCCtotBA       # [uM]
    PKACII_LCCtot = 0.025   # [uM]
    PP1_LCC = 0.025         # [uM]
    PP2A_LCC = 0.025        # [uM]
    k_PKA_LCC = 54e-3       # [1/ms]
    Km_PKA_LCC = 21         # [uM]
    k_PP1_LCC = 8.52e-3     # [1/ms] RABBIT, MOUSE
    Km_PP1_LCC = 3          # [uM]
    k_PP2A_LCC = 10.1e-3    # [1/ms]
    Km_PP2A_LCC = 3         # [uM]
    PKACII_LCC = (PKACII_LCCtot/PKAIItot)*PKACII
    LCCa = LCCtot - LCCap
    LCCa_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa)
    LCCa_dephosph = epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap)
    LCCb = LCCtot - LCCbp
    LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb)
    LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp)
    
    LCC_eqs = [
        D(LCCap) ~ (LCCa_phosph - LCCa_dephosph),
        D(LCCbp) ~ (LCCb_phosph - LCCb_dephosph)
    ]

    ## RyR module (not included in Yang-Saucerman)
    PKAIIryrtot = 0.034         # [uM]
    PP1ryr = 0.034              # [uM]
    PP2Aryr = 0.034             # [uM]
    kcat_pka_ryr = 54e-3        # [1/ms]
    Km_pka_ryr = 21             # [uM]
    kcat_pp1_ryr = 8.52e-3      # [1/ms]
    Km_pp1_ryr = 7              # [uM]
    kcat_pp2a_ryr = 10.1e-3     # [1/ms]
    Km_pp2a_ryr = 4.1           # [uM]
    RyRtot = 0.135              # [uM]

    PKACryr = (PKAIIryrtot/PKAIItot)*PKACII
    RyR = RyRtot-RyRp
    RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyR/(Km_pka_ryr+epsilon*RyR)
    RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyRp/(Km_pp1_ryr+epsilon*RyRp)
    RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyRp/(Km_pp2a_ryr+epsilon*RyRp)

    RyRp_eqs = [
        D(RyRp) ~ (RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A)
    ]


    ## TnI module
    PP2A_TnI = 0.67         # [uM]
    k_PKA_TnI = 54e-3       # [1/ms]
    Km_PKA_TnI = 21         # [uM]
    k_PP2A_TnI = 10.1e-3    # [1/ms]
    Km_PP2A_TnI = 4.1       # [uM]
    TnItot = 70             # [uM] 
    TnI = TnItot - TnIp
    TnI_phosph = k_PKA_TnI*PKACI*TnI/(Km_PKA_TnI+TnI)
    TnI_dephosph = k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp)

    TnI_eqs = [
        D(TnIp) ~ (TnI_phosph - TnI_dephosph)
    ]
    

    ## Iks module (not present in mouse)
    IKstot = 0.025

    IKs_eqs = [
        D(KS79) ~ 0,  # ydot(27) not ODE
        D(KS80) ~ 0,  # ydot(28) not ODE
        D(KSp) ~ 0  # ydot(29)
    ]
    

    ## CFTR module (included 04/30/10)
    ICFTRtot = 0.025  

    CFTR_eqs = [
        D(CFTRp) ~ 0  #CFTRphos - CFTRdephos  # ydot(30)
    ]

    ## Ikur module (included 04/10/12) MOUSE
    PKAII_KURtot = 0.025    # [uM]
    PP1_KURtot = 0.025      # [uM]
    k_pka_KUR = 54e-3       # [1/ms]
    Km_pka_KUR = 21         # [uM]
    k_pp1_KUR = 8.52e-3     # [1/ms]
    Km_pp1_KUR = 7          # [uM]
    IKurtot = 0.025         # [uM]
    KURn = IKurtot - KURp   # Non-phos = tot - phos
    PKAC_KUR = (PKAII_KURtot/PKAIItot)*PKACII     # (PKA_KURtot/PKAIItot)*PKAIIact
    KURphos = epsilon*KURn*PKAC_KUR*k_pka_KUR/(Km_pka_KUR+epsilon*KURn) 
    KURdephos = PP1_KURtot*k_pp1_KUR*epsilon*KURp/(Km_pp1_KUR+epsilon*KURp) 
    Ikur_eqs = [
        D(KURp) ~ (KURphos - KURdephos)
    ]


    ## ECC     
    LCCtotDyad = 31.4*0.9       # [uM] - Total Dyadic [LCC] - (umol/l dyad)
    RyRtot = 382.6              # [uM] - Total RyR (in Dyad)
    plb_val=106                 # MOUSE
    LCCtotBA = 0.025            # [uM] - [umol/L cytosol]
    PLMtotBA = 48               # [uM] - [umol/L cytosol] MOUSE
    RyRtotBA = 0.135            # [uM] - [umol/L cytosol]
    TnItotBA = 70               # [uM] - [umol/L cytosol]
    IKurtotBA = 0.025           # [uM] - [umol/L cytosol] MOUSE
    PLBtot = plb_val            # [uM] - Total [PLB] in cytosolic units
    PLBtotBA = plb_val          # [uM] - [umol/L cytosol]

    LCC_CKp = LCC_CKdyadp ./ LCCtotDyad
    RyR_CKp = RyR2815p ./ RyRtot
    PLB_CKp = PLBT17p ./ PLBtot
    LCCa_PKAp = LCCap ./ LCCtotBA
    LCCb_PKAp = LCCbp ./ LCCtotBA
    PLB_PKAn = (PLBtotBA - PLBp) ./ PLBtotBA
    RyR_PKAp = RyRp ./ RyRtotBA
    TnI_PKAp = TnIp ./ TnItotBA
    IKur_PKAp = KURp ./ IKurtotBA
    PLM_PKAp = PLMp ./ PLMtotBA

    ## Flags
    CKIIOE = 0
    CKIIflag = CKIIOE
    ICa_MarkovFlag = 1
    INa_MarkovFlag = 0
    ItoFlag = 1
    NaClampFlag = 0
    PLMkoFlag = 0
    StrophFlag = 0
    CaffeineFlag = 0
    DigitalisFlag = 0

    # Na loading parameters ON (set flag to 1, 0 otherwise)
    NaGainFlag=0 # WT (default 0)
    if CKIIflag==1
        NaGainFlag=1 # CaMKII-OE (default 1)
    end

    # CaMKII-Na-Ca-CaMKII loop closed (set flag to 1, 0 otherwise)
    loop = 0 # (default 0)


    ## Model Parameters
    freq = 1                   # [Hz] - CHANGE DEPENDING ON FREQUENCY
    cycleLength = 1e3/freq     # [ms]
    # Constants
    R = 8314       # [J/kmol*K]  
    Frdy = 96485   # [C/mol]  
    Temp = 310     # [K] 310 K (37 C) for BT / 295 K (22 C) for RT
    FoRT = Frdy/R/Temp
    Qpow = (Temp-310)/10

    # Cell geometry
    Acell = 20e3 # [um^2] MOUSE
    Cmem = Acell*1e-14 # [F] 200 pF membrane capacitance MOUSE

    # Fractional currents in compartments
    Fjunc = 17/(17+31)*7/17+31/(17+31)*2/31 
    Fsl = 1-Fjunc
    Fjunc_nak = 1.6*17/(1.6*17+31)*7/17+31/(1.6*17+31)*2/31 
    Fsl_nak = 1-Fjunc_nak
    Fjunc_ncx = Fjunc 
    Fsl_ncx = 1-Fjunc_ncx
    Fjunc_CaL = 0.9 
    Fsl_CaL = 1-Fjunc_CaL
    cellLength = 100 # cell length [um]
    cellRadius = 10.25 # cell radius [um]
    junctionLength = 15e-3 # junc length [um]
    junctionRadius = 160e-3 # junc radius [um]
    distSLcyto = 0.45 # dist. SL to cytosol [um]
    distJuncSL = 0.3 # dist. junc to SL [um] MOUSE
    DcaJuncSL = 1.64e-6 # Dca junc to SL [cm^2/sec]
    DcaSLcyto = 1.22e-6 # Dca SL to cyto [cm^2/sec]
    DnaJuncSL = 1.09e-5 # Dna junc to SL [cm^2/sec]
    DnaSLcyto = 1.79e-5 # Dna SL to cyto [cm^2/sec] 
    Vcell = pi*cellRadius^2*cellLength*1e-15 # [L]
    Vmyo = 0.65*Vcell 
    Vsr = 0.035*Vcell 
    Vsl = 0.02*Vcell 
    Vjunc = 0.0539*0.01*Vcell
    SAsl = Fsl*Acell # [um^2]  MOUSE
    Njunc = (Fjunc*Acell)/(pi*junctionRadius^2)
    SAjunc = Njunc*pi*2*junctionLength*junctionRadius # [um^2] MOUSE
    J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10 # [L/msec] [m^2/sec] MOUSE
    J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10 # MOUSE
    J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10 # MOUSE
    J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10 # MOUSE

    # Fixed ion concentrations
    Cli = 15   # Intracellular Cl  [mM]
    Clo = 150  # Extracellular Cl  [mM]
    Ko = 5.4   # Extracellular K   [mM]
    Nao = 140  # Extracellular Na  [mM]
    Cao = 1    # Extracellular Ca  [mM]
    Mgi = 1    # Intracellular Mg  [mM]

    # Nernst Potentials
    ena_junc = (1/FoRT)*log(Nao/Naj)     # [mV]
    ena_sl = (1/FoRT)*log(Nao/Nasl)       # [mV]
    ek = (1/FoRT)*log(Ko/Ki)	           # [mV]
    eca_junc = (1/FoRT/2)*log(Cao/Ca_j)   # [mV]
    eca_sl = (1/FoRT/2)*log(Cao/Ca_sl)     # [mV]
    ecl = (1/FoRT)*log(Cli/Clo)            # [mV]


    # Na transport parameters
    GNa = 10 # [mS/uF]
    GNaB = 4.5*0.297e-3 # [mS/uF]
    IbarNaK = 5 # [uA/uF]
    if NaGainFlag == 1
        GNaB=GNaB*4
        IbarNaK=IbarNaK*0.9
    end
    if DigitalisFlag == 1
        IbarNaK = IbarNaK*0.5
    end
    KmNaip = 19 # [mM]
    KmKo = 1.5 # [mM]
    Q10NaK = 1.63
    Q10KmNai = 1.39
    if PLMkoFlag==1
        PLM_PKAp=1
        GNaB=GNaB*48/20
        IbarNaK=IbarNaK*0.8
    end
    if StrophFlag==1
        IbarNaK = 0
    end

    # INa Markov Model parameters
    GNa2 = 10.64      # [mS/uF]
    P1a1=3.802
    P2a1=0.1027
    P3a1=2.5
    P4a1=17
    P5a1=0.20
    P6a1=150
    P4a2=15
    P5a2=0.23
    P4a3=12
    P5a3=0.25
    P1b1=0.1917
    P2b1=20.3
    P1b2=0.2
    P2b2=2.5
    P1b3=0.22
    P2b3=7.5
    P1a4=0.188495
    P2a4=16.6
    P3a4=0.393956
    P4a4=7
    P1a5=7e-7
    P2a5=7.2 # TG 7.8
    P1b5=0.0044 # TG 0.00467
    P2b5=2e-5
    P1a6=100
    P1b6=8.9554e-7 # TG 6.5449e-7
    P2b6=11.3944
    P1a7=0.487e-4 # TG 0.3377e-3
    P2a7=23.2696
    P1b7=0.2868e-3 # TG 1.868e-4
    P2b7=35.9898
    P1a8=0.1e-7 # TG 6.5e-6
    P1b8=9.8e-3 # TG 3.8e-3
    if CKIIflag == 1 # MOUSE - CaMKII-OE
        P2a5=7.8
        P1b5=0.00467
        P1b6=6.5449e-7
        P1a7=0.3377e-3
        P1b7=1.868e-4
        P1a8=6.5e-6
        P1b8=3.8e-3
    end

    # K currents parameters
    Kcoeff = 1 # K current modulation
    pNaK = 0.01833
    GtoSlow = 0 # [mS/uF] NO ItoSlow in MOUSE
    GtoFast = 0.44 # [mS/uF] changed from rabbit (0.02)
    if CKIIflag == 1 # MOUSE
        GtoFast=GtoFast*2/3 # chronic CaMKII-OE effect
    end
    Gkur1 = 1.1*0.16 # fast
    Gkur2 = 0.14 # slow
    Gss = 0.15 # [mS/uF] only in MOUSE
    gkr = 0.03*sqrt(Ko/5.4)
    gkp = 0.001

    # Cl current parameters
    GClCa = 0.109625 # [mS/uF]
    GClB = 9e-3 # [mS/uF]
    KdClCa = 100e-3 # [mM]

    # LTCC parameters
    K_Ica = 1.65 # MOUSE
    pNa = K_Ica*1.5e-8 # [cm/sec]
    pCa = K_Ica*5.4e-4 # [cm/sec] - Ca permeability
    pK = K_Ica*2.7e-7 # [cm/sec]
    KmCa = 0.6e-3 # [mM]
    Q10CaL = 1.8

    # Ca transport parameters
    IbarNCX = 1 # [uA/uF] changed from rabbit (9)
    if CKIIflag == 1
        IbarNCX = 1.5*IbarNCX
    end
    KmCai = 3.59e-3    # [mM]
    KmCao = 1.3        # [mM]
    KmNai = 12.29      # [mM]
    KmNao = 87.5       # [mM]
    ksat = 0.27        # [none]  
    nu = 0.35          # [none]
    Kdact = 1/2*0.256e-3 # [mM] changed from rabbit
    Q10NCX = 1.57      # [none]
    IbarSLCaP = 0.0673 # [uA/uF]
    KmPCa = 0.5e-3     # [mM]
    GCaB = 3*2.513e-4 # [uA/uF] changed from rabbit (2.513e-4)
    Q10SLCaP = 2.35    # [none]

    # SR flux parameters
    Q10SRCaP = 2.6          # [none]
    Vmax_SRCaP = 1.15*1.15*2.86e-4 # [mM/msec] (mmol/L cytosol/msec) changed
    Kmf = 0.3e-3 # [mM] changed from rabbit (0.246e-3) # from Yang-Saucerman
    Kmr = 2.1 # [mM]L cytosol changed from rabbit (1.7) # from Yang-Saucerman
    hillSRCaP = 1.787       # [mM]
    ks = 25                 # [1/ms]      
    koCa = 10               # [mM^-2 1/ms]
    kom = 0.06              # [1/ms]     
    kiCa = 0.5              # [1/mM/ms]
    kim = 0.005             # [1/ms]
    ec50SR = 0.45+0.05      # [mM] changed from rabbit (0.45)
    if CaffeineFlag==1
        koCa=koCa*7.5
        GCaB=0
        Vmax_SRCaP=0
    end

    # Buffering parameters
    Bmax_Naj = 7.561       # [mM] 
    Bmax_Nasl = 1.65       # [mM]
    koff_na = 1e-3         # [1/ms]
    kon_na = 0.1e-3        # [1/mM/ms]
    Bmax_TnClow = 70e-3    # [mM]                      # TnC low affinity
    koff_tncl = 19.6e-3    # [1/ms] 
    kon_tncl = 32.7        # [1/mM/ms]
    Bmax_TnChigh = 140e-3  # [mM]                      # TnC high affinity 
    koff_tnchca = 0.032e-3 # [1/ms] 
    kon_tnchca = 2.37      # [1/mM/ms]
    koff_tnchmg = 3.33e-3  # [1/ms] 
    kon_tnchmg = 3e-3      # [1/mM/ms]
    Bmax_myosin = 140e-3   # [mM]                      # Myosin buffering
    koff_myoca = 0.46e-3   # [1/ms]
    kon_myoca = 13.8       # [1/mM/ms]
    koff_myomg = 0.057e-3  # [1/ms]
    kon_myomg = 0.0157     # [1/mM/ms]
    Bmax_SR = 19*.9e-3     # [mM] 
    koff_sr = 60e-3        # [1/ms]
    kon_sr = 100           # [1/mM/ms]
    Bmax_SLlowsl = 37.38e-3*Vmyo/Vsl        # [mM]    # SL buffering
    Bmax_SLlowj = 4.62e-3*Vmyo/Vjunc*0.1    # [mM]    
    koff_sll = 1300e-3     # [1/ms]
    kon_sll = 100          # [1/mM/ms]
    Bmax_SLhighsl = 13.35e-3*Vmyo/Vsl       # [mM] 
    Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1  # [mM] 
    koff_slh = 30e-3       # [1/ms]
    kon_slh = 100          # [1/mM/ms]
    Bmax_Csqn = 2.7        # [mM] 
    koff_csqn = 65         # [1/ms] 
    kon_csqn = 100         # [1/mM/ms] 

    # PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
    fracTnIpo = 0.062698
    fPKA_TnI = (1.61-0.61*(1-TnI_PKAp)/(1-fracTnIpo))
    koff_tncl = koff_tncl*fPKA_TnI


    ## I_Na: Fast Na Current
    # Max INa alterations with CaMKII hyperactivity as in Hund & Rudy 2008
    if CKIIflag == 1 # acute effects
        inashift = -3.25
        alphaCKII = -0.18
        if NaGainFlag == 1
            deltGbarNal_CKII = 3  # MOUSE
        else
            deltGbarNal_CKII = 0  # no Na Gain in OE
        end
    else
        inashift = 0
        alphaCKII = 0
        deltGbarNal_CKII = 0
    end

    if loop == 1
        RyRp_WT_mean = 0.2101
        RyRp_OE_mean = 0.7387
        RyRp_OEloop_min = 0.7033
        delta_loop=(3/(RyRp_OE_mean-RyRp_WT_mean))*RyR_CKp-(3/(RyRp_OE_mean-RyRp_WT_mean))*RyRp_WT_mean
        NaVsCaMKIIclamp=0
        if NaVsCaMKIIclamp==1
            delta_loop=(3/(RyRp_OE_mean-RyRp_WT_mean))*RyRp_OEloop_min-(3/(RyRp_OE_mean-RyRp_WT_mean))*RyRp_WT_mean
        end
        GNaB=(4.5)*0.297e-3*(1+delta_loop)
        if CKIIflag == 1 # OE
            if NaGainFlag == 1 # acute
                deltGbarNal_CKII=delta_loop
            else
                deltGbarNal_CKII=0
            end
        else # WT
            deltGbarNal_CKII=0
        end          
    end

    am = 0.32*(Vm+47.13)/(1-exp(-0.1*(Vm+47.13)))
    bm = 0.08*exp(-(Vm)/11)
    ah = ifelse((Vm-inashift) >= -40, 0, 0.135*exp((80+(Vm-inashift))/-6.8)) 
    aj = ifelse((Vm-inashift) >= -40, 0, (1+alphaCKII)*((-1.2714e5*exp(0.2444*(Vm-inashift))-3.474e-5*exp(-0.04391*(Vm-inashift)))*((Vm-inashift)+37.78)/(1+exp(0.311*((Vm-inashift)+79.23)))))
    bh = ifelse((Vm-inashift) >= -40, 0.66*1/(0.13*(1+exp(-((Vm-inashift)+10.66)/11.1))), 1.1*3.56*exp(0.079*(Vm-inashift-2))+3.1e5*exp(0.35*(Vm-inashift-2))) # MOUSE
    bj = ifelse((Vm-inashift) >= -40, 0.3*exp(-2.535e-7*(Vm-inashift))/(1+exp(-0.1*((Vm-inashift)+32))), 0.1212*exp(-0.01052*(Vm-inashift))/(1+exp(-0.1378*((Vm-inashift)+40.14))))

    INa_fast_eqs =[
        D(Na_m) ~ am*(1-Na_m)-bm*Na_m,
        D(Na_h) ~ ah*(1-Na_h)-bh*Na_h,
        D(Na_j) ~ aj*(1-Na_j)-bj*Na_j
    ]

    I_Na_junc1 = Fjunc*GNa*Na_m^3*Na_h*Na_j*(Vm-ena_junc)
    I_Na_sl1 = Fsl*GNa*Na_m^3*Na_h*Na_j*(Vm-ena_sl)


    ## I_Na,L: Late INa current (as in Hund & Rudy 2008)
    GbarNal = 0.0065*(1+deltGbarNal_CKII)*2 # deltGbar assigned in 'Fast INa' section
    # h-gate (note: m-gate is same as INa m-gate -> using Na_m for this)
    hlss = 1/(1+exp((Vm+91)/6.1))
    tauhl = 600 # ms

    Na_h = [
        D(Na_late_h) ~ (hlss-Na_late_h)/tauhl
    ]

    I_Nalj = Fjunc*GbarNal*Na_m^3*Na_late_h*(Vm-ena_junc)
    I_Nalsl = Fsl*GbarNal*Na_m^3*Na_late_h*(Vm-ena_sl)


    ## I_Na: alternative Markov Model - unused
    I2Na = (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3))
    # Transition rates
    alphaNa1 = P1a1/(P2a1*exp(-(Vm+P3a1)/P4a1)+P5a1*exp(-(Vm+P3a1)/P6a1))
    alphaNa2 = P1a1/(P2a1*exp(-(Vm+P3a1)/P4a2)+P5a2*exp(-(Vm+P3a1)/P6a1))
    alphaNa3 = P1a1/(P2a1*exp(-(Vm+P3a1)/P4a3)+P5a3*exp(-(Vm+P3a1)/P6a1))
    betaNa1 = P1b1*exp(-(Vm+P3a1)/P2b1)
    betaNa2 = P1b2*exp(-(Vm-P2b2)/P2b1)
    betaNa3 = P1b3*exp(-(Vm-P2b3)/P2b1)
    alphaNa4 = 1/(P1a4*exp(-(Vm+P4a4)/P2a4)+P3a4)
    alphaNa5 = P1a5*exp(-(Vm+P4a4)/P2a5)
    betaNa5 = (P1b5+P2b5*(Vm+P4a4))
    betaNa6 = P1b6*exp(-Vm/P2b6)
    alphaNa7 = P1a7*exp(Vm/P2a7)
    betaNa7 = P1b7*exp(-Vm/P2b7)
    alphaNa8 = P1a8
    betaNa8 = P1b8
    betaNa4 = (alphaNa3*alphaNa4*alphaNa5)/(betaNa3*betaNa5)
    alphaNa6 = alphaNa4/P1a6

    INa_eqs = [
        D(CNa3) ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LCNa3+betaNa1*CNa2+alphaNa5*ICNa3-(alphaNa1+betaNa5+alphaNa8)*CNa3), 0),
        D(CNa2) ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LCNa2+alphaNa1*CNa3+betaNa2*CNa1+alphaNa5*ICNa2-(betaNa1+alphaNa2+betaNa5+alphaNa8)*CNa2), 0),
        D(CNa1) ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LCNa1+alphaNa2*CNa2+betaNa3*ONa+alphaNa5*IFNa-(betaNa2+alphaNa3+betaNa5+alphaNa8)*CNa1), 0),
        D(ONa)  ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LONa+alphaNa3*CNa1+betaNa4*IFNa-(betaNa3+alphaNa4+alphaNa8)*ONa), 0),
        D(IFNa) ~ ifelse(INa_MarkovFlag == 1, (alphaNa4*ONa+betaNa5*CNa1+betaNa6*I1Na+alphaNa2*ICNa2-(betaNa4+alphaNa5+alphaNa6+betaNa2)*IFNa), 0),
        D(I1Na) ~ ifelse(INa_MarkovFlag == 1, (alphaNa6*IFNa+betaNa7*I2Na-(betaNa6+alphaNa7)*I1Na), 0),
        D(ICNa2)~ ifelse(INa_MarkovFlag == 1, (alphaNa1*ICNa3+betaNa2*IFNa+betaNa5*CNa2-(betaNa1+alphaNa2+alphaNa5)*ICNa2), 0),
        D(ICNa3) ~ ifelse(INa_MarkovFlag == 1, (betaNa1*ICNa2+betaNa5*CNa3-(alphaNa1+alphaNa5)*ICNa3), 0),
        D(LONa) ~ ifelse(INa_MarkovFlag == 1, (alphaNa3*LCNa1+alphaNa8*ONa-(betaNa8+betaNa3)*LONa), 0),
        D(LCNa1) ~ ifelse(INa_MarkovFlag == 1, (alphaNa8*CNa1+alphaNa2*LCNa2+betaNa3*LONa-(betaNa8+betaNa2+alphaNa3)*LCNa1), 0),
        D(LCNa2) ~ ifelse(INa_MarkovFlag == 1, (betaNa2*LCNa1+alphaNa8*CNa2+alphaNa1*LCNa3-(betaNa8+betaNa1+alphaNa2)*LCNa2), 0),
        D(LCNa3) ~ ifelse(INa_MarkovFlag == 1, (alphaNa8*CNa3+betaNa1*LCNa2-(betaNa8+alphaNa1)*LCNa3), 0)
    ]

    I_Na_junc2 = Fjunc*GNa2*(ONa+LONa)*(Vm-ena_junc)
    I_Na_sl2 = Fsl*GNa2*(ONa+LONa)*(Vm-ena_sl)

    ## I_Na: compute total current (fast and late components)
    if INa_MarkovFlag == 1
        I_Na_junc = I_Na_junc2
        I_Na_sl = I_Na_sl2
    else
        I_Na_junc = I_Na_junc1 + I_Nalj
        I_Na_sl = I_Na_sl1 + I_Nalsl
    end
    I_Na = I_Na_junc + I_Na_sl


    ## I_nabk: Na Background Current
    I_nabk_junc = Fjunc*GNaB*(Vm-ena_junc)
    I_nabk_sl = Fsl*GNaB*(Vm-ena_sl)


    ## I_nak: Na/K Pump Current
    sigma = (exp(Nao/67.3)-1)/7
    fnak = 1/(1+0.1245*exp(-0.1*Vm*FoRT)+0.0365*sigma*exp(-Vm*FoRT))
    fracPKA_PLMo = 0.116738
    fracPKA_PLMiso = 0.859251
    kPKA_PLM=KmNaip*(1-0.7019)/(fracPKA_PLMiso/fracPKA_PLMo-1)
    KmNaip_PKA=-kPKA_PLM+kPKA_PLM*(PLM_PKAp/fracPKA_PLMo)
    KmNaip = KmNaip-KmNaip_PKA
    I_nak_junc = Fjunc_nak*IbarNaK*fnak*Ko /(1+(KmNaip/Naj)^4) /(Ko+KmKo)
    I_nak_sl = Fsl_nak*IbarNaK*fnak*Ko /(1+(KmNaip/Nasl)^4) /(Ko+KmKo)
    I_nak = I_nak_junc+I_nak_sl


    ## IK_eqs
    #I_kur - IK,slow
    xurss = 1/(1+exp(-(Vm+15)/14))
    yurss = 1/(1+exp((Vm+48)/6.2))
    tauxur = 0.95+0.05*exp(-0.08*Vm)
    tauxur2 = 1+7/(1+exp(-(Vm+45)/8))+20*exp(-((Vm+35)/10)^2)
    tauyur1 = 400+900*exp(-((Vm+55)/16)^2)-250/(1+exp(-(Vm+60)/8))
    tauyur2 = 400+900*exp(-((Vm+55)/16)^2)+550/(1+exp(-(Vm+60)/8))

    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    fracIKurp0 = 0.437635  # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    fracIKurpISO = 0.718207 # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20-1)/(fracIKurpISO/fracIKurp0-1)
    fracIKuravail = (1-a_Kur)+a_Kur*(IKur_PKAp/fracIKurp0)
    I_kur1 = Kcoeff*fracIKuravail*Gkur1*IKs_x*IKs1_y*(Vm-ek)
    I_kur2 = Kcoeff*Gkur2*IKs_x*IKs2_y*(Vm-ek)
    I_kur = I_kur1 + I_kur2

    # I_ss
    xssss = xurss
    tauxss = 70*exp(-((Vm+43)/30)^2)+14
    I_ss = Kcoeff*Gss*Iss*(Vm-ek)

    # I_kr: Rapidly Activating K Current
    xrss = 1/(1+exp(-(Vm+50)/7.5))
    tauxr = 1/(1.38e-3*(Vm+7)/(1-exp(-0.123*(Vm+7)))+6.1e-4*(Vm+10)/(exp(0.145*(Vm+10))-1))
    rkr = 1/(1+exp((Vm+33)/22.4))
    I_kr = Kcoeff*gkr*Ikr*rkr*(Vm-ek)

    # I_ks: Slowly Activating K Current
    I_ks_junc = 0
    I_ks_sl = 0
    I_ks = I_ks_junc+I_ks_sl

    # I_kp: Plateau K current
    kp_kp = 1/(1+exp(7.488-Vm/5.98))
    I_kp_junc = Kcoeff*Fjunc*gkp*kp_kp*(Vm-ek)
    I_kp_sl = Kcoeff*Fsl*gkp*kp_kp*(Vm-ek)
    I_kp = I_kp_junc+I_kp_sl

    IK_eqs = [
        D(IKs_x) ~ (xurss-IKs_x)/tauxur,
        D(IKs) ~ (xurss-IKs)/tauxur2,
        D(IKs1_y) ~ (yurss-IKs1_y)/tauyur1,
        D(IKs2_y) ~ (yurss-IKs2_y)/tauyur2,
        D(Iss) ~ (xssss-Iss)/tauxss,
        D(Ikr) ~ (xrss-Ikr)/tauxr
    ]


    ## I_to: Transient Outward K Current (slow and fast components)
    # Itos (ABSENT IN MOUSE)
    xtoss = 1/(1+exp(-(Vm+3.0)/13))
    ytoss = 1/(1+exp((Vm+48)/5))
    rtoss = 1/(1+exp((Vm+33.5)/10))
    tauxtos = 0.08+0.7*exp(-((Vm+25)/30)^2)
    if ItoFlag == 0 # (not used)
        # Shannon Versions
        taurtos = 2.8e3/(1+exp((Vm+60.0)/10))+220
        tauytos = 100+400/(1+exp((Vm+25)/5))
    elseif ItoFlag == 1 && CKIIflag == 0 # WT
        # Grandi Versions
        Py = 182 
        Pr1 = 8085 
        Pr2 = 313
        taurtos = Pr1/(1+exp((Vm+33.5)/10))+Pr2
        tauytos = 100+400/(1+exp((Vm+25)/5))
    elseif ItoFlag == 1 && CKIIflag == 1
        Py = 15 
        Pr1 = 3600 
        Pr2 = 500
        taurtos = Pr1/(1+exp((Vm+33.5)/10))+Pr2
        tauytos = 100+35/(1+exp((Vm+25)/5))
    end
    I_tos = 0*Kcoeff*GtoSlow*Itos_x*Itos_y*(Vm-ek)

    # Itof
    xtofs = 1/(1+exp(-(Vm+3.0)/13))
    ytofs = 1/(1+exp((Vm+48)/5))
    tauxtof = 0.08+0.7*exp(-((Vm+25)/30)^2)
    tauytof = 10+32*exp(-((Vm+55)/16)^2)+8/(1+exp(-(Vm+60)/8))
    if CKIIflag == 1
        tauytof = 5+32*exp(-((Vm+55)/16)^2)+12/(1+exp(-(Vm+60)/8))
    end
    I_tof = Kcoeff*GtoFast*Itof_x*Itof_y*(Vm-ek) 
    I_to = I_tos + I_tof

    Ito_eqs = [
        D(Itos_x) ~ (xtoss-Itos_x)/tauxtos,
        D(Itos_y) ~ (ytoss-Itos_y)/tauytos,
        D(Itos_r)~ (rtoss-Itos_r)/taurtos,
        D(Itof_x) ~ (xtofs-Itof_x)/tauxtof,
        D(Itof_y) ~ (ytofs-Itof_y)/tauytof
    ]


    ## I_k1: Time-Independent K Current (I_ki)
    aki = 1.02/(1+exp(0.2385*(Vm-ek-59.215)))
    bki =(0.49124*exp(0.08032*(Vm+5.476-ek))+exp(0.06175*(Vm-ek-594.31)))/(1 + exp(-0.5143*(Vm-ek+4.753)))
    kiss = aki/(aki+bki)
    if CKIIflag == 1
        I_ki = 1/2*0.3*sqrt(Ko/5.4)*kiss*(Vm-ek)*Kcoeff
    else
        I_ki = 0.3*sqrt(Ko/5.4)*kiss*(Vm-ek)*Kcoeff
    end


    ## I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
    I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/Ca_j)*(Vm-ecl)
    I_ClCa_sl = Fsl*GClCa/(1+KdClCa/Ca_sl)*(Vm-ecl)
    I_ClCa = I_ClCa_junc+I_ClCa_sl
    I_Clbk = GClB*(Vm-ecl)


    ## Original H-H formulation for LCC - unused if ICa_MarkovFlag = 1
    dss = 1/(1+exp(-(Vm+14.5)/6.0))
    taud = dss*(1-exp(-(Vm+14.5)/6.0))/(0.035*(Vm+14.5))
    fss = 1/(1+exp((Vm+35.06)/3.6))+0.6/(1+exp((50-Vm)/20))
    tauf = 1/(0.0197*exp( -(0.0337*(Vm+14.5))^2 )+0.02)

    ICa_HH_eqs = [
        D(ICa_HH4) ~ 0,
        D(ICa_HH5) ~ 0,
        D(ICa_HH6) ~ 0,
        D(ICa_HH7) ~ 0
    ]

    ibarca_j = pCa*4*(Vm*Frdy*FoRT) * (0.341*Ca_j*exp(2*Vm*FoRT)-0.341*Cao) /(exp(2*Vm*FoRT)-1)
    ibarca_sl = pCa*4*(Vm*Frdy*FoRT) * (0.341*Ca_sl*exp(2*Vm*FoRT)-0.341*Cao) /(exp(2*Vm*FoRT)-1)
    ibark = pK*(Vm*Frdy*FoRT)*(0.75*Ki*exp(Vm*FoRT)-0.75*Ko) /(exp(Vm*FoRT)-1)
    ibarna_j = pNa*(Vm*Frdy*FoRT) *(0.75*Naj*exp(Vm*FoRT)-0.75*Nao)  /(exp(Vm*FoRT)-1)
    ibarna_sl = pNa*(Vm*Frdy*FoRT) *(0.75*Nasl*exp(Vm*FoRT)-0.75*Nao)  /(exp(Vm*FoRT)-1)

    I_Ca_junc1 = (Fjunc_CaL*ibarca_j*ICa_HH4*ICa_HH5*(1-ICa_HH6)*Q10CaL^Qpow)*0.45
    I_Ca_sl1 = (Fsl_CaL*ibarca_sl*ICa_HH4*ICa_HH5*(1-ICa_HH7)*Q10CaL^Qpow)*0.45
    I_CaK1 = (ibark*ICa_HH4*ICa_HH5*(Fjunc_CaL*(1-ICa_HH6)+Fsl_CaL*(1-ICa_HH7))*Q10CaL^Qpow)*0.45
    I_CaNa_junc1 = (Fjunc_CaL*ibarna_j*ICa_HH4*ICa_HH5*(1-ICa_HH6)*Q10CaL^Qpow)*0.45
    I_CaNa_sl1 = (Fsl_CaL*ibarna_sl*ICa_HH4*ICa_HH5*(1-ICa_HH7)*Q10CaL^Qpow)*.45


    ## LCC MARKOV MODEL - based on Mahajan et al. (2008)
    # LCC Current Fixed Parameters
    taupo = 1          # [ms] - Time constant of activation
    TBa = 450          # [ms] - Time constant
    s1o = 0.0221
    k1o = 0.03
    kop = 2.5e-3       # [mM]
    cpbar = 8e-3       # [mM]
    tca = 78.0312
    ICa_scale = 5.25
    recoveryReduc = 3

    # PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
    fracLCCbp0 = 0.250657
    fracLCCbpISO = 0.525870
    a_favail=(1.56-1)/(fracLCCbpISO/fracLCCbp0-1)
    favail = (1-a_favail)+a_favail*(LCCb_PKAp/fracLCCbp0)
    ICa_scale =  ICa_scale*favail
    SSAshift=0
    SSIshift=0
    # Voltage- and Ca-dependent Parameters
    poss = 1/(1+exp(-(Vm+SSAshift)/8))
    fcaj = 1/(1+(kop/Ca_j)^3) 
    Rv = 10 + 4954*exp(Vm/15.6)
    PrLCC = 1-1/(1+exp(-(Vm+40)/4))  
    PsLCC = 1/(1+exp(-(Vm+40+SSIshift)/11.32))
    TCaj = (tca + 0.1*(1+(Ca_j/cpbar)^2))/(1+(Ca_j/cpbar)^2)
    tauCaj = (Rv-TCaj)*PrLCC + TCaj
    tauBa = (Rv-TBa)*PrLCC + TBa

    # Tranisition Rates (20 rates)
    alphaLCC = poss/taupo
    betaLCC = (1-poss)/taupo
    r1 = 0.3                               # [1/ms] - Opening rate
    r2 = 3                                 # [1/ms] - closing rate
    s1 = s1o*fcaj 
    s1p = 0.00195                          # [ms] - Inactivation rate
    k1 = k1o*fcaj
    k1p = 0.00413                          # [ms] - Inactivation rate
    k2 = 1e-4                              # [ms] - Inactivation rate
    k2p = 0.00224                          # [ms] - Inactivation rate
    s2 = s1*(k2/k1)*(r1/r2)
    s2p = s1p*(k2p/k1p)*(r1/r2)
    k3 = exp(-(Vm+40)/3)/(3*(1+exp(-(Vm+40)/3)))
    k3p = k3
    k5 = (1-PsLCC)/tauCaj
    k6 = (fcaj*PsLCC)/tauCaj
    k5p = (1-PsLCC)/tauBa

    # Recovery terms
    k5 = k5/recoveryReduc
    k5p = k5p/recoveryReduc
    k6p = PsLCC/tauBa
    k4 = k3*(alphaLCC/betaLCC)*(k1/k2)*(k5/k6)
    k4p = k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p)

    # State transitions for MODE 1 junctional LCCs
    Po_LCCj_m1 = 1.0-C2_m1j-C1_m1j-I1Ca_m1j-I2Ca_m1j-I1Ba_m1j-I2Ba_m1j    # O_m1j

    ICaMar_m1j_eqs = [
        D(C2_m1j) ~ betaLCC*C1_m1j + k5*I2Ca_m1j + k5p*I2Ba_m1j - (k6+k6p+alphaLCC)*C2_m1j,                      # C2_m1j
        D(C1_m1j) ~ alphaLCC*C2_m1j + k2*I1Ca_m1j + k2p*I1Ba_m1j + r2*Po_LCCj_m1 - (r1+betaLCC+k1+k1p)*C1_m1j,   # C1_m1j
        D(I1Ca_m1j) ~ k1*C1_m1j + k4*I2Ca_m1j + s1*Po_LCCj_m1 - (k2+k3+s2)*I1Ca_m1j,                              # I1Ca_m1j
        D(I2Ca_m1j) ~ k3*I1Ca_m1j + k6*C2_m1j - (k4+k5)*I2Ca_m1j,                                                 # I2Ca_m1j
        D(I1Ba_m1j) ~ k1p*C1_m1j + k4p*I2Ba_m1j + s1p*Po_LCCj_m1 - (k2p+k3p+s2p)*I1Ba_m1j,                        # I1Ba_m1j
        D(I2Ba_m1j) ~ k3p*I1Ba_m1j + k6p*C2_m1j - (k5p+k4p)*I2Ba_m1j                                             # I2Ba_m1j
    ]

    ibarca_jm1 = (4*pCa*Vm*Frdy*FoRT)*(0.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1)
    I_Ca_junc_m1 = (Fjunc_CaL*ibarca_jm1*Po_LCCj_m1*Q10CaL^Qpow)*ICa_scale

    # Re-define all parameters as mode 2 specific parameters
    s1om2 = 0.0221
    k1om2 = 0.03
    kopm2 = 2.5e-3
    cpbarm2 = 8e-3
    tcam2 = 78.0312
    possm2 = 1/(1+exp(-(Vm+SSAshift)/8))
    fcajm2 = 1/(1+(kopm2/Ca_j)^3)
    Rvm2 = 10 + 4954*exp(Vm/15.6)
    PrLCCm2 = 1-1/(1+exp(-(Vm+40)/4))
    PsLCCm2 = 1/(1+exp(-(Vm+40+SSIshift)/11.32))
    TCajm2 = (tcam2 + 0.1*(1+(Ca_j/cpbarm2)^2))/(1+(Ca_j/cpbarm2)^2)
    tauCajm2 = (Rvm2-TCajm2)*PrLCCm2 + TCajm2
    tauBam2 = (Rvm2-TBa)*PrLCCm2 + TBa
    alphaLCCm2 = possm2/taupo
    betaLCCm2 = (1-possm2)/taupo
    r1m2 = 0.3                               # [1/ms] - Opening rate
    r2m2 = 3/8                               # [1/ms] - closing rate
    s1m2 = s1om2*fcajm2 
    s1pm2 = 0.00195                           # [ms] - Inactivation rate
    k1m2 = k1om2*fcajm2 
    k1pm2 = 0.00413                           # [ms] - Inactivation rate
    k2m2 = 1e-4                              # [ms] - Inactivation rate
    k2pm2 = 0.00224                           # [ms] - Inactivation rate
    s2m2 = s1m2*(k2m2/k1m2)*(r1m2/r2m2)
    s2pm2 = s1pm2*(k2pm2/k1pm2)*(r1m2/r2m2)
    k3m2 = exp(-(Vm+40)/3)/(3*(1+exp(-(Vm+40)/3)))
    k3pm2 = k3m2
    k5m2 = (1-PsLCCm2)/tauCajm2
    k6m2 = (fcajm2*PsLCCm2)/tauCajm2
    k5pm2 = (1-PsLCCm2)/tauBam2
    k5m2 = k5m2/recoveryReduc      # reduced for recovery
    k5pm2 = k5pm2/recoveryReduc    # reduced for recovery    
    k6pm2 = PsLCCm2/tauBam2
    k4m2 = k3m2*(alphaLCCm2/betaLCCm2)*(k1m2/k2m2)*(k5m2/k6m2)
    k4pm2 = k3pm2*(alphaLCCm2/betaLCCm2)*(k1pm2/k2pm2)*(k5pm2/k6pm2)

    # State transitions for MODE 2 junctional LCCs
    Po_LCCj_m2 = 1.0-C2_m2j-C1_m2j-I1Ca_m2j-I2Ca_m2j-I1Ba_m2j-I2Ba_m2j                                                           # O_m2j

    ICaMar_m2j_eqs=[
        D(C2_m2j) ~ betaLCCm2*C1_m2j + k5m2*I2Ca_m2j + k5pm2*I2Ba_m2j - (k6m2+k6pm2+alphaLCCm2)*C2_m2j,                          # C2_m2j
        D(C1_m2j) ~ alphaLCCm2*C2_m2j + k2m2*I1Ca_m2j + k2pm2*I1Ba_m2j + r2m2*Po_LCCj_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*C1_m2j,   # C1_m2j
        D(I1Ca_m2j) ~ k1m2*C1_m2j + k4m2*I2Ca_m2j + s1m2*Po_LCCj_m2 - (k2m2+k3m2+s2m2)*I1Ca_m2j,                                  # I1Ca_m2j
        D(I2Ca_m2j) ~ k3m2*I1Ca_m2j + k6m2*C2_m2j - (k4m2+k5m2)*I2Ca_m2j,                                                         # I2Ca_m2j
        D(I1Ba_m2j) ~ k1pm2*C1_m2j + k4pm2*I2Ba_m2j + s1pm2*Po_LCCj_m2 - (k2pm2+k3pm2+s2pm2)*I1Ba_m2j,                            # I1Ba_m2j
        D(I2Ba_m2j) ~ k3pm2*I1Ba_m2j + k6pm2*C2_m2j - (k5pm2+k4pm2)*I2Ba_m2j                                                     # I2Ba_m2j
    ]

    ibarca_jm2 = (4*pCa*Vm*Frdy*FoRT)*(.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1)
    I_Ca_junc_m2 = (Fjunc_CaL*ibarca_jm2*(Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale

    # CaMKII AND PKA-DEPENDENT SHIFTING OF DYADIC LCCS TO MODE 2
    fracLCCap0 = 0.219577
    frac_fpkam2 = (0.15*fracLCCap0)/(1-fracLCCap0)
    fpkam2 = (0.15+frac_fpkam2)*LCCa_PKAp - frac_fpkam2
    fckiim2 = LCC_CKp*0.1 # Assumes max phosphorylation results in 10# mode 2 channels (max LCC_CKp = 1)
    # Sum up total fraction of CKII and PKA-shifted mode 2 channels
    junc_mode2 = fckiim2 + fpkam2
    # Total junctional ICa
    I_Ca_junc2 = (1-junc_mode2)*I_Ca_junc_m1 + junc_mode2*I_Ca_junc_m2
    # SUB-SARCOLEMMAL LCCs
    # Re-assign necessary params to be Casl sensitive
    fcasl = 1/(1+(kop/Ca_sl)^3)
    TCasl = (tca + 0.1*(1+(Ca_sl/cpbar))^2)/(1+(Ca_sl/cpbar)^2)
    tauCasl = (Rv-TCasl)*PrLCC + TCasl
    # Re-assign necessary rates to be Casl sensitive
    s1sl = s1o*fcasl
    k1sl = k1o*fcasl
    s2sl = s1sl*(k2/k1sl)*(r1/r2)
    s2psl = s1p*(k2p/k1p)*(r1/r2)
    k5sl = (1-PsLCC)/tauCasl
    k5sl = k5sl/recoveryReduc  # Reduced for recovery
    k6sl = (fcasl*PsLCC)/tauCasl
    k4sl = k3*(alphaLCC/betaLCC)*(k1sl/k2)*(k5sl/k6sl)
    k4psl = k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p)

    # State transitions for 'mode 1' sarcolemmal LCCs
    Po_LCCsl_m1 = 1-C2_m1sl-C1_m1sl-I1Ca_m1sl-I2Ca_m1sl-I1Ba_m1sl-I2Ba_m1sl;                                                # O_m1sl
    ICaMar_m1sl_eqs = [
        D(C2_m1sl) ~ betaLCC*C1_m1sl + k5sl*I2Ca_m1sl + k5p*I2Ba_m1sl - (k6sl+k6p+alphaLCC)*C2_m1sl,                      # C2_m1sl
        D(C1_m1sl) ~ alphaLCC*C2_m1sl + k2*I1Ca_m1sl + k2p*I1Ba_m1sl + r2*Po_LCCsl_m1 - (r1+betaLCC+k1sl+k1p)*C1_m1sl,    # C1_m1sl
        D(I1Ca_m1sl) ~ k1sl*C1_m1sl + k4sl*I2Ca_m1sl + s1sl*Po_LCCsl_m1 - (k2+k3+s2sl)*I1Ca_m1sl,                         # I1Ca_m1sl
        D(I2Ca_m1sl) ~ k3*I1Ca_m1sl + k6sl*C2_m1sl - (k4sl+k5sl)*I2Ca_m1sl,                                               # I2Ca_m1sl
        D(I1Ba_m1sl) ~ k1p*C1_m1sl + k4psl*I2Ba_m1sl + s1p*Po_LCCsl_m1 - (k2p+k3p+s2psl)*I1Ba_m1sl,                       # I1Ba_m1sl
        D(I2Ba_m1sl) ~ k3p*I1Ba_m1sl + k6p*C2_m1sl - (k5p+k4psl)*I2Ba_m1sl                                               # I2Ba_m1sl
    ]
    ibarca_slm1 = (4*pCa*Vm*Frdy*FoRT)*(.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1)
    I_Casl_m1 = (Fsl_CaL*ibarca_slm1*Po_LCCsl_m1*Q10CaL^Qpow)*ICa_scale

    # Adjust closing rate for 'mode 2' sarcolemmal LCCs
    r2slm2 = r2m2
    s2slm2 = s1sl*(k2/k1sl)*(r1/r2slm2)
    s2pslm2 = s1p*(k2p/k1p)*(r1/r2slm2)
    # State transitions for mode 2 sarcolemmal LCCs
    Po_LCCsl_m2 = 1-C2_m2sl-C1_m2sl-I1Ca_m2sl-I2Ca_m2sl-I1Ba_m2sl-I2Ba_m2sl     # O_m2sl
    ICaMar_m2sl_eqs = [
        D(C2_m2sl) ~ betaLCC*C1_m2sl + k5sl*I2Ca_m2sl + k5p*I2Ba_m2sl - (k6sl+k6p+alphaLCC)*C2_m2sl,                        # C2_m2sl
        D(C1_m2sl) ~ alphaLCC*C2_m2sl + k2*I1Ca_m2sl + k2p*I1Ba_m2sl + r2slm2*Po_LCCsl_m2 - (r1+betaLCC+k1sl+k1p)*C1_m2sl,  # C1_m2sl
        D(I1Ca_m2sl) ~ k1sl*C1_m2sl + k4sl*I2Ca_m2sl + s1sl*Po_LCCsl_m2 - (k2+k3+s2slm2)*I1Ca_m2sl,                       # I1Ca_m2sl
        D(I2Ca_m2sl) ~ k3*I1Ca_m2sl + k6sl*C2_m2sl - (k4sl+k5sl)*I2Ca_m2sl,                                               # I2Ca_m2sl
        D(I1Ba_m2sl) ~ k1p*C1_m2sl + k4psl*I2Ba_m2sl + s1p*Po_LCCsl_m2 - (k2p+k3p+s2pslm2)*I1Ba_m2sl,                     # I1Ba_m2sl
        D(I2Ba_m2sl) ~ k3p*I1Ba_m2sl + k6p*C2_m2sl - (k5p+k4psl)*I2Ba_m2sl                                                # I2Ba_m2sl
    ]

    ibarca_slm2 = (4*pCa*Vm*Frdy*FoRT)*(.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1)
    I_Casl_m2 = (Fsl_CaL*ibarca_slm2*Po_LCCsl_m2*Q10CaL^Qpow)*ICa_scale
    # Sum mode 1 and mode 2 sl channels for total sl current
    fckiim2_sl = 0
    sl_mode2 = fckiim2_sl + fpkam2
    I_Ca_sl2 = (1-sl_mode2)*I_Casl_m1 + sl_mode2*I_Casl_m2
    # Na and K currents through LCC
    I_CaKj2 = ibark*Fjunc_CaL*((1-junc_mode2)*Po_LCCj_m1 + junc_mode2*Po_LCCj_m2)*Q10CaL^Qpow*ICa_scale
    I_CaKsl2 = ibark*Fsl_CaL*((1-sl_mode2)*Po_LCCsl_m1 + sl_mode2*Po_LCCsl_m2)*Q10CaL^Qpow*ICa_scale
    I_CaK2 = I_CaKj2+I_CaKsl2
    I_CaNa_junc2 = (Fjunc_CaL*ibarna_j*((1-junc_mode2)*Po_LCCj_m1+junc_mode2*Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale
    I_CaNa_sl2 = Fsl_CaL*ibarna_sl*((1-sl_mode2)*Po_LCCsl_m1 + sl_mode2*Po_LCCsl_m2)*Q10CaL^Qpow*ICa_scale
    # These are now able to switch depending on whether or not the flag to
    # switch to Markov model of ICa is ON
    I_Ca_junc = (1-ICa_MarkovFlag)*I_Ca_junc1 + ICa_MarkovFlag*I_Ca_junc2
    I_Ca_sl = (1-ICa_MarkovFlag)*I_Ca_sl1 + ICa_MarkovFlag*I_Ca_sl2
    I_Ca = I_Ca_junc+I_Ca_sl
    I_CaNa_junc = (1-ICa_MarkovFlag)*(I_CaNa_junc1) + (ICa_MarkovFlag)*(I_CaNa_junc2)
    I_CaNa_sl = (1-ICa_MarkovFlag)*(I_CaNa_sl1) + (ICa_MarkovFlag)*(I_CaNa_sl2)
    I_CaNa = I_CaNa_junc + I_CaNa_sl
    I_CaK = (1-ICa_MarkovFlag)*(I_CaK1) + ICa_MarkovFlag*(I_CaK2)
    # Collect all currents through LCC
    I_Catot = I_Ca+I_CaK+I_CaNa
    ICa_total = [
        D(influx_LTCC) ~ -I_Ca*Cmem/(Vmyo*2*Frdy)*1e3
    ]


    ## I_ncx: Na/Ca Exchanger flux
    Ka_junc = 1/(1+(Kdact/Ca_j)^3)
    Ka_sl = 1/(1+(Kdact/Ca_sl)^3)
    s1_junc = exp(nu*Vm*FoRT)*Naj^3*Cao
    s1_sl = exp(nu*Vm*FoRT)*Nasl^3*Cao
    s2_junc = exp((nu-1)*Vm*FoRT)*Nao^3*Ca_j
    s3_junc = (KmCai*Nao^3*(1+(Naj/KmNai)^3)+KmNao^3*Ca_j+ KmNai^3*Cao*(1+Ca_j/KmCai)+KmCao*Naj^3+Naj^3*Cao+Nao^3*Ca_j)*(1+ksat*exp((nu-1)*Vm*FoRT))
    s2_sl = exp((nu-1)*Vm*FoRT)*Nao^3*Ca_sl
    s3_sl = (KmCai*Nao^3*(1+(Nasl/KmNai)^3) + KmNao^3*Ca_sl+KmNai^3*Cao*(1+Ca_sl/KmCai)+KmCao*Nasl^3+Nasl^3*Cao+Nao^3*Ca_sl)*(1+ksat*exp((nu-1)*Vm*FoRT))
    I_ncx_junc = Fjunc_ncx*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc
    I_ncx_sl = Fsl_ncx*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl
    I_ncx = I_ncx_junc+I_ncx_sl
    I_ncx_eqs = [
        D(influx_NCX) ~ 2*I_ncx*Cmem/(Vmyo*2*Frdy)*1e3 #uM/ms
    ]


    ## I_pca: Sarcolemmal Ca Pump Current
    I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*Ca_j^1.6/(KmPCa^1.6+Ca_j^1.6)
    I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*Ca_sl^1.6/(KmPCa^1.6+Ca_sl^1.6)
    I_pca = I_pca_junc+I_pca_sl
    Ipca_eqs = [
        D(influx_PMCA) ~ -I_pca*Cmem/(Vmyo*2*Frdy)*1e3
    ]


    ## I_cabk: Ca Background Current
    I_cabk_junc = Fjunc*GCaB*(Vm-eca_junc)
    I_cabk_sl = Fsl*GCaB*(Vm-eca_sl)
    I_cabk = I_cabk_junc+I_cabk_sl
    ICabp_eqs = [
        D(influx_ICa) ~ -I_cabk*Cmem/(Vmyo*2*Frdy)*1e3
    ]


    ## I_CFTR or I_cl_(cAMP) - Cystic Fibrosis Transmembrane Conductance Reg.
    Icftr = 0   # NO Icftr in MOUSE


    ## RyR model - SR release fluxes and leak
    # CaMKII and PKA-dependent phosphoregulation of RyR Po
    fCKII_ec50SR = 1.16 - 4/5*RyR_CKp
    ec50SR = fCKII_ec50SR*ec50SR
    MaxSR = 15
    MinSR = 1
    kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/Ca_sr)^2.5)
    koSRCa = koCa/kCaSR
    kiSRCa = kiCa*kCaSR
    kleak = 2*5.348e-6
    fCKII_RyR = (10*RyR_CKp - 1)
    frac_RyRo = 0.204276
    a_RyR = (2-1)/(1/frac_RyRo-1)
    fPKA_RyR = 1-a_RyR+a_RyR*(RyR_PKAp/frac_RyRo)
    koSRCa = (fCKII_RyR + fPKA_RyR - 1)*koSRCa

    # ODEs for RyR states and SR release through open RyRs
    RI = 1-RyR_R-RyR_O-RyR_I
    RyR_eqs = [
        D(RyR_R) ~ (kim*RI-kiSRCa*Ca_j*RyR_R)-(koSRCa*Ca_j^2*RyR_R-kom*RyR_O),
        D(RyR_O) ~ (koSRCa*Ca_j^2*RyR_R-kom*RyR_O)-(kiSRCa*Ca_j*RyR_O-kim*RyR_I),
        D(RyR_I) ~ (kiSRCa*Ca_j*RyR_O-kim*RyR_I)-(kom*RyR_I-koSRCa*Ca_j^2*RI)
    ]
    J_SRCarel = ks*RyR_O*(Ca_sr-Ca_j)
    # Passive RyR leak - includes CaMKII regulation of leak flux
    kleak = (1/2 + 5*RyR_CKp/2)*kleak
    J_SRleak = kleak*(Ca_sr-Ca_j)


    ## SERCA model - SR uptake fluxes
    # CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
    fCKII_PLB = (1-0.5*PLB_CKp)
    fracPKA_PLBo = 1-0.079755
    fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*(100-55.31)/100 + 55.31/100

    # Select smaller value (resulting in max reduction of Kmf)
    Kmf = ifelse(fCKII_PLB < fPKA_PLB, Kmf*fCKII_PLB, Kmf*fPKA_PLB)     #fCKII_PLB
    J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((Cai/Kmf)^hillSRCaP-(Ca_sr/Kmr)^hillSRCaP)/(1+(Cai/Kmf)^hillSRCaP+(Ca_sr/Kmr)^hillSRCaP)


    ## Na and Ca Buffering
    Buffering_eqs = [
        D(NaBj) ~ kon_na*Naj*(Bmax_Naj-NaBj)-koff_na*NaBj,        # NaBj      [mM/ms]
        D(NaBsl) ~ kon_na*Nasl*(Bmax_Nasl-NaBsl)-koff_na*NaBsl,       # NaBsl     [mM/ms]
        D(TnCL) ~ kon_tncl*Cai*(Bmax_TnClow-TnCL)-koff_tncl*TnCL,            # TnCL      [mM/ms]
        D(TnCHc) ~ kon_tnchca*Cai*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchca*TnCHc, # TnCHc     [mM/ms]
        D(TnCHm) ~ kon_tnchmg*Mgi*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchmg*TnCHm,   # TnCHm     [mM/ms]
        D(CaM) ~ 0,
        D(Myosin_ca) ~ kon_myoca*Cai*(Bmax_myosin-Myosin_ca-Myosin_mg)-koff_myoca*Myosin_ca,    # Myosin_ca [mM/ms]
        D(Myosin_mg) ~ kon_myomg*Mgi*(Bmax_myosin-Myosin_ca-Myosin_mg)-koff_myomg*Myosin_mg,      # Myosin_mg [mM/ms]
        D(SRB) ~ kon_sr*Cai*(Bmax_SR-SRB)-koff_sr*SRB,                    # SRB       [mM/ms]
        D(SLLj) ~ kon_sll*Ca_j*(Bmax_SLlowj-SLLj)-koff_sll*SLLj,       # SLLj      [mM/ms]
        D(SLLsl) ~ kon_sll*Ca_sl*(Bmax_SLlowsl-SLLsl)-koff_sll*SLLsl,      # SLLsl     [mM/ms]
        D(SLHj) ~ kon_slh*Ca_j*(Bmax_SLhighj-SLHj)-koff_slh*SLHj,      # SLHj      [mM/ms]
        D(SLHsl) ~ kon_slh*Ca_sl*(Bmax_SLhighsl-SLHsl)-koff_slh*SLHsl      # SLHsl     [mM/ms]
    ]
    # Cytosolic Ca Buffers
    J_CaB_cytosol = D(TnCL)+D(TnCHc)+D(CaM)+D(Myosin_ca)+D(SRB)
    # Junctional and SL Ca Buffers
    J_CaB_junction = D(SLLj)+D(SLHj)
    J_CaB_sl = D(SLLsl)+D(SLHsl)


    ## Ion concentrations
    # Na Concentrations
    I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc   # [uA/uF]
    I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl   #[uA/uF]
    # K Concentration
    I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur+I_ss     # [uA/uF]
    # Ca Concentrations
    I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc # [uA/uF]
    I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl           # [uA/uF]

    Ion_eqs = [
        D(Csqn) ~ (kon_csqn*Ca_sr*(Bmax_Csqn-Csqn)-koff_csqn*Csqn),
        D(Ca_sr) ~ (J_serca*Vmyo/Vsr-(J_SRleak*Vmyo/Vsr+J_SRCarel))-D(Csqn),
        D(Naj) ~ -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(Nasl-Naj)-D(NaBj),
        D(Nasl) ~ (-I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(Naj-Nasl)+J_na_slmyo/Vsl*(Nai-Nasl))-D(NaBsl),
        D(Nai) ~ J_na_slmyo/Vmyo*(Nasl-Nai),
        D(Ki) ~ D(Ca_j)+ 1e-3*JCaDyad,
        D(Ca_j) ~ (-I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(Ca_sl-Ca_j)-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc),
        D(Ca_sl) ~ (-I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(Ca_j-Ca_sl)+ J_ca_slmyo/Vsl*(Cai-Ca_sl)-J_CaB_sl+ 1e-3*JCaSL),
        D(Cai) ~ (-J_serca-J_CaB_cytosol+J_ca_slmyo/Vmyo*(Ca_sl-Cai)+1e-3*JCaCyt)
    ]


    ## Simulation type
    protocol = 1        # "pace"
    if CaffeineFlag==1
        protocol = 2    # "vcRest"
    end
    if StrophFlag == 1
        protocol = 3    # "none"
    end
    # AP Waveform for AP clamp
    if protocol == 1
        I_app = ifelse(mod(t,cycleLength) <= 5, 9.5, 0.0)
    elseif protocol == 2 
        V_clamp = -83   #Vm#-80
        R_clamp = 0.01
        I_app = (V_clamp-Vm)/R_clamp
    else
        I_app = 0.0
    end  

    ## Membrane Potential
    I_Na_tot = I_Na_tot_junc + I_Na_tot_sl                 # [uA/uF]
    I_Cl_tot = I_ClCa+I_Clbk+Icftr                         # [uA/uF]
    I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl                   # [uA/uF]
    I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot             # [uA/uF]

    MemPot_eqs = [
        D(Vm) ~ -(I_tot-I_app)
    ]


    return vcat(INa_fast_eqs, INa_eqs, IK_eqs, Ito_eqs, ICa_HH_eqs, ICaMar_m1j_eqs, ICa_total, Na_h,
                ICaMar_m2j_eqs, ICaMar_m1sl_eqs, ICaMar_m2sl_eqs, I_ncx_eqs, Ipca_eqs, ICabp_eqs, RyR_eqs, 
                Buffering_eqs, Ion_eqs, MemPot_eqs, CaMdyad_eqs, CaMSL_eqs, CaMcyt_eqs, CaMKII_eqs, bar_eqs, 
                cAMP_eqs, PKA_eqs, PP1_eqs, PLB_eqs, PLM_eqs, LCC_eqs, RyRp_eqs, TnI_eqs, IKs_eqs, CFTR_eqs, Ikur_eqs)
end



eq_morotti = get_Morotti_equations()

@named osys = ODESystem(eq_morotti)

osys = structural_simplify(osys)

@variables t Cai(t)

##Chemical Reaction
ca_model = @reaction_network begin
    ##(d*50e-9, d), 0 <--> Ca
    ##  Two Ca2+ ions bind to C or N-lobe.
    (k_1C_on*($Cai)^2*k_2C_on/(k_1C_off+k_2C_on*($Cai)),k_1C_off*k_2C_off/(k_1C_off+k_2C_on*($Cai))), CaM0 <--> Ca2CaM_C
    (k_1N_on*($Cai)^2*k_2N_on/(k_1N_off+k_2N_on*($Cai)), k_1N_off*k_2N_off/(k_1N_off+k_2N_on*($Cai))), CaM0 <--> Ca2CaM_N
    (k_1C_on*($Cai)^2*k_2C_on/(k_1C_off+k_2C_on*($Cai)), k_1C_off*k_2C_off/(k_1C_off+k_2C_on*($Cai))), Ca2CaM_C <--> Ca4CaM
    (k_1N_on*($Cai)^2*k_2N_on/(k_1N_off+k_2N_on*($Cai)), k_1N_off*k_2N_off/(k_1N_off+k_2N_on*($Cai))), Ca2CaM_N <--> Ca4CaM
    ##  Two Ca2+ ions bind to C or N-lobe of CaM-CaMKII complex.
    (k_K1C_on*($Cai)^2*k_K2C_on/(k_K1C_off+k_K2C_on*($Cai)), k_K1C_off*k_K2C_off/(k_K1C_off+k_K2C_on*($Cai))), CaM0_CaMK <--> Ca2CaM_C_CaMK
    (k_K1N_on*($Cai)^2*k_K2N_on/(k_K1N_off+k_K2N_on*($Cai)), k_K1N_off*k_K2N_off/(k_K1N_off+k_K2N_on*($Cai))), CaM0_CaMK <--> Ca2CaM_N_CaMK
    (k_K1C_on*($Cai)^2*k_K2C_on/(k_K1C_off+k_K2C_on*($Cai)), k_K1C_off*k_K2C_off/(k_K1C_off+k_K2C_on*($Cai))), Ca2CaM_C_CaMK <--> Ca4CaM_CaMK
    (k_K1N_on*($Cai)^2*k_K2N_on/(k_K1N_off+k_K2N_on*($Cai)), k_K1N_off*k_K2N_off/(k_K1N_off+k_K2N_on*($Cai))), Ca2CaM_N_CaMK <--> Ca4CaM_CaMK
    ##  Binding of Ca to CaM-CaMKIIP.
    (k_K1C_on*k_K2C_on/(k_K1C_off+k_K2C_on*($Cai))*($Cai)^2, k_K1C_off*k_K2C_off/(k_K1C_off+k_K2C_on*($Cai))), CaM0_CaMKP <--> Ca2CaM_C_CaMKP
    (k_K1N_on*k_K2N_on/(k_K1N_off+k_K2N_on*($Cai))*($Cai)^2, k_K1N_off*k_K2N_off/(k_K1N_off+k_K2N_on*($Cai))), CaM0_CaMKP <--> Ca2CaM_N_CaMKP
    (k_K1C_on*k_K2C_on/(k_K1C_off+k_K2C_on*($Cai))*($Cai)^2, k_K1C_off*k_K2C_off/(k_K1C_off+k_K2C_on*($Cai))), Ca2CaM_C_CaMKP <--> Ca4CaM_CaMKP
    (k_K1N_on*k_K2N_on/(k_K1N_off+k_K2N_on*($Cai))*($Cai)^2, k_K1N_off*k_K2N_off/(k_K1N_off+k_K2N_on*($Cai))), Ca2CaM_N_CaMKP <--> Ca4CaM_CaMKP
    ##  Binding of CaM to CaMKII or CaMII-P
    (kCaM0_on, kCaM0_off), CaM0 + CaMK <--> CaM0_CaMK
    (kCaM2C_on, kCaM2C_off), Ca2CaM_C + CaMK <--> Ca2CaM_C_CaMK
    (kCaM2N_on, kCaM2N_off), Ca2CaM_N + CaMK <--> Ca2CaM_N_CaMK
    (kCaM4_on, kCaM4_off), Ca4CaM + CaMK <--> Ca4CaM_CaMK
    (kCaM0P_on, kCaM0P_off), CaM0 + CaMKP <--> CaM0_CaMKP
    (kCaM2CP_on, kCaM2CP_off), Ca2CaM_C + CaMKP <--> Ca2CaM_C_CaMKP
    (kCaM2NP_on, kCaM2NP_off), Ca2CaM_N + CaMKP <--> Ca2CaM_N_CaMKP
    (kCaM4P_on, kCaM4P_off), Ca4CaM + CaMKP <--> Ca4CaM_CaMKP
    ##  Phosphorylation CaMXCaMKII -> CaMXCaMKIIP. 
    k_phosCaM*(CaMKP+CaMKP2+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)/CaMKII_T, Ca2CaM_C_CaMK --> Ca2CaM_C_CaMKP
    k_phosCaM*(CaMKP+CaMKP2+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)/CaMKII_T, Ca2CaM_N_CaMK --> Ca2CaM_N_CaMKP
    k_phosCaM*(CaMKP+CaMKP2+CaM0_CaMK+Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP)/CaMKII_T, Ca4CaM_CaMK --> Ca4CaM_CaMKP
    ##  Dephosphorylation CaMKP -> CaMK
    k_dephospho, CaMKP --> CaMK
    ##  Second phosphorylation state (P2) CaMKP <-> CaMKP2
    (k_P1_P2, k_P2_P1), CaMKP <--> CaMKP2
end

###########################  Parameters  ###########################
CaMT = 30e-6 #Total calmodulin concentration.
CaMKII_T = 70e-6 #Total CaMKII concentration.

binding_To_PCaMK = 0.1
decay_CaM = 3 # seconds
phospho_rate = 1
phosphatase = 1

rn_osys = convert(ODESystem, ca_model)
@named sys = extend(osys, rn_osys)


@unpack Na_m, Na_h, Na_j, ICa_HH4, ICa_HH5, ICa_HH6, ICa_HH7, Itos_x, Itos_y, Itof_x, Itof_y, Ikr, IKs, RyR_R, RyR_O, RyR_I, NaBj, NaBsl,  
        TnCL, TnCHc, TnCHm, CaM, Myosin_ca, Myosin_mg, SRB, SLLj, SLLsl, SLHj, SLHsl, Csqn, Ca_sr, Naj, Nasl, Nai, Ki, Ca_j, Ca_sl, Cai, Vm, 
        Itos_r, influx_LTCC, influx_PMCA, influx_NCX, influx_ICa, Na_late_h, CNa2, CNa1, ONa, IFNa, I1Na, CNa3, ICNa2, ICNa3, LONa, LCNa1, LCNa2, LCNa3, 
        C2_m1j, C1_m1j, I1Ca_m1j, I2Ca_m1j, I1Ba_m1j, I2Ba_m1j, C2_m2j, C1_m2j, I1Ca_m2j, I2Ca_m2j, I1Ba_m2j, I2Ba_m2j, C2_m1sl, C1_m1sl, I1Ca_m1sl, 
        I2Ca_m1sl, I1Ba_m1sl, I2Ba_m1sl, C2_m2sl, C1_m2sl, I1Ca_m2sl, I2Ca_m2sl, I1Ba_m2sl, I2Ba_m2sl, IKs_x, IKs1_y, Iss, IKs2_y,  # ecc_ODEfile
        CaM_dyad, Ca2CaM_dyad, Ca4CaM_dyad, CaMB_dyad, Ca2CaMB_dyad, Ca4CaMB_dyad, Pb2_dyad, Pb_dyad, 
        Pt_dyad, Pt2_dyad, Pa_dyad, Ca4CaN_dyad, CaMCa4CaN_dyad, Ca2CaMCa4CaN_dyad, Ca4CaMCa4CaN_dyad,                              # camdyad_ODEfile
        CaM_sl, Ca2CaM_sl, Ca4CaM_sl, CaMB_sl, Ca2CaMB_sl, Ca4CaMB_sl, Pb2_sl, Pb_sl, 
        Pt_sl, Pt2_sl, Pa_sl, Ca4CaN_sl, CaMCa4CaN_sl, Ca2CaMCa4CaN_sl, Ca4CaMCa4CaN_sl,                                            # camsl_ODEfile
        CaM_cyt, Ca2CaM_cyt, Ca4CaM_cyt, CaMB_cyt, Ca2CaMB_cyt, Ca4CaMB_cyt, Pb2_cyt, Pb_cyt, 
        Pt_cyt, Pt2_cyt, Pa_cyt, Ca4CaN_cyt, CaMCa4CaN_cyt, Ca2CaMCa4CaN_cyt, Ca4CaMCa4CaN_cyt,                                     # camcyt_ODEfile
        LCC_PKAp, LCC_CKdyadp, RyR2809p, RyR2815p, PLBT17p, LCC_CKslp,                                                              # camkii_ODEfile
        LR, LRG, RG, b1AR_S464, b1AR_S301, GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, PDEp, cAMPtot, RC_I, RCcAMP_I, 
        RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII,                         # bar_ODEfile
        PKACII_PKI, I1p_PP1, I1ptot, PLBp, PLMp, LCCap, LCCbp, RyRp, TnIp, KS79, KS80, KSp, CFTRp, KURp,
        # CaMKII Model
        CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, 
        Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2, k_1C_on, k_1C_off, k_2C_on, k_2C_off,
        k_1N_on, k_1N_off, k_2N_on, k_2N_off, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off, k_K1N_on, k_K1N_off,
        k_K2N_on, k_K2N_off, kCaM0_on, kCaM2C_on, kCaM2N_on, kCaM4_on, kCaM0_off, kCaM2C_off, kCaM2N_off, kCaM4_off, 
        kCaM0P_on, kCaM2CP_on, kCaM2NP_on, kCaM4P_on, kCaM0P_off, kCaM2CP_off, kCaM2NP_off, kCaM4P_off, k_phosCaM, 
        k_dephospho, k_P1_P2, k_P2_P1, CaMKII_T = sys

oprob = ODEProblem(sys, [
        Na_m => 1.94e-3, Na_h => 0.981, Na_j => 0.987, ICa_HH4 => 7.02e-6, ICa_HH5 => 1.00068, ICa_HH6 => 2.7e-2, 
        ICa_HH7 => 1.6e-2, Itos_x => 2.02e-3, Itos_y => 0.99, Itof_x => 2.02e-3, Itof_y => 0.9992, Ikr => 1.11e-2, 
        IKs => 7.37e-3, RyR_R => 0.698, RyR_O => 4.24e-6, RyR_I => 1.84e-6, NaBj => 3.993, NaBsl => 0.87, TnCL => 9.26e-3, 
        TnCHc => 0.118, TnCHm => 1.03e-2, CaM => 2.53e-4, Myosin_ca => 1.989e-3, Myosin_mg => 0.138, SRB => 2.26e-3, 
        SLLj => 2.2e-2, SLLsl => 1.35e-2, SLHj => 0.127, SLHsl => 0.142, Csqn => 1.177, Ca_sr => 0.503, Naj => 11.182, 
        Nasl => 11.182, Nai => 11.182, Ki => 134.99, Ca_j => 5.34e-4, Ca_sl => 1.46e-4, Cai => 9.12e-5, Vm => -83.632, 
        Itos_r => 0.946, influx_LTCC => 5.59e4, influx_PMCA => -3.38e4, influx_NCX => -3.096e5, influx_ICa => 2.875e5, 
        Na_late_h => 0.222, CNa2 => 0.105, CNa1 => 1.92e-3, ONa => 4.15e-5, IFNa => 0.303, I1Na => 0.566, CNa3 => 1.01e-2, 
        ICNa2 => 7.01e-5, ICNa3 => 8.62e-8, LONa => 1.47e-4, LCNa1 => 2.64e-6, LCNa2 => 1.82e-8, LCNa3 => 2.17e-11, 
        C2_m1j => 0.939, C1_m1j => 2.71e-5, I1Ca_m1j => 9.17e-5, I2Ca_m1j => 6.71e-4, I1Ba_m1j => 4.99e-5, I2Ba_m1j => 5.97e-2,
        C2_m2j => 0.939, C1_m2j => 2.71e-5, I1Ca_m2j => 9.13e-5, I2Ca_m2j => 6.69e-4, I1Ba_m2j => 4.99e-5, I2Ba_m2j => 5.97e-2, 
        C2_m1sl => 0.94, C1_m1sl => 2.71e-5, I1Ca_m1sl => 3.88e-6, I2Ca_m1sl => 2.81e-5, I1Ba_m1sl => 5.00e-5, I2Ba_m1sl => 5.98e-2,
        C2_m2sl => 0.94, C1_m2sl => 2.71e-5, I1Ca_m2sl => 4.17e-6, I2Ca_m2sl => 3.02e-5, I1Ba_m2sl => 5.00e-5, I2Ba_m2sl => 5.98e-2, 
        IKs_x => 7.37e-3, IKs1_y => 0.99, Iss => 7.37e-3, IKs2_y => 0.995, CaM_dyad => 388.68, Ca2CaM_dyad => 13.02, 
        Ca4CaM_dyad => 8.68e-3, CaMB_dyad => 0.0, Ca2CaMB_dyad => 0.0, Ca4CaMB_dyad => 0.0, Pb2_dyad => 0.67, Pb_dyad => 7.09e-2, 
        Pt_dyad => 2.42e-5, Pt2_dyad => 8.73e-9, Pa_dyad => 3.37e-9, Ca4CaN_dyad => 7.35e-5, CaMCa4CaN_dyad => 2.43e-3, 
        Ca2CaMCa4CaN_dyad => 0.013, Ca4CaMCa4CaN_dyad => 3.6, CaM_sl => 4.42e-2, Ca2CaM_sl => 7.34e-5, Ca4CaM_sl => 8.89e-9, 
        CaMB_sl => 2.44, Ca2CaMB_sl => 11.86, Ca4CaMB_sl => 4.38e-4, Pb2_sl => 1.47e-5, Pb_sl => 6.31e-6, Pt_sl => 6.60e-8, 
        Pt2_sl => 7.37e-13, Pa_sl => 4.37e-9, Ca4CaN_sl => 5.22e-4, CaMCa4CaN_sl => 1.98e-6, Ca2CaMCa4CaN_sl => 5.02e-6, 
        Ca4CaMCa4CaN_sl => 1.43e-3, CaM_cyt => 4.4e-2, Ca2CaM_cyt => 4.11e-5, Ca4CaM_cyt => 6.17e-10, CaMB_cyt => 4.179, 
        Ca2CaMB_cyt => 1.11, Ca4CaMB_cyt => 1.61e-5, Pb2_cyt => 8.23e-6, Pb_cyt => 4.15e-8, Pt_cyt => 2.29e-13, 
        Pt2_cyt => 2.47e-18, Pa_cyt => 1.53e-14, Ca4CaN_cyt => 1.17e-4, CaMCa4CaN_cyt => 4.39e-7, Ca2CaMCa4CaN_cyt => 1.59e-7, 
        Ca4CaMCa4CaN_cyt => 1.49e-6, LCC_PKAp => 16.454, LCC_CKdyadp =>16.934 , RyR2809p => 297.36, RyR2815p => 76.985, 
        PLBT17p => 0.614, LCC_CKslp => 8.66e-6, LR => -6.7e-36, LRG => 2.46e-34, RG => 4.8e-4, b1AR_S464 => 5.97e-35, 
        b1AR_S301 => 6.48e-4, GsaGTPtot => 9.6e-3, GsaGDP => 6.21e-4, Gsby => 0.01, AC_GsaGTP => 1.42e-3, PDEp => 2.22e-3, 
        cAMPtot => 1.023, RC_I => 0.804, RCcAMP_I => 0.142, RCcAMPcAMP_I => 4.48e-3, RcAMPcAMP_I => 0.229, PKACI => 8.55e-2, 
        PKACI_PKI => 0.144, RC_II => 0.051, RCcAMP_II => 8.99e-3, RCcAMPcAMP_II => 2.84e-4, RcAMPcAMP_II => 5.77e-2, 
        PKACII => 2.15e-2, PKACII_PKI => 3.62e-2, I1p_PP1 => 7.27e-2, I1ptot => 7.28e-2, PLBp => 8.454, PLMp => 5.6, 
        LCCap => 5.49e-3, LCCbp => 6.27e-3, RyRp => 2.76e-2, TnIp => 4.389, KS79 => 1.53e-3, KS80 => 1.53e-3, KSp => 1.84e-3, 
        CFTRp => 4.06e-3, KURp => 1.09e-2,
        CaM0 => 2.82e-5, Ca2CaM_C => 1.01e-8, Ca2CaM_N => 1.40e-9, Ca4CaM => 4.78e-13,
        CaM0_CaMK => 1.29e-6, Ca2CaM_C_CaMK => 9.13e-8, Ca2CaM_N_CaMK => 3.74e-9, Ca4CaM_CaMK => 5.92e-10,
        CaM0_CaMKP => 2.36e-7, Ca2CaM_C_CaMKP => 1.13e-7, Ca2CaM_N_CaMKP => 1.54e-9, Ca4CaM_CaMKP => 7.82e-10,
        CaMK => 6.73e-5, CaMKP => 6.57e-7, CaMKP2 => 2.66e-7], 
        (0.0, 10000.0), 
        [k_1C_on => 5e3, k_1C_off => 50e-3, k_2C_on => 10e3, k_2C_off => 10e-3, 
        k_1N_on => 100e3, k_1N_off => 2000e-3, k_2N_on => 200e3, k_2N_off => 500e-3,
        k_K1C_on => 44e3, k_K1C_off => 33e-3, k_K2C_on => 44e3, k_K2C_off => 0.8e-3,
        k_K1N_on => 76e3, k_K1N_off => 300e-3, k_K2N_on => 76e3, k_K2N_off => 20e-3,
        kCaM0_on => 3.8, kCaM2C_on => 0.92e3, kCaM2N_on => 0.12e3, kCaM4_on => 30e3, 
        kCaM0_off => 5.5e-3, kCaM2C_off => 6.8e-3, kCaM2N_off => 1.7e-3, kCaM4_off => 1.5e-3, 
        kCaM0P_on => 3.8*binding_To_PCaMK, kCaM2CP_on => 0.92e3*binding_To_PCaMK, 
        kCaM2NP_on => 0.12e3*binding_To_PCaMK, kCaM4P_on => 30e3*binding_To_PCaMK,
        kCaM0P_off => 1e-3/decay_CaM, kCaM2CP_off => 1e-3/decay_CaM, kCaM2NP_off => 1e-3/decay_CaM, kCaM4P_off => 1e-3/decay_CaM,
        k_phosCaM => 30e-3 * phospho_rate, k_dephospho => 1e-3/6 * phosphatase, k_P1_P2 => 1e-3/60, k_P2_P1 => 1e-3/6*0.25, CaMKII_T => 70e-6])

sol = solve(oprob, CVODE_BDF(), abstol = 1e-10, reltol = 1e-10, tstops = 0:1000:10000)

plot(sol, idxs=Cai, linewidth=1.5, title="Calcium Transient", xlabel="Time(s)", ylabel="[Ca2+](M)", label="ISO=0.1", xlim=(0,1050))

plot(sol, idxs=Cai, linewidth=1.5, title="Calcium Transient (Control)", xlabel="Time(s)", ylabel="[Ca2+](M)", label="Rodas5, tol=1e-12", xlim=(0,10.15e3))

plot(sol, idxs=Vm, linewidth=1.5, title="Action Potential (Control)", xlabel="Time(ms)", ylabel="Voltage (mV)",ylim=(-90,60),xlim=(0,10.15e3),label="Rodas5, tol=1e-12")
