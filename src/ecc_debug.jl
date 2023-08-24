#=
    This module describes the EC coupling starting from the framework of the
    Shannon-Bers model, and this file was built upon the code developeded
    by Yang and Saucerman.
    Reference: Yang JH & Saucerman JJ. (2012). Phospholemman is a negative
    feed-forward regulator of Ca2+ in beta-adrenergic signaling,
    accelerating beta-adrenergic inotropy. Journal of Molecular and Cellular
    Cardiology 52, 1048-1055.
=#
using DifferentialEquations
using ModelingToolkit

export get_ecc_equations


# Flags
const freq = 1                    # [Hz] - CHANGE DEPENDING ON FREQUENCY
const cycleLength = 1e3/freq      # [ms]
const ICa_MarkovFlag = 1          # Set ICa_MarkovFlag to 1 for Markov ICa, 0 otherwise
const INa_MarkovFlag = 0          # Set INa_MarkovFlag to 1 for Markov INa, 0 otherwise
const NaClampFlag = 0             # Na clamp (set flag to 1, 0 otherwise)
const PLMkoFlag = 0               # PLM KO (set flag to 1, 0 otherwise)
const StrophFlag = 0              # Strophanthidin (dyastolic Na influx) (set flag to 1, 0 otherwise)
const CaffeineFlag = 0            # Caffeine-induced Ca transient (set flag to 1, 0 otherwise)
const DigitalisFlag = 0           # Digitalis (set flag to 1, 0 otherwise)
const CKIIOE = 0            # "WT", "KO" = 0  "OE" = 1
const CKIIflag = CKIIOE     # Set CKIIflag to 1 for CaMKII-OE, 0 otherwise
const ItoFlag = 1           # Set Ito to use either origl params (=0) or Grandi Params (=1)
# CaMKII-Na-Ca-CaMKII loop closed (set flag to 1, 0 otherwise)
const loop = 0 # (default 0)

## Na loading and CaMKII-Na-Ca-CaMKII loop properties
if CKIIflag == 1                # Na loading parameters ON (set flag to 1, 0 otherwise)
    const NaGainFlag = 1  # CaMKII-OE (default 1)
    const inashift = -3.25
else
    const NaGainFlag = 0   # WT (default 0)
    const inashift = 0
end

## Simulation type

const protocol = 1        # "pace"

if CaffeineFlag==1
    const protocol = 2    # "vcRest"
end
if StrophFlag == 1
    const protocol = 3    # "none"
end


function get_ecc_equations()
    @variables t
    D = Differential(t)

    # Variables 1-13
    @variables Na_m(t) Na_h(t) Na_j(t) ICa_HH4(t) ICa_HH5(t) ICa_HH6(t) ICa_HH7(t) Itos_x(t) Itos_y(t) Itof_x(t) Itof_y(t) Ikr(t) IKs(t)
    # Variables 14-25 (22 not used)
    @variables RyR_R(t) RyR_O(t) RyR_I(t) NaBj(t) NaBsl(t) TnCL(t) TnCHc(t) TnCHm(t) CaM(t) Myosin_ca(t) Myosin_mg(t) SRB(t)
    # Variables 26-38  
    @variables SLLj(t) SLLsl(t) SLHj(t) SLHsl(t) Csqn(t) Ca_sr(t) Naj(t) Nasl(t) Nai(t) Ki(t) Ca_j(t) Ca_sl(t) Cai(t)       
    # Variables 39, 40, 43-47 (41-42 not used) / 84-87
    @variables Vm(t) Itos_r(t) influx_LTCC(t) influx_PMCA(t) influx_NCX(t) influx_ICa(t) Na_late_h(t) IKs_x(t) IKs1_y(t) Iss(t) IKs2_y(t)
    # Variables 48-59
    @variables CNa2(t) CNa1(t) ONa(t) IFNa(t) I1Na(t) CNa3(t) ICNa2(t) ICNa3(t) LONa(t) LCNa1(t) LCNa2(t) LCNa3(t) 
    # Variables 60-71
    @variables C2_m1j(t) C1_m1j(t) I1Ca_m1j(t) I2Ca_m1j(t) I1Ba_m1j(t) I2Ba_m1j(t) C2_m2j(t) C1_m2j(t) I1Ca_m2j(t) I2Ca_m2j(t) I1Ba_m2j(t) I2Ba_m2j(t)
    # Variables 72-83
    @variables C2_m1sl(t) C1_m1sl(t) I1Ca_m1sl(t) I2Ca_m1sl(t) I1Ba_m1sl(t) I2Ba_m1sl(t) C2_m2sl(t) C1_m2sl(t) I1Ca_m2sl(t) I2Ca_m2sl(t) I1Ba_m2sl(t) I2Ba_m2sl(t)

    ## Adjusting Variables
    @variables LCC_CKdyadp(t) RyR2815p(t) PLBT17p(t)  # fractional CaMKII-dependent LCC dyad / RyR / PLB phosphorylation
    @variables PLMp(t) PLBp(t) LCCap(t) LCCbp(t) RyRp(t) TnIp(t) KURp(t) # BAR (PKA phosphorylation) module

    @parameters LCCtotDyad = 31.4*.9        # [uM] - Total Dyadic [LCC] - (umol/l dyad)
    @parameters RyRtot = 382.6              # [uM] - Total RyR (in Dyad)
    @parameters plb_val=106                 # MOUSE
    @parameters LCCtotBA = 0.025            # [uM] - [umol/L cytosol]
    @parameters PLMtotBA = 48               # [uM] - [umol/L cytosol] MOUSE
    @parameters RyRtotBA = 0.135            # [uM] - [umol/L cytosol]
    @parameters TnItotBA = 70               # [uM] - [umol/L cytosol]
    @parameters IKurtotBA = 0.025           # [uM] - [umol/L cytosol] MOUSE
    PLBtot = plb_val                        # [uM] - Total [PLB] in cytosolic units
    PLBtotBA = plb_val                      # [uM] - [umol/L cytosol]

    LCC_CKp = LCC_CKdyadp / LCCtotDyad 
    RyR_CKp = RyR2815p / RyRtot 
    PLB_CKp = PLBT17p / PLBtot 
    LCCa_PKAp = LCCap / LCCtotBA 
    LCCb_PKAp = LCCbp / LCCtotBA 
    PLB_PKAn = (PLBtotBA - PLBp) / PLBtotBA 
    RyR_PKAp = RyRp / RyRtotBA 
    TnI_PKAp = TnIp / TnItotBA 
    IKur_PKAp = KURp / IKurtotBA 
    PLM_PKAp = PLMp / PLMtotBA


    # Parameters
    # Constants
    @parameters R = 8314        # [J/kmol*K]  
    @parameters Frdy = 96485    # [C/mol]  
    @parameters Temp = 310      # [K] 310 K (37 C) for BT / 295 K (22 C) for RT
    FoRT = Frdy/R/Temp
    Qpow = (Temp-310)/10

    # Cell geometry
    @parameters Acell = 20e3        # [um^2] MOUSE
    Cmem = Acell*1e-14  # [F] 200 pF membrane capacitance MOUSE

    # Fractional currents in compartments
    @parameters Fjunc = 17/(17+31)*7/17+31/(17+31)*2/31
    Fsl = 1-Fjunc
    @parameters Fjunc_nak = 1.6*17/(1.6*17+31)*7/17+31/(1.6*17+31)*2/31 
    Fsl_nak = 1-Fjunc_nak
    Fjunc_ncx = Fjunc
    Fsl_ncx = 1-Fjunc_ncx
    @parameters Fjunc_CaL = 0.9 
    Fsl_CaL = 1-Fjunc_CaL

    @parameters cellLength = 100                # cell length [um]
    @parameters cellRadius = 10.25              # cell radius [um]
    @parameters junctionLength = 15e-3          # junc length [um]
    @parameters junctionRadius = 160e-3         # junc radius [um]
    @parameters distSLcyto = 0.45               # dist. SL to cytosol [um]
    @parameters distJuncSL = 0.3                # dist. junc to SL [um] MOUSE
    @parameters DcaJuncSL = 1.64e-6             # Dca junc to SL [cm^2/sec]
    @parameters DcaSLcyto = 1.22e-6             # Dca SL to cyto [cm^2/sec]
    @parameters DnaJuncSL = 1.09e-5             # Dna junc to SL [cm^2/sec]
    @parameters DnaSLcyto = 1.79e-5             # Dna SL to cyto [cm^2/sec] 
    Vcell = pi*cellRadius^2*cellLength*1e-15 # [L]
    Vmyo = 0.65*Vcell 
    Vsr = 0.035*Vcell 
    Vsl = 0.02*Vcell 
    Vjunc = 0.0539*.01*Vcell
    SAsl = Fsl*Acell                                    # [um^2]  MOUSE
    Njunc = (Fjunc*Acell)/(pi*junctionRadius^2)         # [-]
    SAjunc = Njunc*pi*2*junctionLength*junctionRadius   # [um^2] MOUSE

    J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10     # [L/msec] [m^2/sec] MOUSE
    J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10        # MOUSE
    J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10     # MOUSE
    J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10        # MOUSE

    # Fixed ion concentrations     
    @parameters Cli = 15    # Intracellular Cl  [mM]
    @parameters Clo = 150   # Extracellular Cl  [mM]
    @parameters Ko = 5.4    # Extracellular K   [mM]
    @parameters Nao = 140   # Extracellular Na  [mM]
    @parameters Cao = 1     # Extracellular Ca  [mM] MOUSE # 1.8 mM in RABBIT
    @parameters Mgi = 1     # Intracellular Mg  [mM]

    # Nernst Potentials
    ena_junc = (1/FoRT)*log(Nao/Naj)        # [mV]
    ena_sl = (1/FoRT)*log(Nao/Nasl)         # [mV]
    ek = (1/FoRT)*log(Ko/Ki)                # [mV]
    eca_junc = (1/FoRT/2)*log(Cao/Ca_j)     # [mV]
    eca_sl = (1/FoRT/2)*log(Cao/Ca_sl)      # [mV]
    ecl = (1/FoRT)*log(Cli/Clo)             # [mV]

    ## Na transport parameters
    @parameters GNa = 10                # [mS/uF] changed from rabbit (16)
    @parameters GNaB = 4.5*0.297e-3     # [mS/uF] changed from rabbit
    @parameters IbarNaK = 5             # [uA/uF] changed from rabbit (1.90719)
    if NaGainFlag == 1
        GNaB=GNaB*4
        IbarNaK = IbarNaK*0.9
    end
    if DigitalisFlag == 1
        IbarNaK = IbarNaK*0.5 # 50# block
    end
    @parameters KmNaip = 19     # [mM] changed from rabbit (11)
    @parameters KmKo = 1.5      # [mM]
    @parameters Q10NaK = 1.63
    @parameters  Q10KmNai = 1.39
    if PLMkoFlag == 1
        PLM_PKAp=1
        GNaB=GNaB*48/20
        IbarNaK = IbarNaK * 0.8
    end
    if StrophFlag == 1
        IbarNaK = 0
    end

    # INa Markov Model parameters
    @parameters GNa2 = 10.64        # [mS/uF]
    @parameters P1a1 = 3.802
    @parameters P2a1 = 0.1027
    @parameters P3a1 = 2.5
    @parameters P4a1 = 17
    @parameters P5a1 = 0.20
    @parameters P6a1 = 150
    @parameters P4a2 = 15
    @parameters P5a2 = 0.23
    @parameters P4a3 = 12
    @parameters P5a3 = 0.25
    @parameters P1b1 = 0.1917
    @parameters P2b1 = 20.3
    @parameters P1b2 = 0.2
    @parameters P2b2 = 2.5
    @parameters P1b3 = 0.22
    @parameters P2b3 = 7.5
    @parameters P1a4 = 0.188495
    @parameters P2a4 = 16.6
    @parameters P3a4 = 0.393956
    @parameters P4a4 = 7
    @parameters P1a5 = 7e-7
    @parameters P2a5 = 7.2 # TG 7.8
    @parameters P1b5 = 0.0044 # TG 0.00467
    @parameters P2b5 = 2e-5
    @parameters P1a6 = 100
    @parameters P1b6 = 8.9554e-7 # TG 6.5449e-7
    @parameters P2b6 = 11.3944
    @parameters P1a7 = 0.487e-4 # TG 0.3377e-3
    @parameters P2a7 = 23.2696
    @parameters P1b7 = 0.2868e-3 # TG 1.868e-4
    @parameters P2b7 = 35.9898
    @parameters P1a8 = 0.1e-7 # TG 6.5e-6
    @parameters P1b8 = 9.8e-3 # TG 3.8e-3
    if CKIIflag == 1 # MOUSE - CaMKII-OE
        @parameters P2a5 = 7.8
        @parameters P1b5 = 0.00467
        @parameters P1b6 = 6.5449e-7
        @parameters P1a7 = 0.3377e-3
        @parameters P1b7 = 1.868e-4
        @parameters P1a8 = 6.5e-6
        @parameters P1b8 = 3.8e-3
    end

    ## K currents parameters
    @parameters Kcoeff = 1 # K current modulation
    @parameters pNaK = 0.01833
    @parameters GtoSlow = 0 # [mS/uF] changed from rabbit (0.06): NO ItoSlow in MOUSE
    @parameters GtoFast = 0.44 # [mS/uF] changed from rabbit (0.02)
    if CKIIflag == 1 # MOUSE
        GtoFast = GtoFast * 2/3 # chronic CaMKII-OE effect
    end

    #Gkur = 0.3  # [mS/uF] only in MOUSE
    @parameters Gkur1 = 1.1 * 0.16 # fast
    @parameters Gkur2 = 0.14 # slow
    @parameters Gss = 0.15 # [mS/uF] only in MOUSE
    gkr = 0.03 * sqrt(Ko/5.4)
    @parameters gkp = 0.001

    # Cl current parameters
    @parameters GClCa = 0.109625 # [mS/uF]
    @parameters GClB = 9e-3 # [mS/uF]
    @parameters KdClCa = 100e-3 # [mM]

    ## LTCC parameters
    @parameters K_Ica = 1.65 # MOUSE
    pNa = K_Ica * 1.5e-8 # [cm/sec]
    pCa = K_Ica * 5.4e-4 # [cm/sec] - Ca permeability
    pK = K_Ica * 2.7e-7 # [cm/sec]
    @parameters KmCa = 0.6e-3 # [mM]
    @parameters Q10CaL = 1.8

    ## Ca transport parameters
    @parameters IbarNCX = 1 # [uA/uF] changed from rabbit (9)
    if CKIIflag == 1
        IbarNCX = 1.5 * IbarNCX
    end
    @parameters KmCai = 3.59e-3         # [mM]
    @parameters KmCao = 1.3             # [mM]
    @parameters KmNai = 12.29           # [mM]
    @parameters KmNao = 87.5            # [mM]
    @parameters ksat = 0.27             # [none]  
    @parameters nu = 0.35               # [none]
    @parameters Kdact = 1/2 * 0.256e-3  # [mM] changed from rabbit
    @parameters Q10NCX = 1.57           # [none]
    @parameters IbarSLCaP = 0.0673      # [uA/uF]
    @parameters KmPCa = 0.5e-3          # [mM]
    @parameters GCaB = 3 * 2.513e-4     # [uA/uF] changed from rabbit (2.513e-4)
    @parameters Q10SLCaP = 2.35         # [none]

    # SR flux parameters
    @parameters Q10SRCaP = 2.6                  # [none]
    @parameters Vmax_SRCaP = 1.15*1.15*2.86e-4  # [mM/msec] (mmol/L cytosol/msec) changed
    @parameters Kmf = 0.3e-3                    # [mM] changed from rabbit (0.246e-3) # from Yang-Saucerman
    @parameters Kmr = 2.1                       # [mM]L cytosol changed from rabbit (1.7) # from Yang-Saucerman
    @parameters hillSRCaP = 1.787               # [mM]
    @parameters ks = 25                         # [1/ms]      
    @parameters koCa = 10                       # [mM^-2 1/ms]
    @parameters kom = 0.06                      # [1/ms]     
    @parameters kiCa = 0.5                      # [1/mM/ms]
    @parameters kim = 0.005                     # [1/ms]
    @parameters ec50SR = 0.45 + 0.05            # [mM] changed from rabbit (0.45)

    if CaffeineFlag == 1
        koCa = koCa * 7.5
        @parameters GCaB = 0
        @parameters Vmax_SRCaP = 0
    end

    ## Buffering parameters
    @parameters Bmax_Naj = 7.561        # [mM] 
    @parameters Bmax_Nasl = 1.65        # [mM]
    @parameters koff_na = 1e-3          # [1/ms]
    @parameters kon_na = 0.1e-3         # [1/mM/ms]
    @parameters Bmax_TnClow = 70e-3     # [mM]                      # TnC low affinity
    @parameters koff_tncl = 19.6e-3     # [1/ms] 
    @parameters kon_tncl = 32.7         # [1/mM/ms]
    @parameters Bmax_TnChigh = 140e-3   # [mM]                      # TnC high affinity 
    @parameters koff_tnchca = 0.032e-3  # [1/ms] 
    @parameters kon_tnchca = 2.37       # [1/mM/ms]
    @parameters koff_tnchmg = 3.33e-3   # [1/ms] 
    @parameters kon_tnchmg = 3e-3       # [1/mM/ms]
    @parameters Bmax_myosin = 140e-3    # [mM]                      # Myosin buffering
    @parameters koff_myoca = 0.46e-3    # [1/ms]
    @parameters kon_myoca = 13.8        # [1/mM/ms]
    @parameters koff_myomg = 0.057e-3   # [1/ms]
    @parameters kon_myomg = 0.0157      # [1/mM/ms]
    @parameters Bmax_SR = 19 * 0.9e-3   # [mM] 
    @parameters koff_sr = 60e-3         # [1/ms]
    @parameters kon_sr = 100            # [1/mM/ms]
    Bmax_SLlowsl = 37.38e-3 * Vmyo/Vsl          # [mM]    # SL buffering
    Bmax_SLlowj = 4.62e-3 * Vmyo/Vjunc * 0.1    # [mM]    
    @parameters koff_sll = 1300e-3      # [1/ms]
    @parameters kon_sll = 100           # [1/mM/ms]
    Bmax_SLhighsl = 13.35e-3 * Vmyo/Vsl         # [mM] 
    Bmax_SLhighj = 1.65e-3 * Vmyo/Vjunc * 0.1   # [mM] 
    @parameters koff_slh = 30e-3        # [1/ms]
    @parameters kon_slh = 100           # [1/mM/ms]
    @parameters Bmax_Csqn = 2.7         # 140e-3*Vmyo/Vsr  [mM] 
    @parameters koff_csqn = 65          # [1/ms] 
    @parameters kon_csqn = 100          # [1/mM/ms] 

    # PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
    @parameters fracTnIpo = 0.062698  # Derived quantity (TnI_PKAp(baseline)/TnItot)
    fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
    koff_tncl = koff_tncl * fPKA_TnI

    ## I_Na: Fast Na Current
    # Max INa alterations with CaMKII hyperactivity as in Hund & Rudy 2008
    if CKIIflag == 1 # acute effects
        @parameters alphaCKII = -0.18
        if NaGainFlag == 1
            deltGbarNal_CKII = 3  # MOUSE
        else
            deltGbarNal_CKII = 0  # no Na Gain in OE
        end
    else
        @parameters alphaCKII = 0
        deltGbarNal_CKII = 0
    end
    
    if loop == 1
        @parameters RyRp_WT_mean=0.2101 
        @parameters RyRp_OE_mean=0.7387 # Derived (1 Hz, no loop)
        @parameters RyRp_OEloop_min=0.7033 # Derived (1 Hz, OE loop)
        delta_loop=(3/(RyRp_OE_mean-RyRp_WT_mean))*RyR_CKp-(3/(RyRp_OE_mean-RyRp_WT_mean))*RyRp_WT_mean
        @parameters NaVsCaMKIIclamp=0 # if 1, CaMKII Clamp on NaV
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
    
    am = 0.32 * (Vm + 47.13)/(1-exp(-0.1*(Vm + 47.13)))
    bm = 0.08 * exp(-(Vm)/11)
    ah = ifelse((Vm-inashift) >= -40, 0, 0.135*exp((80+(Vm-inashift))/-6.8)) 
    aj = ifelse((Vm-inashift) >= -40, 0, (1+alphaCKII)*((-1.2714e5*exp(0.2444*(Vm-inashift))-3.474e-5*exp(-0.04391*(Vm-inashift)))*((Vm-inashift)+37.78)/(1+exp(0.311*((Vm-inashift)+79.23)))))
    bh = ifelse((Vm-inashift) >= -40, 0.66*1/(0.13*(1+exp(-((Vm-inashift)+10.66)/11.1))), 1.1*3.56*exp(0.079*(Vm-inashift-2))+3.1e5*exp(0.35*(Vm-inashift-2))) # MOUSE
    bj = ifelse((Vm-inashift) >= -40, 0.3*exp(-2.535e-7*(Vm-inashift))/(1+exp(-0.1*((Vm-inashift)+32))), 0.1212*exp(-0.01052*(Vm-inashift))/(1+exp(-0.1378*((Vm-inashift)+40.14))))
    I_Na_junc1 = Fjunc*GNa*(Na_m^3)*Na_h*Na_j*(Vm-ena_junc)
    I_Na_sl1 = Fsl*GNa*(Na_m^3)*Na_h*Na_j*(Vm-ena_sl)
    
    
    INa_fast_eqs = [
        D(Na_m) ~ am*(1-Na_m)-bm*Na_m,          # du[1]
        D(Na_h) ~ ah*(1-Na_h)-bh*Na_h,          # du[2]
        D(Na_j) ~ aj*(1-Na_j)-bj*Na_j           # du[3]
    ]
        

    ## I_Na,L: Late INa current (as in Hund & Rudy 2008)

    @parameters tauhl = 600 # ms
    # h-gate (note: m-gate is same as INa m-gate -> using y(1) for this)
    GbarNal = 0.0065*(1+deltGbarNal_CKII)*2 # deltGbar assigned in 'Fast INa' section
    hlss = 1/(1+exp((Vm+91)/6.1))
    I_Nalj = Fjunc*GbarNal*(Na_m^3)*Na_late_h*(Vm-ena_junc)
    I_Nalsl = Fsl*GbarNal*(Na_m^3)*Na_late_h*(Vm-ena_sl)


    ## I_Na: alternative Markov Model (unused) + compute total current (fast and late components)
    
    alphaNa8 = P1a8
    betaNa8 = P1b8    
    I2Na = (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3))
    # Transition rates
    alphaNa1 = P1a1/(P2a1*exp(-(Vm+P3a1)/P4a1)+P5a1*exp(-(Vm+P3a1)/P6a1))
    alphaNa2 = P1a1/(P2a1*exp(-(Vm+P3a1)/P4a2)+P5a2*exp(-(Vm+P3a1)/P6a1))
    alphaNa3 = P1a1/(P2a1*exp(-(Vm+P3a1)/P4a3)+P5a3*exp(-(Vm+P3a1)/P6a1))
    betaNa1 = P1b1*exp(-(Vm+P3a1)/P2b1) # shift
    betaNa2 = P1b2*exp(-(Vm-P2b2)/P2b1)
    betaNa3 = P1b3*exp(-(Vm-P2b3)/P2b1)
    alphaNa4 = 1/(P1a4*exp(-(Vm+P4a4)/P2a4)+P3a4)
    alphaNa5 = P1a5*exp(-(Vm+P4a4)/P2a5)
    betaNa5 = (P1b5+P2b5*(Vm+P4a4))
    betaNa6 = P1b6*exp(-Vm/P2b6)
    alphaNa7 = P1a7*exp(Vm/P2a7)
    betaNa7 = P1b7*exp(-Vm/P2b7)
    betaNa4 = (alphaNa3*alphaNa4*alphaNa5)/(betaNa3*betaNa5)
    alphaNa6 = alphaNa4/P1a6
    I_Na_junc2 = Fjunc*GNa2*(ONa+LONa)*(Vm-ena_junc)
    I_Na_sl2 = Fsl*GNa2*(ONa+LONa)*(Vm-ena_sl)
    # I_Na: compute total current (fast and late components)
    I_Na_junc = ifelse(INa_MarkovFlag == 1, I_Na_junc2, I_Na_junc1 + I_Nalj)
    I_Na_sl = ifelse(INa_MarkovFlag == 1, I_Na_sl2, I_Na_sl1 + I_Nalsl)


    INa_eqs = [
        D(Na_late_h) ~ (hlss-Na_late_h)/tauhl,
        D(CNa3) ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LCNa3+betaNa1*CNa2+alphaNa5*ICNa3-(alphaNa1+betaNa5+alphaNa8)*CNa3), 0),
        D(CNa2) ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LCNa2+alphaNa1*CNa3+betaNa2*CNa1+alphaNa5*ICNa2-(betaNa1+alphaNa2+betaNa5+alphaNa8)*CNa2), 0),
        D(CNa1) ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LCNa1+alphaNa2*CNa2+betaNa3*ONa+alphaNa5*IFNa-(betaNa2+alphaNa3+betaNa5+alphaNa8)*CNa1), 0),
        D(ONa)  ~ ifelse(INa_MarkovFlag == 1, (betaNa8*LONa+alphaNa3*CNa1+betaNa4*IFNa-(betaNa3+alphaNa4+alphaNa8)*ONa), 0),
        D(IFNa) ~ ifelse(INa_MarkovFlag == 1, (alphaNa4*ONa+betaNa5*CNa1+betaNa6*I1Na+alphaNa2*ICNa2-(betaNa4+alphaNa5+alphaNa6+betaNa2)*IFNa), 0),
        D(I1Na) ~ ifelse(INa_MarkovFlag == 1, (alphaNa6*IFNa+betaNa7*I2Na-(betaNa6+alphaNa7)*I1Na), 0),
        D(ICNa2) ~ ifelse(INa_MarkovFlag == 1, (alphaNa1*ICNa3+betaNa2*IFNa+betaNa5*CNa2-(betaNa1+alphaNa2+alphaNa5)*ICNa2), 0),
        D(ICNa3) ~ ifelse(INa_MarkovFlag == 1, (betaNa1*ICNa2+betaNa5*CNa3-(alphaNa1+alphaNa5)*ICNa3), 0),
        D(LONa)  ~ ifelse(INa_MarkovFlag == 1, (alphaNa3*LCNa1+alphaNa8*ONa-(betaNa8+betaNa3)*LONa), 0),
        D(LCNa1) ~ ifelse(INa_MarkovFlag == 1, (alphaNa8*CNa1+alphaNa2*LCNa2+betaNa3*LONa-(betaNa8+betaNa2+alphaNa3)*LCNa1), 0),
        D(LCNa2) ~ ifelse(INa_MarkovFlag == 1, (betaNa2*LCNa1+alphaNa8*CNa2+alphaNa1*LCNa3-(betaNa8+betaNa1+alphaNa2)*LCNa2), 0),
        D(LCNa3) ~ ifelse(INa_MarkovFlag == 1, (alphaNa8*CNa3+betaNa1*LCNa2-(betaNa8+alphaNa1)*LCNa3), 0)
    ]

    ## I_nabk: Na Background Current
    I_nabk_junc = Fjunc*GNaB*(Vm-ena_junc)
    I_nabk_sl = Fsl*GNaB*(Vm-ena_sl)


    ## I_nak: Na/K Pump Current
    @parameters sigma = (exp(Nao/67.3)-1)/7
    @parameters fracPKA_PLMo = 0.116738     # Derived quantity (PLM_PKAp(baseline)/PLMtot)
    @parameters fracPKA_PLMiso = 0.859251   # Derived quantity (PLM_PKAp(ISO)/PLMtot)
    kPKA_PLM=KmNaip*(1-0.7019)/(fracPKA_PLMiso/fracPKA_PLMo-1)  # PLM_PKAp ISO
    KmNaip_PKA = -kPKA_PLM+kPKA_PLM*(PLM_PKAp/fracPKA_PLMo)
    KmNaip = KmNaip-KmNaip_PKA
    fnak = 1/(1+0.1245*exp(-0.1*Vm*FoRT)+0.0365*sigma*exp(-Vm*FoRT))
    I_nak_junc = Fjunc_nak*IbarNaK*fnak*Ko /(1+(KmNaip/Naj)^4) /(Ko+KmKo)
    I_nak_sl = Fsl_nak*IbarNaK*fnak*Ko /(1+(KmNaip/Nasl)^4) /(Ko+KmKo)
    I_nak = I_nak_junc+I_nak_sl
    

    ## IK
    #I_kur - IK,slow
    xurss = 1/(1+exp(-(Vm+15)/14)) 
    yurss = 1/(1+exp((Vm+48)/6.2)) 
    tauxur = 0.95+0.05*exp(-0.08*Vm)
    tauxur2 = 1+7/(1+exp(-(Vm+45)/8))+20*exp(-((Vm+35)/10)^2) #8+20*exp(-((Vm+35)/10)^2) 
    tauyur1 = 400+900*exp(-((Vm+55)/16)^2)-250/(1+exp(-(Vm+60)/8))  # fast
    tauyur2 = 400+900*exp(-((Vm+55)/16)^2)+550/(1+exp(-(Vm+60)/8))  # slow
    
    # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    @parameters fracIKurp0 = 0.437635       # Derived quantity (IKur_PKAp(baseline)/IKurtot)
    @parameters fracIKurpISO = 0.718207     # Derived quantity (IKur_PKAp(ISO)/IKurtot)
    a_Kur = (1.20-1)/(fracIKurpISO/fracIKurp0-1) 
    fracIKuravail = (1-a_Kur)+a_Kur*(IKur_PKAp/fracIKurp0)  # +20# with 0.1 uM ISO
    I_kur1 = Kcoeff*fracIKuravail*Gkur1*IKs_x*IKs1_y*(Vm-ek)  # IKslow1
    I_kur2 = Kcoeff*Gkur2*IKs_x*IKs2_y*(Vm-ek)  # IKslow2 # no PKA effect
    I_kur = I_kur1 + I_kur2 
    
    # I_ss
    xssss = xurss
    tauxss = 70*exp(-((Vm+43)/30)^2)+14  # Iss
    I_ss = Kcoeff*Gss*Iss*(Vm-ek)  #store
    
    # I_kr: Rapidly Activating K Current
    xrss = 1/(1+exp(-(Vm+50)/7.5)) 
    tauxr = 1/(1.38e-3*(Vm+7)/(1-exp(-0.123*(Vm+7)))+6.1e-4*(Vm+10)/(exp(0.145*(Vm+10))-1)) 
    rkr = 1/(1+exp((Vm+33)/22.4)) 
    I_kr = Kcoeff*gkr*Ikr*rkr*(Vm-ek) 
    
    IK_eqs = [
        D(IKs_x) ~ (xurss-IKs_x)/tauxur,    # IKslow1
        D(IKs) ~ (xurss-IKs)/tauxur2,       # IKslow2
        D(IKs1_y) ~ (yurss-IKs1_y)/tauyur1, 
        D(IKs2_y) ~ (yurss-IKs2_y)/tauyur2,
        D(Iss) ~ (xssss-Iss)/tauxss,
        D(Ikr) ~ (xrss-Ikr)/tauxr 
    ] 

    # I_ks: Slowly Activating K Current
    @parameters I_ks_junc = 0   # No IKs in mouse
    @parameters I_ks_sl = 0     # No IKs in mouse
    I_ks = I_ks_junc+I_ks_sl 

    # I_kp: Plateau K current
    kp_kp = 1/(1+exp(7.488-Vm/5.98)) 
    I_kp_junc = Kcoeff*Fjunc*gkp*kp_kp*(Vm-ek) 
    I_kp_sl = Kcoeff*Fsl*gkp*kp_kp*(Vm-ek) 
    I_kp = I_kp_junc+I_kp_sl 


    ## I_to: Transient Outward K Current (slow and fast components)
    xtoss = 1/(1+exp(-(Vm+3.0)/13))
    ytoss = 1/(1+exp((Vm+48)/5))
    rtoss = 1/(1+exp((Vm+33.5)/10)) # Rto not used in MOUSE model
    tauxtos = 0.08+0.7*exp(-((Vm+25)/30)^2)
    Pr1 = ifelse(CKIIflag == 0, 8085, 3600)
    Pr2 = ifelse(CKIIflag == 0, 313, 500)
    taurtos = ifelse(ItoFlag == 0, 2.8e3/(1+exp((Vm+60.0)/10))+220, Pr1/(1+exp((Vm+33.5)/10))+Pr2) # no Rto
    tauytos = ifelse(CKIIflag == 0, 100+400/(1+exp((Vm+25)/5)), 100+35/(1+exp((Vm+25)/5)))
    I_tos = 0*Kcoeff*GtoSlow*Itos_x*Itos_y*(Vm-ek)    # N0 Itos in MOUSE # [uA/uF]
    xtofs = 1/(1+exp(-(Vm+3.0)/13))           # = xtoss
    ytofs = 1/(1+exp((Vm+48)/5))              # = ytoss
    tauxtof = 0.08+0.7*exp(-((Vm+25)/30)^2)   # = tauxtos 
    tauytof = ifelse(CKIIflag == 1, 5+32*exp(-((Vm+55)/16)^2)+12/(1+exp(-(Vm+60)/8)), 10+32*exp(-((Vm+55)/16)^2)+8/(1+exp(-(Vm+60)/8)))
    I_tof = Kcoeff*GtoFast*Itof_x*Itof_y*(Vm-ek)
    I_to = I_tos + I_tof

    # Itos (ABSENT IN MOUSE)
    # Shannon Versions
    Ito_eqs = [
        D(Itos_x) ~ (xtoss-Itos_x)/tauxtos,     # du[8]
        D(Itos_y) ~ (ytoss-Itos_y)/tauytos,     # du[9]
        D(Itos_r) ~ (rtoss-Itos_r)/taurtos,     # Rto not used in MOUSE model, du[40]
        # Itof
        D(Itof_x) ~ (xtofs-Itof_x)/tauxtof,     # du[10]
        D(Itof_y) ~ (ytofs-Itof_y)/tauytof      # du[11]
    ]
            

    ## I_k1: Time-Independent K Current (I_k1)
    aki = 1.02/(1+exp(0.2385*(Vm-ek-59.215)))
    bki = (0.49124*exp(0.08032*(Vm+5.476-ek))+exp(0.06175*(Vm-ek-594.31)))/(1 + exp(-0.5143*(Vm-ek+4.753)))
    kiss = aki/(aki+bki)
    I_k1 = ifelse(CKIIflag == 1, 1/2*0.3*sqrt(Ko/5.4)*kiss*(Vm-ek)*Kcoeff, 0.3*sqrt(Ko/5.4)*kiss*(Vm-ek)*Kcoeff)

    
    ## I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
    I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/Ca_j)*(Vm-ecl)
    I_ClCa_sl = Fsl*GClCa/(1+KdClCa/Ca_sl)*(Vm-ecl)
    I_ClCa = I_ClCa_junc+I_ClCa_sl
    I_Clbk = GClB*(Vm-ecl)
    

    ## Original H-H formulation for LCC - unused if ICa_MarkovFlag = 1
    ibarca_j = pCa*4*(Vm*Frdy*FoRT) * (0.341*Ca_j*exp(2*Vm*FoRT)-0.341*Cao) /(exp(2*Vm*FoRT)-1)
    ibarca_sl = pCa*4*(Vm*Frdy*FoRT) * (0.341*Ca_sl*exp(2*Vm*FoRT)-0.341*Cao) /(exp(2*Vm*FoRT)-1)
    ibark = pK*(Vm*Frdy*FoRT)*(0.75*Ki*exp(Vm*FoRT)-0.75*Ko) /(exp(Vm*FoRT)-1)
    ibarna_j = pNa*(Vm*Frdy*FoRT) *(0.75*Naj*exp(Vm*FoRT)-0.75*Nao)  /(exp(Vm*FoRT)-1)
    ibarna_sl = pNa*(Vm*Frdy*FoRT) *(0.75*Nasl*exp(Vm*FoRT)-0.75*Nao)  /(exp(Vm*FoRT)-1)

    I_Ca_junc1 = (Fjunc_CaL*ibarca_j*ICa_HH4*ICa_HH5*(1-ICa_HH6)*Q10CaL^Qpow)*0.45
    I_Ca_sl1 = (Fsl_CaL*ibarca_sl*ICa_HH4*ICa_HH5*(1-ICa_HH7)*Q10CaL^Qpow)*0.45
    I_CaK1 = (ibark*ICa_HH4*ICa_HH5*(Fjunc_CaL*(1-ICa_HH6)+Fsl_CaL*(1-ICa_HH7))*Q10CaL^Qpow)*0.45
    I_CaNa_junc1 = (Fjunc_CaL*ibarna_j*ICa_HH4*ICa_HH5*(1-ICa_HH6)*Q10CaL^Qpow)*0.45
    I_CaNa_sl1 = (Fsl_CaL*ibarna_sl*ICa_HH4*ICa_HH5*(1-ICa_HH7)*Q10CaL^Qpow)*0.45

    ICa_HH_eqs = [
        D(ICa_HH4) ~ 0,     # du[4] = (dss-y(4))/taud
        D(ICa_HH5) ~ 0,     # du[5] = (fss-y(5))/(tauf)
        D(ICa_HH6) ~ 0,     # du[6] = (1.7)*y(36)*(1-y(6))-11.9e-3*y(6) # fCa_junc  
        D(ICa_HH7) ~ 0      # du[7] = 1.7*y(37)*(1-y(7))-11.9e-3*y(7) # fCa_sl
    ]
    
    
    ## LCC MARKOV MODEL - based on Mahajan et al. (2008)
    # LCC Current Fixed Parameters
    @parameters taupo = 1          # [ms] - Time constant of activation
    @parameters TBa = 450          # [ms] - Time constant
    @parameters s1o = 0.0221
    @parameters k1o = 0.03
    @parameters kop = 2.5e-3       # [mM]
    @parameters cpbar = 8e-3       # [mM]
    @parameters tca = 78.0312
    @parameters ICa_scale = 5.25
    @parameters recoveryReduc = 3
    
    # PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
    @parameters fracLCCbp0 = 0.250657 # Derived quantity - (LCCbp(baseline)/LCCbtot)
    @parameters fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
    a_favail=(1.56-1)/(fracLCCbpISO/fracLCCbp0-1) # fracLCCbp ISO (x1.56 o.1 ISO)
    favail = (1-a_favail)+a_favail*(LCCb_PKAp/fracLCCbp0)  # Test (max x2.52 100# phosph)
    ICa_scale = ICa_scale*favail 
    
    @parameters SSAshift=0 
    @parameters SSIshift=0
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
    @parameters r1 = 0.3                            # [1/ms] - Opening rate
    @parameters r2 = 3                              # [1/ms] - closing rate
    s1 = s1o*fcaj 
    @parameters s1p = 0.00195                       # [ms] - Inactivation rate
    k1 = k1o*fcaj 
    @parameters k1p = 0.00413                       # [ms] - Inactivation rate
    @parameters k2 = 1e-4                           # [ms] - Inactivation rate
    @parameters k2p = 0.00224                       # [ms] - Inactivation rate
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
    Po_LCCj_m1 = 1.0-C2_m1j-C1_m1j-I1Ca_m1j-I2Ca_m1j-I1Ba_m1j-I2Ba_m1j 

    ICaMar_m1j_eqs = [
        D(C2_m1j) ~ (betaLCC*C1_m1j + k5*I2Ca_m1j + k5p*I2Ba_m1j - (k6+k6p+alphaLCC)*C2_m1j),                     # du[60]
        D(C1_m1j) ~ (alphaLCC*C2_m1j + k2*I1Ca_m1j + k2p*I1Ba_m1j + r2*Po_LCCj_m1 - (r1+betaLCC+k1+k1p)*C1_m1j),  # du[61]
        D(I1Ca_m1j) ~ (k1*C1_m1j + k4*I2Ca_m1j + s1*Po_LCCj_m1 - (k2+k3+s2)*I1Ca_m1j),                            # du[62]
        D(I2Ca_m1j) ~ (k3*I1Ca_m1j + k6*C2_m1j - (k4+k5)*I2Ca_m1j),                                               # du[63]
        D(I1Ba_m1j) ~ (k1p*C1_m1j + k4p*I2Ba_m1j + s1p*Po_LCCj_m1 - (k2p+k3p+s2p)*I1Ba_m1j),                      # du[64]
        D(I2Ba_m1j) ~ (k3p*I1Ba_m1j + k6p*C2_m1j - (k5p+k4p)*I2Ba_m1j)                                            # du[65]
    ]

    ibarca_jm1 = (4*pCa*Vm*Frdy*FoRT)*(0.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1) 
    I_Ca_junc_m1 = (Fjunc_CaL*ibarca_jm1*Po_LCCj_m1*Q10CaL^Qpow)*ICa_scale

    
    # Re-define all parameters as mode 2 specific parameters
    @parameters s1om2 = 0.0221 
    @parameters k1om2 = 0.03 
    @parameters kopm2 = 2.5e-3 
    @parameters cpbarm2 = 8e-3 
    @parameters tcam2 = 78.0312 

    possm2 = 1/(1+exp(-(Vm+SSAshift)/8)) 
    fcajm2 = 1/(1+(kopm2/Ca_j)^3)  # Depends on junctional Ca
    Rvm2 = 10 + 4954*exp(Vm/15.6) 
    PrLCCm2 = 1-1/(1+exp(-(Vm+40)/4)) 
    PsLCCm2 = 1/(1+exp(-(Vm+40+SSIshift)/11.32)) 
    TCajm2 = (tcam2 + 0.1*(1+(Ca_j/cpbarm2)^2))/(1+(Ca_j/cpbarm2)^2)  # Caj dependent
    tauCajm2 = (Rvm2-TCajm2)*PrLCCm2 + TCajm2  # Caj dependence
    tauBam2 = (Rvm2-TBa)*PrLCCm2 + TBa 

    alphaLCCm2 = possm2/taupo 
    betaLCCm2 = (1-possm2)/taupo 
    r1m2 = 0.3                                # [1/ms] - Opening rate
    r2m2 = 3/8  # [1/ms] - closing rate,  changed from rabbit (3/10) - MOUSE
    s1m2 = s1om2*fcajm2  
    s1pm2 = .00195                            # [ms] - Inactivation rate
    k1m2 = k1om2*fcajm2  
    k1pm2 = .00413                            # [ms] - Inactivation rate
    k2m2 = 1e-4                               # [ms] - Inactivation rate
    k2pm2 = .00224                            # [ms] - Inactivation rate
    s2m2 = s1m2*(k2m2/k1m2)*(r1m2/r2m2) 
    s2pm2 = s1pm2*(k2pm2/k1pm2)*(r1m2/r2m2) 
    k3m2 = exp(-(Vm+40)/3)/(3*(1+exp(-(Vm+40)/3))) 
    k3pm2 = k3m2 
    k5m2 = (1-PsLCCm2)/tauCajm2 
    k6m2 = (fcajm2*PsLCCm2)/tauCajm2 
    k5pm2 = (1-PsLCCm2)/tauBam2 
    k5m2 = k5m2/recoveryReduc       # reduced for recovery
    k5pm2 = k5pm2/recoveryReduc     # reduced for recovery    
    k6pm2 = PsLCCm2/tauBam2 
    k4m2 = k3m2*(alphaLCCm2/betaLCCm2)*(k1m2/k2m2)*(k5m2/k6m2) 
    k4pm2 = k3pm2*(alphaLCCm2/betaLCCm2)*(k1pm2/k2pm2)*(k5pm2/k6pm2) 
    
    # State transitions for MODE 2 junctional LCCs                                      
    Po_LCCj_m2 = 1.0-C2_m2j-C1_m2j-I1Ca_m2j-I2Ca_m2j-I1Ba_m2j-I2Ba_m2j                                                              # O_m2j                                    

    ICaMar_m2j_eqs = [
        D(C2_m2j) ~ (betaLCCm2*C1_m2j + k5m2*I2Ca_m2j + k5pm2*I2Ba_m2j - (k6m2+k6pm2+alphaLCCm2)* C2_m2j),                            # du[66]
        D(C1_m2j) ~ (alphaLCCm2* C2_m2j + k2m2*I1Ca_m2j + k2pm2*I1Ba_m2j + r2m2*Po_LCCj_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*C1_m2j),     # du[67]
        D(I1Ca_m2j) ~ (k1m2*C1_m2j + k4m2*I2Ca_m2j + s1m2*Po_LCCj_m2 - (k2m2+k3m2+s2m2)*I1Ca_m2j),                                    # du[68]
        D(I2Ca_m2j) ~ (k3m2*I1Ca_m2j + k6m2* C2_m2j - (k4m2+k5m2)*I2Ca_m2j),                                                          # du[69]
        D(I1Ba_m2j) ~ (k1pm2*C1_m2j + k4pm2*I2Ba_m2j + s1pm2*Po_LCCj_m2 - (k2pm2+k3pm2+s2pm2)*I1Ba_m2j),                              # du[70]
        D(I2Ba_m2j) ~ (k3pm2*I1Ba_m2j + k6pm2* C2_m2j - (k5pm2+k4pm2)*I2Ba_m2j)                                                       # du[71]
    ]

    ibarca_jm2 = (4*pCa*Vm*Frdy*FoRT)*(.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1) 
    I_Ca_junc_m2 = (Fjunc_CaL*ibarca_jm2*(Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale 

    # CaMKII AND PKA-DEPENDENT SHIFTING OF DYADIC LCCS TO MODE 2
    fracLCCap0 = 0.219577  # Derived
    frac_fpkam2 = (0.15*fracLCCap0)/(1-fracLCCap0) 
    fpkam2 = (0.15+frac_fpkam2)*LCCa_PKAp - frac_fpkam2  # Assumes max (100#) phosphorylation results in 15# mode 2 channels
    fckiim2 = LCC_CKp*0.1  # Assumes max phosphorylation results in 10# mode 2 channels (max LCC_CKp = 1)
    # Sum up total fraction of CKII and PKA-shifted mode 2 channels
    junc_mode2 = fckiim2 + fpkam2 
    # Total junctional ICa
    I_Ca_junc2 = (1-junc_mode2)*I_Ca_junc_m1 + junc_mode2*I_Ca_junc_m2 

    # SUB-SARCOLEMMAL LCCs
    # Re-assign necessary params to be Casl sensitive
    fcasl = 1/(1+(kop/Ca_sl)^3)     # Depends on sl Ca
    TCasl = (tca + 0.1*(1+(Ca_sl/cpbar))^2)/(1+(Ca_sl/cpbar)^2) 
    tauCasl = (Rv-TCasl)*PrLCC + TCasl 
    # Re-assign necessary rates to be Casl sensitive
    s1sl = s1o*fcasl 
    k1sl = k1o*fcasl 
    s2sl = s1sl*(k2/k1sl)*(r1/r2) 
    s2psl = s1p*(k2p/k1p)*(r1/r2) 
    k5sl = (1-PsLCC)/tauCasl 
    k5sl = k5sl/recoveryReduc   # Reduced for recovery
    k6sl = (fcasl*PsLCC)/tauCasl 
    k4sl = k3*(alphaLCC/betaLCC)*(k1sl/k2)*(k5sl/k6sl) 
    k4psl = k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p) 

    # State transitions for 'mode 1' sarcolemmal LCCs
    Po_LCCsl_m1 = 1-C2_m1sl-C1_m1sl-I1Ca_m1sl-I2Ca_m1sl-I1Ba_m1sl-I2Ba_m1sl                                             # O_m1sl

    ICaMar_m1sl_eqs = [
        D(C2_m1sl) ~ (betaLCC*C1_m1sl + k5sl*I2Ca_m1sl + k5p*I2Ba_m1sl - (k6sl+k6p+alphaLCC)*C2_m1sl),                    # du[72]
        D(C1_m1sl) ~ (alphaLCC*C2_m1sl + k2*I1Ca_m1sl + k2p*I1Ba_m1sl + r2*Po_LCCsl_m1 - (r1+betaLCC+k1sl+k1p)*C1_m1sl),  # du[73]
        D(I1Ca_m1sl) ~ (k1sl*C1_m1sl + k4sl*I2Ca_m1sl + s1sl*Po_LCCsl_m1 - (k2+k3+s2sl)*I1Ca_m1sl),                       # du[74]
        D(I2Ca_m1sl) ~ (k3*I1Ca_m1sl + k6sl*C2_m1sl - (k4sl+k5sl)*I2Ca_m1sl),                                             # du[75]
        D(I1Ba_m1sl) ~ (k1p*C1_m1sl + k4psl*I2Ba_m1sl + s1p*Po_LCCsl_m1 - (k2p+k3p+s2psl)*I1Ba_m1sl),                     # du[76]
        D(I2Ba_m1sl) ~ (k3p*I1Ba_m1sl + k6p*C2_m1sl - (k5p+k4psl)*I2Ba_m1sl)                                              # du[77]
    ]

    ibarca_slm1 = (4*pCa*Vm*Frdy*FoRT)*(.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1) 
    I_Casl_m1 = (Fsl_CaL*ibarca_slm1*Po_LCCsl_m1*Q10CaL^Qpow)*ICa_scale 

    # Adjust closing rate for 'mode 2' sarcolemmal LCCs
    r2slm2 = r2m2 
    s2slm2 = s1sl*(k2/k1sl)*(r1/r2slm2) 
    s2pslm2 = s1p*(k2p/k1p)*(r1/r2slm2) 
    
    # State transitions for mode 2 sarcolemmal LCCs
    Po_LCCsl_m2 = 1-C2_m2sl-C1_m2sl-I1Ca_m2sl-I2Ca_m2sl-I1Ba_m2sl-I2Ba_m2sl                                                 # O_m2sl
    ibarca_slm2 = (4*pCa*Vm*Frdy*FoRT)*(.001*exp(2*Vm*FoRT)-0.341*Cao)/(exp(2*Vm*FoRT)-1) 
    I_Casl_m2 = (Fsl_CaL*ibarca_slm2*Po_LCCsl_m2*Q10CaL^Qpow)*ICa_scale 
    
    # Sum mode 1 and mode 2 sl channels for total sl current
    @parameters fckiim2_sl = 0  # Set to zero since SL LCCp by CaMKII is negligible
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
    I_Ca = I_Ca_junc+I_Ca_sl    # Total Ca curren throuhgh LCC
    I_CaNa_junc = (1-ICa_MarkovFlag)*(I_CaNa_junc1) + (ICa_MarkovFlag)*(I_CaNa_junc2) 
    I_CaNa_sl = (1-ICa_MarkovFlag)*(I_CaNa_sl1) + (ICa_MarkovFlag)*(I_CaNa_sl2)
    I_CaK = (1-ICa_MarkovFlag)*(I_CaK1) + ICa_MarkovFlag*(I_CaK2)   # Total K current through LCC
    
    ICaMar_m2sl_eqs = [
        D(C2_m2sl) ~ (betaLCC*C1_m2sl + k5sl*I2Ca_m2sl + k5p*I2Ba_m2sl - (k6sl+k6p+alphaLCC)*C2_m2sl),                        # du[78]
        D(C1_m2sl) ~ (alphaLCC*C2_m2sl + k2*I1Ca_m2sl + k2p*I1Ba_m2sl + r2slm2*Po_LCCsl_m2 - (r1+betaLCC+k1sl+k1p)*C1_m2sl),  # du[79]
        D(I1Ca_m2sl) ~ (k1sl*C1_m2sl + k4sl*I2Ca_m2sl + s1sl*Po_LCCsl_m2 - (k2+k3+s2slm2)*I1Ca_m2sl),                         # du[80]
        D(I2Ca_m2sl) ~ (k3*I1Ca_m2sl + k6sl*C2_m2sl - (k4sl+k5sl)*I2Ca_m2sl),                                                 # du[81]
        D(I1Ba_m2sl) ~ (k1p*C1_m2sl + k4psl*I2Ba_m2sl + s1p*Po_LCCsl_m2 - (k2p+k3p+s2pslm2)*I1Ba_m2sl),                       # du[82]
        D(I2Ba_m2sl) ~ (k3p*I1Ba_m2sl + k6p*C2_m2sl - (k5p+k4psl)*I2Ba_m2sl),                                                 # du[83]
        D(influx_LTCC) ~ (-I_Ca*Cmem/(Vmyo*2*Frdy)*1e3)                                                                       # du[43]
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
        D(influx_NCX) ~ (2*I_ncx*Cmem/(Vmyo*2*Frdy)*1e3)     # uM/ms, du[45]
    ]
    

    ## I_pca: Sarcolemmal Ca Pump Current / I_cabk: Ca Background Current
    # I_pca: Sarcolemmal Ca Pump Current
    I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*Ca_j^1.6/(KmPCa^1.6+Ca_j^1.6) 
    I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*Ca_sl^1.6/(KmPCa^1.6+Ca_sl^1.6) 
    I_pca = I_pca_junc+I_pca_sl 
    # I_cabk: Ca Background Current
    I_cabk_junc = Fjunc*GCaB*(Vm-eca_junc) 
    I_cabk_sl = Fsl*GCaB*(Vm-eca_sl) 
    I_cabk = I_cabk_junc+I_cabk_sl 

    ICabp_eqs = [
        D(influx_PMCA) ~ (-I_pca*Cmem/(Vmyo*2*Frdy)*1e3), # du[44]
        D(influx_ICa) ~ (-I_cabk*Cmem/(Vmyo*2*Frdy)*1e3) # du[46]
    ]


    ## I_CFTR or I_cl_(cAMP) - Cystic Fibrosis Transmembrane Conductance Reg.
    @parameters Icftr = 0       # gCFTR*(y(39) - ecl)           # NO Icftr in MOUSE
    

    ## RyR model - SR release fluxes and leak
    # CaMKII and PKA-dependent phosphoregulation of RyR Po
    fCKII_ec50SR = 1.16 - 4/5*RyR_CKp 
    ec50SR = fCKII_ec50SR*ec50SR  # MOUSE - 60# 

    @parameters MaxSR = 15  
    @parameters MinSR = 1 
    kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/Ca_sr)^2.5) 
    koSRCa = koCa/kCaSR 
    kiSRCa = kiCa*kCaSR 
    @parameters kleak = 2*5.348e-6  # [1/ms] changed from rabbit (5.348e-6)

    #fCKII_RyR = (20*RyR_CKp/3 - 1/3)  # 1 at basal condition - RABBIT
    fCKII_RyR = (10*RyR_CKp - 1)  # 1 at basal condition - MOUSE

    #fPKA_RyR = RyR_PKAp*1.025 + 0.9750  # 1 with NO ISO
    frac_RyRo = 0.204276  # Derived (RyR_PKAp(basal)/RyRtot)
    a_RyR = (2-1)/(1/frac_RyRo-1)  # Max effect: fPKA_RyR=2
    fPKA_RyR = 1-a_RyR+a_RyR*(RyR_PKAp/frac_RyRo) 
    koSRCa = (fCKII_RyR + fPKA_RyR - 1)*koSRCa 

    # ODEs for RyR states and SR release through open RyRs
    RI = 1-RyR_R-RyR_O-RyR_I 
    RyR_eqs = [
        D(RyR_R) ~ ((kim*RI-kiSRCa*Ca_j*RyR_R)-(koSRCa*Ca_j^2*RyR_R-kom*RyR_O)),      # R du[14]
        D(RyR_O) ~ ((koSRCa*Ca_j^2*RyR_R-kom*RyR_O)-(kiSRCa*Ca_j*RyR_O-kim*RyR_I)),   # O du[15]
        D(RyR_I) ~ ((kiSRCa*Ca_j*RyR_O-kim*RyR_I)-(kom*RyR_I-koSRCa*Ca_j^2*RI))       # I du[16]
    ]
    J_SRCarel = ks*RyR_O*(Ca_sr-Ca_j)           # [mmol/L SR/ ms]

    # Passive RyR leak - includes CaMKII regulation of leak flux
    #kleak = (1/3 + 10*RyR_CKp/3)*kleak  # RABBIT
    kleak = (1/2 + 5*RyR_CKp/2)*kleak  # MOUSE (reduced CaMKII effect on leak)
    J_SRleak = kleak*(Ca_sr-Ca_j)  # [mmol/L cyt/ms]


    ## SERCA model - SR uptake fluxes
    # CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
    fCKII_PLB = (1-.5*PLB_CKp)  # Max effect: fCKII_PLB=0.5
    @parameters fracPKA_PLBo = 1-0.079755
    fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*(100-55.31)/100 + 55.31/100 

    # Select smaller value (resulting in max reduction of Kmf)
    Kmf = ifelse(fCKII_PLB < fPKA_PLB, Kmf*fCKII_PLB, Kmf*fPKA_PLB)     #fCKII_PLB
    J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*((Cai/Kmf)^hillSRCaP-(Ca_sr/Kmr)^hillSRCaP)/(1+(Cai/Kmf)^hillSRCaP+(Ca_sr/Kmr)^hillSRCaP) # [mM/msec] 


    ## Na and Ca Buffering
    Buffering_eqs = [
        D(NaBj) ~ (kon_na*Naj*(Bmax_Naj-NaBj)-koff_na*NaBj),                                      # du[17],  [mM/ms]      
        D(NaBsl) ~ (kon_na*Nasl*(Bmax_Nasl-NaBsl)-koff_na*NaBsl),                                 # du[18],  [mM/ms]      
        # Cytosolic Ca Buffers0
        D(TnCL) ~ (kon_tncl*Cai*(Bmax_TnClow-TnCL)-koff_tncl*TnCL),                               # du[19],  [mM/ms]      
        D(TnCHc) ~ (kon_tnchca*Cai*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchca*TnCHc),                 # du[20],  [mM/ms]      
        D(TnCHm) ~ (kon_tnchmg*Mgi*(Bmax_TnChigh-TnCHc-TnCHm)-koff_tnchmg*TnCHm),                 # du[21],  [mM/ms]      
        D(CaM) ~ 0,                                                                             # du[22],  [mM/ms]
        D(Myosin_ca) ~ (kon_myoca*Cai*(Bmax_myosin-Myosin_ca-Myosin_mg)-koff_myoca*Myosin_ca),    # du[23],  [mM/ms]      
        D(Myosin_mg) ~ (kon_myomg*Mgi*(Bmax_myosin-Myosin_ca-Myosin_mg)-koff_myomg*Myosin_mg),    # du[24],  [mM/ms]      
        D(SRB) ~ (kon_sr*Cai*(Bmax_SR-SRB)-koff_sr*SRB),                                          # du[25],  [mM/ms]
        # Junctional and SL Ca Buffers
        D(SLLj) ~ (kon_sll*Ca_j*(Bmax_SLlowj-SLLj)-koff_sll*SLLj),                                # du[26],  [mM/ms]      
        D(SLLsl) ~ (kon_sll*Ca_sl*(Bmax_SLlowsl-SLLsl)-koff_sll*SLLsl),                           # du[27],  [mM/ms]      
        D(SLHj) ~ (kon_slh*Ca_j*(Bmax_SLhighj-SLHj)-koff_slh*SLHj),                               # du[28],  [mM/ms]      
        D(SLHsl) ~ (kon_slh*Ca_sl*(Bmax_SLhighsl-SLHsl)-koff_slh*SLHsl)                           # du[29],  [mM/ms]    
    ]
    J_CaB_cytosol = D(TnCL)+D(TnCHc)+D(CaM)+D(Myosin_ca)+D(SRB)
    J_CaB_junction = D(SLLj)+D(SLHj)
    J_CaB_sl = D(SLLsl)+D(SLHsl)
    

    ## Ion concentrations

    @variables JCaDyad(t) JCaSL(t) JCaCyt(t)
    # Na Concentrations
    I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc    # [uA/uF]
    I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl    #[uA/uF]
    # K Concentration
    I_K_tot = I_to+I_kr+I_ks+I_k1-2*I_nak+I_CaK+I_kp+I_kur+I_ss      # [uA/uF]
    # Ca Concentrations
    I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc  # [uA/uF]
    I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl            # [uA/uF]

    Ion_eqs = [
        D(Csqn) ~ (kon_csqn*Ca_sr*(Bmax_Csqn-Csqn)-koff_csqn*Csqn),
        D(Ca_sr) ~ (J_serca*Vmyo/Vsr-(J_SRleak*Vmyo/Vsr+J_SRCarel))-D(Csqn),
        D(Naj) ~ (-I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(Nasl-Naj))-D(NaBj),
        D(Nasl) ~ (-I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(Naj-Nasl)+J_na_slmyo/Vsl*(Nai-Nasl))-D(NaBsl),
        D(Nai) ~ ifelse(NaClampFlag == 1, 0, J_na_slmyo/Vmyo*(Nasl-Nai)),
        D(Ki) ~ -D(Ca_j) + 1e-3*JCaDyad,
        D(Ca_j) ~ (-I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(Ca_sl-Ca_j)-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc),
        D(Ca_sl) ~ (-I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(Ca_j-Ca_sl)+ J_ca_slmyo/Vsl*(Cai-Ca_sl)-J_CaB_sl+ 1e-3*JCaSL),
        D(Cai) ~ (-J_serca-J_CaB_cytosol+J_ca_slmyo/Vmyo*(Ca_sl-Cai)+1e-3*JCaCyt)
    ]

    #=
    Ion_eqs = [
        D(Csqn) ~ (kon_csqn*Ca_sr*(Bmax_Csqn-Csqn)-koff_csqn*Csqn),                                                                               # du[30], [mM/ms]
        D(Ca_sr) ~ (J_serca*Vmyo/Vsr-(J_SRleak*Vmyo/Vsr+J_SRCarel))-D(Csqn),                                                                      # du[31], [mM/ms]
        D(Naj) ~ (-I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(Nasl-Naj))-D(NaBj),                                                         # du[32] 
        D(Nasl) ~ (-I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(Naj-Nasl)+J_na_slmyo/Vsl*(Nai-Nasl))-D(NaBsl),                                   # du[33] 
        D(Nai) ~ ifelse(NaClampFlag == 1, 0, J_na_slmyo/Vmyo*(Nasl-Nai)),                                                                       # du[34], [mM/msec]
        D(Ki) ~ -D(Ca_j) + 1e-3*JCaDyad,                                                                                                        # du[35], [mM/msec]
        D(Ca_j) ~ (-I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(Ca_sl-Ca_j)-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc),   # du[36]
        D(Ca_sl) ~ (-I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(Ca_j-Ca_sl)+ J_ca_slmyo/Vsl*(Cai-Ca_sl)-J_CaB_sl+ 1e-3*JCaSL),                # du[37]   
        D(Cai) ~ (-J_serca-J_CaB_cytosol+J_ca_slmyo/Vmyo*(Ca_sl-Cai)+1e-3*JCaCyt)                                                                 # du[38]
    ]
=#

    # AP Waveform for AP clamp

    if protocol == 1
        I_app = ifelse(mod(t,cycleLength) <= 5, 9.5, 0.0)
    elseif protocol == 2 
        @parameters V_clamp = -83   #y(39)#-80
        @parameters R_clamp = 0.01
        I_app = (V_clamp-Vm)/R_clamp
    else
        I_app = 0.0
    end          

    ## Membrane Potential

    I_Na_tot = I_Na_tot_junc + I_Na_tot_sl                  # [uA/uF]
    I_Cl_tot = I_ClCa+I_Clbk+Icftr                          # [uA/uF]
    I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl                    # [uA/uF]
    I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot              # [uA/uF]

    MemPot_eqs = [
        D(Vm) ~ -(I_tot-I_app)                                  # du[39]
    ]
    
    #eqs = append!(var_eqs, INa_fast_eqs, INa_eqs, IK_eqs, Ito_eqs, ICa_HH_eqs, ICaMar_m1j_eqs, 
    #            ICaMar_m2j_eqs, ICaMar_m1sl_eqs, ICaMar_m2sl_eqs, I_ncx_eqs, ICabp_eqs, RyR_eqs, Buffering_eqs, Ion_eqs, AP_eqs, MemPot_eqs)

    #@named sys = ODESystem(eqs, t)

    #return sys

    return vcat(INa_fast_eqs, INa_eqs, IK_eqs, Ito_eqs, ICa_HH_eqs, ICaMar_m1j_eqs, 
    ICaMar_m2j_eqs, ICaMar_m1sl_eqs, ICaMar_m2sl_eqs, I_ncx_eqs, ICabp_eqs, RyR_eqs, Buffering_eqs, Ion_eqs, MemPot_eqs)
end



