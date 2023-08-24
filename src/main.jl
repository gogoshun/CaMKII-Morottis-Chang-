using Catalyst
using DifferentialEquations
using Plots
using ModelingToolkit


include("camdyad_ODEfile(M).jl")
include("camsl_ODEfile(M).jl")
include("camcyt_ODEfile(M).jl")
include("camkii_ODEfile(M).jl")
include("ecc_debug.jl")
include("bar_ODEfile(M).jl")


eq_camdyad = get_camdyad_equations()
eq_camsl = get_camsl_equations()
eq_camcyt = get_camcyt_equations()
eq_camkii = get_camkii_equations()
eq_ecc =  get_ecc_equations()
eq_bar = get_bar_equations()

@named osys = ODESystem(vcat(eq_camdyad, eq_camcyt, eq_camsl, eq_camkii, eq_ecc, eq_bar))

osys = structural_simplify(osys)

@variables t Cai(t)

#Chemical Reaction
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

@variables t J_cam_dyadSL(t)
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
        # Intermediate Variables // dyad
        JCaDyad, CaMtotDyad, Bdyad, J_cam_dyadSL, J_ca2cam_dyadSL, J_ca4cam_dyadSL,
        # Intermediate Variables // sl
        JCaSL,
        # Intermediate Variables // cyt
        JCaCyt, J_cam_SLmyo, J_ca2cam_SLmyo, J_ca4cam_SLmyo,
        # CaMKII Model
        CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, 
        Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2, k_1C_on, k_1C_off, k_2C_on, k_2C_off,
        k_1N_on, k_1N_off, k_2N_on, k_2N_off, k_K1C_on, k_K1C_off, k_K2C_on, k_K2C_off, k_K1N_on, k_K1N_off,
        k_K2N_on, k_K2N_off, kCaM0_on, kCaM2C_on, kCaM2N_on, kCaM4_on, kCaM0_off, kCaM2C_off, kCaM2N_off, kCaM4_off, 
        kCaM0P_on, kCaM2CP_on, kCaM2NP_on, kCaM4P_on, kCaM0P_off, kCaM2CP_off, kCaM2NP_off, kCaM4P_off, k_phosCaM, 
        k_dephospho, k_P1_P2, k_P2_P1, CaMKII_T = sys

oprob = ODEProblem(sys, [CaM0 => 2.82e-5, Ca2CaM_C => 1.01e-8, Ca2CaM_N => 1.40e-9, Ca4CaM => 4.78e-13,
        CaM0_CaMK => 1.29e-6, Ca2CaM_C_CaMK => 9.13e-8, Ca2CaM_N_CaMK => 3.74e-9, Ca4CaM_CaMK => 5.92e-10,
        CaM0_CaMKP => 2.36e-7, Ca2CaM_C_CaMKP => 1.13e-7, Ca2CaM_N_CaMKP => 1.54e-9, Ca4CaM_CaMKP => 7.82e-10,
        CaMK => 6.73e-5, CaMKP => 6.57e-7, CaMKP2 => 2.66e-7,
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
        #dyad
        JCaDyad => 0.0, CaMtotDyad => 0.0, Bdyad => 0.0, J_cam_dyadSL => 0.0, J_ca2cam_dyadSL => 0.0, J_ca4cam_dyadSL => 0.0,
        #sl
        JCaSL => 0.0,
        #cyt
        JCaCyt => 0.0, J_cam_SLmyo => 0.0, J_ca2cam_SLmyo => 0.0, J_ca4cam_SLmyo => 0.0,
        CaM0 => 2.82e-5, Ca2CaM_C => 1.01e-8, Ca2CaM_N => 1.40e-9, Ca4CaM => 4.78e-13,
        CaM0_CaMK => 1.29e-6, Ca2CaM_C_CaMK => 9.13e-8, Ca2CaM_N_CaMK => 3.74e-9, Ca4CaM_CaMK => 5.92e-10,
        CaM0_CaMKP => 2.36e-7, Ca2CaM_C_CaMKP => 1.13e-7, Ca2CaM_N_CaMKP => 1.54e-9, Ca4CaM_CaMKP => 7.82e-10,
        CaMK => 6.73e-5, CaMKP => 6.57e-7, CaMKP2 => 2.66e-7], 
        (0.0, 500000.0), 
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

sol = solve(oprob, Rodas5(), progress=true)

plot(sol, idxs=Cai, linewidth=1.5, title="Calcium Transient", xlabel="Time(s)", ylabel="[Ca2+](M)", label="Control")

plot(sol, idxs=[CaM0_CaMK + Ca2CaM_C_CaMK+Ca2CaM_N_CaMK+Ca4CaM_CaMK+CaM0_CaMKP+Ca2CaM_C_CaMKP+Ca2CaM_N_CaMKP+Ca4CaM_CaMKP+CaMKP+CaMKP2], linewidth=3, xlabel="Time(s)", ylabel="Concentration(M)", title="CaMKII Activity", xlim=(190,320), label="activeCaMKII")
