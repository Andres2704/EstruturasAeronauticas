import structural as st
import numpy as np 
from visualization import * 

# Dados de entrada da aeronave e da seção transversal
dados = dict(
    # Dados da aeronave
    Vd          = 72.02,    # Velocidade de mergulho, se 0 calculada com Vc       
    Vc          = 0.0,      # Velocidade de cruzeiro, se 0 calculada com Vd
    rho         = 1.225,    # Densidade na altitude de estudo
    rho0        = 1.225,    # Densidade de referencia 
    CL_max      = 1.46,     # Coef. de sustentação máximo
    CL_min      = -0.45,    # Coef. de sustentação mínimo
    factor      = 1.25,     # Fator de correção da velocidade de stall
    CL_alpha    = 2.8060,   # Derivada de estabilidade CL vs alpha
    CMA         = 1.5,      # Corda média aerodinamica
    W           = 5873,     # Peso da aeronave [N]
    S           = 10.74,    # Área de referencia [m²]
    n_max       = 4,        # Fator de carga máximo
    n_min       = -2.67,    # Fator de cara mínimo
    lamb        = 0.875,    # Afilamento
    env         = 7.12,     # Envergadura [m]
    npoints     = 100,      # Número de pontos para as simulações
    
    # Dados da seção 
    type        = "H",                      # Tipo de viga (H do nosso projeto)
    N_cel       = 2,                        # Número de células da seção 
    E           = [110e9, 70e9, 200e9],     # Módulo de elasticidade de cada material
    G           = [44e9, 27e9, 4e9],        # Módulo de torção 
    t           = [0.02, 0.02, 0.02],       # Espessura da alma, aba e revestimento 
    h1          = 0.3,                      # Altura da alma
    l           = 0.3,                      # Comprimento total da caixa
    b_por       = 80,                       # Porcentagem do comprimento ocupado com aba
    i_ref_mat   = 0,                        # Material de referencia E0 e G0 (0, 1 ou 2)
)

# Dados de entrada nos pontos críticos do diagrama V-n
dados_crit = dict(
    CL_PHAA = 1.825,
    CL_NHAA = -0.450,
    alpha_PHAA = 30*np.pi/180,
    alpha_PLAA = 6*np.pi/180,
    alpha_NLAA = -9*np.pi/180,
    alpha_NHAA = -9*np.pi/180,
    cd_PHAA = 0.292,
    cd_PLAA = 0.037, 
    cd_NHAA = 0.023, 
    cd_NLAA = 0.023,
    Cm_PHAA = -0.428,           # Utilizado apenas com momento uniforme
    Cm_PLAA = -0.230,           # Utilizado apenas com momento uniforme
    Cm_NHAA = 0.043,            # Utilizado apenas com momento uniforme
    Cm_NLAA = 0.043             # Utilizado apenas com momento uniforme
)

# Construindo o objeto estrutura e torcao 
estrutura = st.aero_struct(**dados)
torcao = st.torsion(**dados)

# Propagando a solucao nos pontos críticos do diagram V-n
estrutura.propagate_aero_critic(**dados_crit)
torcao.propagate_torsion_crit(
    type = 'fit',           # Fit = curva Ycp vs Cl | Uniform = momento uniforme
    c = 1.6,                    # Comprimento de referencia
    **dados_crit)

# Plot do esforço cortante
plot_subplot(estrutura.x_coord, 
             estrutura.Vy_PHAA, estrutura.Vy_PLAA, estrutura.Vy_NHAA, estrutura.Vy_NLAA,
             estrutura.Vz_PHAA, estrutura.Vz_PLAA, estrutura.Vz_NHAA, estrutura.Vz_NLAA,
             f = 1000, 
             title = 'Esforço cortante', 
             ylabel = [r"$V_y$ (x) [kN]", r"$V_z$ (x) [kN]"],
             xlabel = ["x [m]", "x [m]"],
             label = ["PHAA", "PLAA", "NHAA", "NLAA"])

# Plot do momento fletor
plot_subplot(estrutura.x_coord, 
             estrutura.My_PHAA, estrutura.My_PLAA, estrutura.My_NHAA, estrutura.My_NLAA,
             estrutura.Mz_PHAA, estrutura.Mz_PLAA, estrutura.Mz_NHAA, estrutura.Mz_NLAA,
             f = 1000, 
             title = 'Momento Fletor', 
             ylabel = [r"$M_y$ (x) [kN.m]", r"$M_z$ (x) [kN.m]"],
             xlabel = ["x [m]", "x [m]"],
             label = ["PHAA", "PLAA", "NHAA", "NLAA"])

# Plot das deflexoes
plot_subplot(estrutura.x_coord, 
             estrutura.v_def_PHAA, estrutura.v_def_PLAA, estrutura.v_def_NHAA, estrutura.v_def_NLAA,
             estrutura.w_def_PHAA, estrutura.w_def_PLAA, estrutura.w_def_NHAA, estrutura.w_def_NLAA,
             f = 1/1000, 
             ylabel = [r"$v$ (x) [mm]", r"$w$ (x) [mm]"],
             xlabel = ["x [m]", "x [m]"],
             label = ["PHAA", "PLAA", "NHAA", "NLAA"])

# Plot do momento de arfagem e torque torsor interno
plot_subplot(estrutura.x_coord, 
             torcao.mx_PHAA, torcao.mx_PLAA, torcao.mx_NHAA, torcao.mx_NLAA,
             torcao.Tx_PHAA, torcao.Tx_PLAA, torcao.Tx_NHAA, torcao.Tx_NLAA,
             f = 1000, 
             title = 'Momento de arfagem distribuido e momento torsor interno',
             ylabel = [r"$m_x$ (x) [kN/m]", r"$T_x$ (x) [kN.m]"],
             xlabel = ["x [m]", "x [m]"],
             label = ["PHAA", "PLAA", "NHAA", "NLAA"])

# Plot do ângulo de deformação 
plot_plot(estrutura.x_coord, 
             torcao.theta_PHAA, torcao.theta_PLAA, torcao.theta_NHAA, torcao.theta_NLAA,
             f = np.pi/180, 
             title = 'Ângulo de torção',
             ylabel = r"$\theta$ (x) [º]",
             xlabel = "x [m]",
             label = ["PHAA", "PLAA", "NHAA", "NLAA"])

plot_sigma(estrutura.sigma_PHAA, l = 0.3, htot = estrutura.htot, 
           i_ref_plot= 0, 
           centroide=[estrutura.ys, estrutura.zs], 
           alpha = estrutura.alpha_axis_PHAA[0], 
           title = "Campo de Tensões σₓ(y,z) - PHAA" )

plot_sigma(estrutura.sigma_PLAA, l = 0.3, htot = estrutura.htot, 
           i_ref_plot= 0, 
           centroide=[estrutura.ys, estrutura.zs], 
           alpha = estrutura.alpha_axis_PLAA[0], 
           title = "Campo de Tensões σₓ(y,z) - PLAA" )

plot_sigma(estrutura.sigma_NHAA, l = 0.3, htot = estrutura.htot, 
           i_ref_plot= 0, 
           centroide=[estrutura.ys, estrutura.zs], 
           alpha = estrutura.alpha_axis_NHAA[0], 
           title = "Campo de Tensões σₓ(y,z) - NHAA" )

plot_sigma(estrutura.sigma_NLAA, l = 0.3, htot = estrutura.htot, 
           i_ref_plot= 0, 
           centroide=[estrutura.ys, estrutura.zs], 
           alpha = estrutura.alpha_axis_NLAA[0], 
           title = "Campo de Tensões σₓ(y,z) - NLAA" )

plot_qfield(torcao.q_field_PHAA, q_val = torcao.q_PHAA[:, 0], 
           l = 0.3, htot = estrutura.htot, i_ref_plot= 0,
           centroide=[estrutura.ys, estrutura.zs], 
           title = "Campo de fluxo de cisalhamento q(y,z) - PHAA" )

plot_qfield(torcao.q_field_PLAA, q_val = torcao.q_PLAA[:, 0], 
            l = 0.3, htot = estrutura.htot, i_ref_plot= 0,
           centroide=[estrutura.ys, estrutura.zs], 
           title = "Campo de fluxo de cisalhamento q(y,z) - PLAA" )

plot_qfield(torcao.q_field_NHAA, q_val = torcao.q_NHAA[:, 0], 
            l = 0.3, htot = estrutura.htot, i_ref_plot= 0,
           centroide=[estrutura.ys, estrutura.zs], 
           title = "Campo de fluxo de cisalhamento q(y,z) - NHAA" )

plot_qfield(torcao.q_field_NLAA, q_val = torcao.q_NLAA[:, 0], 
            l = 0.3, htot = estrutura.htot, i_ref_plot= 0,
           centroide=[estrutura.ys, estrutura.zs], 
           title = "Campo de fluxo de cisalhamento q(y,z) - NLAA" )



plt.show()