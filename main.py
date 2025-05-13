from func_solver import * 
 
if __name__ == "__main__":
    # ================ Dados de entrada ===============
    Vd = 259.3/3.6 # Velocidade de cruzeiro nunca ultrapassada
    rho = 1.225                                         # Densidade na altura desejada [kg/m³]
    Cl_max = 1.46                                       # Máximo coef de sustenação
    Cl_min = -0.45                                      # Min coef de sustentação
    AOA_clmax = 30*np.pi/180                            # AOA de estol positivo
    AOA_clmin = -9*np.pi/180                            # AOA de estol negativo
    Clalpha = (Cl_max - Cl_min)/(AOA_clmax - AOA_clmin) # Clalpha
    W = 598.74*9.81                                     # Peso da aeronave [N]
    S = 10.74                                           # Área de ref [m²]
    n_max = 6/1.5
    n_min = -4/1.5

    # ============================= Plot do diagrama V-n ===========================
    [V_NHAA, V_PHAA] = plot_vn(rho, W, S,Vd = Vd, Clmax = Cl_max, Clmin = Cl_min, Cla = Clalpha, c_bar = 1.5, n_max=n_max, n_min=n_min, plot = False)

    # ===================================== Cálculo da sustentação =================================
    L_max = n_max*W                         # Sustentação max
    L_min = n_min*W                         # Sustentação min
    Cl_PHAA = 1.825                         # Cl no ponto PHAA
    Cl_PLAA = 2*L_max/(rho*Vd**2 * S)       # Cl no ponto PLAA
    Cl_NHAA = -0.450                        # Cl no ponto NHAA
    Cl_NLAA = 2*L_min/(rho*Vd**2 * S)       # Cl no ponto NLAA
    B = 8.23                                # Envergadura + cabine [m]
    C = 1.11                                # Envergadura da Cabine [m]
    b = B - C                               # Envergadura da asa sem a cabine (para visualização da semi-asa) [m]
    lamb = 0.875

    alpha_NHAA = -9*np.pi/180
    alpha_PHAA = 30*np.pi/180
    alpha_NLAA = -9*np.pi/180
    alpha_PLAA = 6*np.pi/180

    print('======================= SUSTENTACAO ===================')
    print('Cl_PLAA = ', Cl_PLAA)
    print('Cl_PHAA = ', Cl_PHAA)
    print('Cl_NLAA = ', Cl_NLAA)
    print('Cl_NHAA = ', Cl_NHAA)
    print('Sustentação PHAA = ', L_max, 'N')
    print('Sustentação PLAA = ', L_max, 'N')
    print('Sustentação NHAA = ', L_min, 'N')
    print('Sustentação NLAA = ', L_min, 'N')

    [x, w_sust_pos] = lift_Schrenk(L_max, b, lamb, n = 100)
    [x, w_sust_neg] = lift_Schrenk(L_min, b, lamb, n = 100)


    # ===================================== CÁLCULO DO ARRASTO ===========================================
    Cd_PHAA = 0.292                           # Coeficiente de arrasto
    Cd_NHAA = 0.023                           # Coeficiente de arrasto
    Cd_PLAA = 0.037                           # Coeficiente de arrasto
    Cd_NLAA = 0.023                           # Coeficiente de arrasto

    w_arrasto_PHAA = drag_linear(rho, V_PHAA, S, Cd_PHAA, b, n = len(x))
    w_arrasto_PLAA = drag_linear(rho, Vd, S, Cd_PLAA, b, n = len(x))
    w_arrasto_NHAA = drag_linear(rho, V_NHAA, S, Cd_NHAA, b, n = len(x))
    w_arrasto_NLAA = drag_linear(rho, Vd, S, Cd_NLAA, b, n = len(x))

    arrasto_total_PHAA = drag_total(rho, V_PHAA, S, Cd_PHAA)
    arrasto_total_PLAA = drag_total(rho, Vd, S, Cd_PLAA)
    arrasto_total_NHAA = drag_total(rho, V_NHAA, S, Cd_NHAA)
    arrasto_total_NLAA = drag_total(rho, Vd, S, Cd_NLAA)

    # ====================================  PLOT DAS FORÇAS AERODINAMICAS NO REF AERO ==============================
    F_distrib_PHAA = np.vstack((w_arrasto_PHAA, w_sust_pos))
    F_distrib_PLAA = np.vstack((w_arrasto_PLAA, w_sust_pos))
    F_distrib_NHAA = np.vstack((w_arrasto_NHAA, w_sust_neg))
    F_distrib_NLAA = np.vstack((w_arrasto_NLAA, w_sust_neg))

    plot_3d_forces(x, F_distrib_PHAA, L_max, arrasto_total_PHAA, title="Distribuição de Forças Aerodinâmicas PHAA")
    plot_3d_forces(x, F_distrib_PLAA, L_max, arrasto_total_PLAA, title="Distribuição de Forças Aerodinâmicas PLAA")
    plot_3d_forces(x, F_distrib_NHAA, L_min, arrasto_total_NHAA, title="Distribuição de Forças Aerodinâmicas NHAA")
    plot_3d_forces(x, F_distrib_NLAA, L_min, arrasto_total_NLAA, title="Distribuição de Forças Aerodinâmicas NLAA")

    print('======================= ARRASTO ===================')
    print('Cd_PHAA = ', Cd_PHAA)
    print('Cd_PLAA = ', Cd_PLAA)
    print('Cd_NHAA = ', Cd_NHAA)
    print('Cd_NLAA = ', Cd_NLAA)
    print('Arrasto PHAA =', arrasto_total_PHAA, 'N')
    print('Arrasto PLAA =', arrasto_total_PLAA, 'N')
    print('Arrasto NHAA =', arrasto_total_NHAA, 'N')
    print('Arrasto NLAA =', arrasto_total_NLAA, 'N')

    print('======================= ANGULOS ===================')
    print('alpha_PHAA = ', alpha_PHAA*np.pi/180, 'Deg')
    print('alpha_PLAA = ', alpha_PLAA*np.pi/180, 'Deg')
    print('alpha_NHAA = ', alpha_NHAA*np.pi/180, 'Deg')
    print('alpha_NLAA = ', alpha_NLAA*np.pi/180, 'Deg')

    # ===================================================================================================
    # =================================== --------------------- =========================================
    # =================================== | ANÁLISE DE FLEXAO | =========================================
    # =================================== --------------------- =========================================
    # ===================================================================================================
    
    # ====================================== CARGAS NO REF DA ESTRUTURA =================================
    F_PHAA = transf_aero2struct(alpha_PHAA, np.concatenate((w_arrasto_PHAA, w_sust_pos)).reshape(2, len(x)))
    F_PLAA = transf_aero2struct(alpha_PLAA, np.concatenate((w_arrasto_PLAA, w_sust_pos)).reshape(2, len(x)))
    F_NHAA = transf_aero2struct(alpha_NHAA, np.concatenate((w_arrasto_NHAA, w_sust_neg)).reshape(2, len(x)))
    F_NLAA = transf_aero2struct(alpha_NLAA, np.concatenate((w_arrasto_NLAA, w_sust_neg)).reshape(2, len(x)))

    # ========================= PROPRIEDADES DE ÁREA PARA VIGA T DE N CELULAS ===========================
    N_cel = 2                                   # Número de células
    E1, E2, E3 = 110, 70, 200                   # Módulo E [GPa] (alma Ti-6Al-4V, aba Al 7075-T6 e revestimento fibra de carbono)
    E0 = E1                                     # VERIFICAR SE ISSO É OK
    t1, t2, t3 = 20/1000, 20/1000, 20/1000      # Espessuras de cada material [m]
    h1 = 0.3                                    # Altura da alma [m]
    l = 0.3                                     # Largura total da caixa de asa [m]
    b = 0.8*l/(N_cel + 1)                       # Largura da aba [m] | A soma das abas é 80% da largura total
    e = (l - (N_cel+1)*b)/(N_cel + 2)           # Espaçamento entre as células
    [As, zs, ys, Izz_s, Iyy_s] = calc_area_prop(N_cel, [E1, E2, E3], [t1, t2, t3], h1, b, e)
    desenhar_vigas_T(N_cel+1, t1, t2, t3, b, h1, e, ponto_centroide=(ys, zs), alpha=None) 

    print('======================= PROPRIEDADES DE ÁREA ===================')
    print('Area Ponderada = ', As, 'm²')
    print(f'Posição do Centroide = {zs}, {ys:.2f} [z (m), y (m)]')
    print('Momento de Inercia Ponderado Izz  = ', Izz_s, 'kg.m²')
    print('Momento de Inercia Ponderado Iyy = ', Iyy_s, 'kg.m²')

    # ====================================== ESFORÇOS CORTANTES =========================================
    Vy_PLAA = shear_force(x, -F_PLAA[0, :])
    Vy_PHAA = shear_force(x, -F_PHAA[0, :])
    Vy_NLAA = shear_force(x, -F_NLAA[0, :])
    Vy_NHAA = shear_force(x, -F_NHAA[0, :])
    Vz_PLAA = shear_force(x, -F_PLAA[1, :])
    Vz_PHAA = shear_force(x, -F_PHAA[1, :])
    Vz_NLAA = shear_force(x, -F_NLAA[1, :])
    Vz_NHAA = shear_force(x, -F_NHAA[1, :])

    plt.figure(dpi = 100, figsize = (13,4))
    plt.subplot(1,2,1)
    plt.plot(x, Vz_PHAA/1000, color = 'r', label = 'PHAA')
    plt.plot(x, Vz_PLAA/1000, color = 'b', label = 'PLAA')
    plt.plot(x, Vz_NHAA/1000, color = 'k', label = 'NHAA')
    plt.plot(x, Vz_NLAA/1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$V_z$ (x) [kN]')
    plt.xlabel('x[m]')

    plt.subplot(1,2,2)
    plt.plot(x, Vy_PHAA/1000, color = 'r', label = 'PHAA')
    plt.plot(x, Vy_PLAA/1000, color = 'b', label = 'PLAA')
    plt.plot(x, Vy_NHAA/1000, color = 'k', label = 'NHAA')
    plt.plot(x, Vy_NLAA/1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$V_y$ (x) [kN]')
    plt.xlabel('x[m]')

    # ====================================== MOMENTO FLETOR ============================================
    My_PLAA = shear_force(x, -Vy_PLAA)
    My_PHAA = shear_force(x, -Vy_PHAA)
    My_NLAA = shear_force(x, -Vy_NLAA)
    My_NHAA = shear_force(x, -Vy_NHAA)
    Mz_PLAA = shear_force(x, Vz_PLAA)
    Mz_PHAA = shear_force(x, Vz_PHAA)
    Mz_NLAA = shear_force(x, Vz_NLAA)
    Mz_NHAA = shear_force(x, Vz_NHAA)

    plt.figure(dpi = 100, figsize = (13,4))
    plt.subplot(1,2,1)
    plt.plot(x, Mz_PHAA/1000, color = 'r', label = 'PHAA')
    plt.plot(x, Mz_PLAA/1000, color = 'b', label = 'PLAA')
    plt.plot(x, Mz_NHAA/1000, color = 'k', label = 'NHAA')
    plt.plot(x, Mz_NLAA/1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$M_z$ (x) [kN.m]')
    plt.xlabel('x[m]')

    plt.subplot(1,2,2)
    plt.plot(x, My_PHAA/1000, color = 'r', label = 'PHAA')
    plt.plot(x, My_PLAA/1000, color = 'b', label = 'PLAA')
    plt.plot(x, My_NHAA/1000, color = 'k', label = 'NHAA')
    plt.plot(x, My_NLAA/1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$M_y$ (x) [kN.m]')
    plt.xlabel('x[m]')
    


    # ======================================== DEFLEXÕES ==========================================
    v_def_NHAA = integrate_deflection(x, Mz_NHAA/(E0*1e9*Izz_s)) # 1e9 para transformar de GPa para Pa
    v_def_NLAA = integrate_deflection(x, Mz_NLAA/(E0*1e9*Izz_s))
    v_def_PHAA = integrate_deflection(x, Mz_PHAA/(E0*1e9*Izz_s))
    v_def_PLAA = integrate_deflection(x, Mz_PLAA/(E0*1e9*Izz_s))

    w_def_NHAA = integrate_deflection(x, -My_NHAA/(E0*1e9*Iyy_s))
    w_def_NLAA = integrate_deflection(x, -My_NLAA/(E0*1e9*Iyy_s))
    w_def_PHAA = integrate_deflection(x, -My_PHAA/(E0*1e9*Iyy_s))
    w_def_PLAA = integrate_deflection(x, -My_PLAA/(E0*1e9*Iyy_s))

    plt.figure(dpi = 100, figsize = (13,4))
    plt.subplot(1,2,1)
    plt.plot(x, v_def_PHAA*1000, color = 'r', label = 'PHAA')
    plt.plot(x, v_def_PLAA*1000, color = 'b', label = 'PLAA')
    plt.plot(x, v_def_NHAA*1000, color = 'k', label = 'NHAA')
    plt.plot(x, v_def_NLAA*1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$v$ (x) [mm]')
    plt.xlabel('x[m]')

    plt.subplot(1,2,2)
    plt.plot(x, w_def_PHAA*1000, color = 'r', label = 'PHAA')
    plt.plot(x, w_def_PLAA*1000, color = 'b', label = 'PLAA')
    plt.plot(x, w_def_NHAA*1000, color = 'k', label = 'NHAA')
    plt.plot(x, w_def_NLAA*1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$w$ (x) [mm]')
    plt.xlabel('x[m]')

    # =============================================== EIXO NEUTRO =========================================================
    alpha_eixo_PHAA = eixo_neutro(Mz_PHAA[0], My_PHAA[0], Izz_s, Iyy_s)
    alpha_eixo_PLAA = eixo_neutro(Mz_PLAA[0], My_PLAA[0], Izz_s, Iyy_s)
    alpha_eixo_NHAA = eixo_neutro(Mz_NHAA[0], My_NHAA[0], Izz_s, Iyy_s)
    alpha_eixo_NLAA = eixo_neutro(Mz_NLAA[0], My_NLAA[0], Izz_s, Iyy_s)

    print('======================== EIXO NEUTRO ========================')
    print('alpha_eixo_PHAA = ', alpha_eixo_PHAA*180/np.pi, 'deg')
    print('alpha_eixo_PLAA = ', alpha_eixo_PLAA*180/np.pi, 'deg')
    print('alpha_eixo_NHAA = ', alpha_eixo_NHAA*180/np.pi, 'deg')
    print('alpha_eixo_NLAA = ', alpha_eixo_NLAA*180/np.pi, 'deg')

    # ================================================ CAMPO DE TENSOES EM X = 0 (RAIZ) ===================================
    sigma_x(N_cel, h1, t1, t2, t3, e, b, E1, E2, E3, E1, My_PHAA[0], Mz_PHAA[0], Iyy_s, Izz_s, point= 'PHAA', centroide=(ys, zs), alpha=alpha_eixo_PHAA)
    sigma_x(N_cel, h1, t1, t2, t3, e, b, E1, E2, E3, E1, My_PLAA[0], Mz_PLAA[0], Iyy_s, Izz_s, point= 'PLAA', centroide=(ys, zs), alpha=alpha_eixo_PLAA)
    sigma_x(N_cel, h1, t1, t2, t3, e, b, E1, E2, E3, E1, My_NHAA[0], Mz_NHAA[0], Iyy_s, Izz_s, point= 'NHAA', centroide=(ys, zs), alpha=alpha_eixo_NHAA)
    sigma_x(N_cel, h1, t1, t2, t3, e, b, E1, E2, E3, E1, My_NLAA[0], Mz_NLAA[0], Iyy_s, Izz_s, point= 'NLAA', centroide=(ys, zs), alpha=alpha_eixo_NLAA)



    # =============================================== MAXIMAS DEFLEXOES ==================================================

    print('======================== MAXIMAS DEFLEXOES ========================')
    print('Número de células = ', N_cel)
    print('v_max_PHAA = ', np.max(np.abs(v_def_PHAA))*1000, 'mm')
    print('v_max_PLAA = ', np.max(np.abs(v_def_PLAA))*1000, 'mm')
    print('v_max_NHAA = ', np.max(np.abs(v_def_NHAA))*1000, 'mm')
    print('v_max_NLAA = ', np.max(np.abs(v_def_NLAA))*1000, 'mm')

    print('w_max_PHAA = ', np.max(np.abs(w_def_PHAA))*1000, 'mm')
    print('w_max_PLAA = ', np.max(np.abs(w_def_PLAA))*1000, 'mm')
    print('w_max_NHAA = ', np.max(np.abs(w_def_NHAA))*1000, 'mm')
    print('w_max_NLAA = ', np.max(np.abs(w_def_NLAA))*1000, 'mm')
    plt.show()