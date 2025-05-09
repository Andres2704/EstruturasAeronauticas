import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import pandas as pd 

def ve_function(v, rho, rho0 = 1.225):
    return v*np.sqrt(rho/rho0)

def dyn_pressure(rho, v):
    return 0.5*rho*v*v

def load_factor(q, Cl, S, W, f = 1.25):
    return f*q*Cl*S/W

def plot_vn(rho, W, S, Vd = 0, Vc = 0, Cla = 0, c_bar = 0, f = 1.25, rho0 = 1.225, Clmax = 0, Clmin = 0, n_max = 0, n_min = 0, plot = True):

   # ===========================================================================
    # ============================= DIAGRAMA V-n manobra ========================
    # ===========================================================================

    # Determinando a velocidade de mergulho ou cruzeiro, com base na entrada
    if (Vc>0):
        Vd = 1.296*Vc
        print('Velocidade de cruzeiro fornecida, neste caso Vd = 1.296Vc = ', np.round(Vd, 3), 'm/s') #modifiquei para 1.296 pq 259.3 = 200*1.296 (igual no manual do avião)
    elif (Vd > 0):
        Vc = 0.771*Vd
        print('Velocidade de mergulho fornecida, neste caso Vc = 0.771Vd = ', np.round(Vc, 3), 'm/s') #msm coisa, 200 = 259.3*0,771

    if (Vc == 0 and Vd == 0):
        print('A velocidade de cruzeiro e mergulho são nulas, favor rever')

    Ved = ve_function(Vd, rho)
    Vec = ve_function(Vc, rho)

    # Calculo de Vs_p e Vs_m para determinar o limite inferior de n = f(v^2)
    Vs_p = np.sqrt(2*W/(f*Clmax*S*rho0))
    Vs_m = np.sqrt(2*W/(abs(Clmin*S*rho0)))

    # Determinando o n_max e n_min caso não sejam fornecidos
    if (n_max == 0 or n_min == 0):
        n_max = 2.1 + 24000/(W + 10000) # FAR-23 - Aeronave de categoria normal
        n_min = -0.4*n_max
        if n_min > -1: n_min = -1

    n_ult_pos = 1.5*n_max
    n_ult_neg = 1.5*n_min
    print('n_max = ', np.round(n_max, 3))
    print('n_ult_+ = ', np.round(n_ult_pos, 3))
    print('n_min = ', np.round(n_min, 3))
    print('n_ult_- = ', np.round(n_ult_neg, 3))

    # Determinando a velocidade em Clmax Nmax (PHAA)
    V_PHAA = np.sqrt(n_max)*Vs_p
    V_NHAA = np.sqrt(abs(n_min))*Vs_m

    print('V_PHAA = ', np.round(V_PHAA, 3), 'm/s')
    print('V_NHAA = ', np.round(V_NHAA, 3), 'm/s')

    # Calculando a curva de estol positivo e negativo
    n_points = 1000
    V_estolpos = np.linspace(0, V_PHAA, n_points)
    V_estolneg = np.linspace(0, V_NHAA, n_points)
    npos = load_factor(dyn_pressure(rho0, V_estolpos), Clmax, S, W)
    nneg = load_factor(dyn_pressure(rho0, V_estolneg), Clmin, S, W, f = 1.0)

    # Plot do diagrama de manobra
    lw = 2
    plt.figure()

    # Plot manobra positiva
    plt.plot(V_estolpos, npos, label = 'Manobra Positiva', color = 'b', lw = lw)
    plt.hlines(n_max, V_PHAA, Vd, color = 'b', lw = lw)
    plt.vlines(Vs_p, 0, 1, color = 'b', lw = lw, ls = ':') # Vs+
    plt.vlines(Ved, 0, n_max, color = 'b', lw = lw)# Vd

    # Plot manobra negativa
    plt.plot(V_estolneg, nneg, label = 'Manobra Negativa', color = 'k', lw = lw)
    plt.hlines(n_min, V_NHAA, Vd, color = 'k', lw = lw)
    plt.hlines(0, 0, Ved, color = 'k', lw = lw, ls = ':')
    plt.vlines(Ved, n_min, 0, color = 'k', lw = lw)
    plt.vlines(Ved, n_min, n_max, color = 'k', ls = ':')
    plt.text(1.005*Ved, 0.5*n_max, r'$V_{d}$', fontsize=10)

    plt.vlines(Vs_m, 0, -1, color = 'k', lw = lw, ls = ':')

    # Anotações

#     plt.text(0.6*Ved, 1.05*n_max, r'$n_{max}$ = ' + str(np.round(n_max, 3)), fontsize=12)
#     plt.text(0.7*Ved, 1.2*n_min, r'$n_{min}$ = ' + str(np.round(n_min, 3)), fontsize=12)

    plt.text(0.95*Vs_p, -0.25, r'$V_{s+}$', fontsize=10)
    plt.text(0.99*Vs_m, 0.5*n_min, r'$V_{s-}$', fontsize=10)

    # Linha de cruzeiro
    idx_closest_vc = np.argmin(np.abs(V_estolneg - Vec))
    nneg_vc = nneg[idx_closest_vc]

    plt.vlines(Vec, nneg_vc, n_max, color = 'k', ls = ':')
    plt.text(1.01*Vec, 0.5*n_max, r'$V_{c}$', fontsize=10)

    # ===========================================================================
    # ============================= DIAGRAMA V-n rajadas ========================
    # ===========================================================================
    if (Cla == 0 or c_bar == 0):
        print('Valor dCl/dAlpha ou c_barra não fornecidos, curva de rajada não será mostrada')
    else:
        V_raj_c = np.linspace(0, Vec, n_points)
        V_raj_d = np.linspace(0, Ved, n_points)

        # Voo subsônico
        mu = 2*(W/S)/(rho*9.81*c_bar*Cla)
        K = 0.88*mu/(5.3 + mu)

        # Rajada na velocidade de cruzeiro Vc
        Ude = 9
        U = K*Ude
        Ue = ve_function(U, rho)

        # Curvas de rajada
        n_rajpos_c = 1 + 0.5*rho0*V_raj_c*S*Cla*Ue/W
        n_rajneg_c = 1 - 0.5*rho0*V_raj_c*S*Cla*Ue/W

        plt.plot(V_raj_c, n_rajpos_c, label = 'Rajada Vc', color = 'g', ls = '--')
        plt.plot(V_raj_c, n_rajneg_c, color = 'g', ls = '--')

        # Rajada na velocidade de mergulho Vd
        Ude = 4.5
        U = K*Ude
        Ue = ve_function(U, rho)

        # Curvas de rajada
        n_rajpos_d = 1 + 0.5*rho0*V_raj_d*S*Cla*Ue/W
        n_rajneg_d = 1 - 0.5*rho0*V_raj_d*S*Cla*Ue/W

        plt.plot(V_raj_d, n_rajpos_d, label = 'Rajada Vd', color = 'm', ls = '--')
        plt.plot(V_raj_d, n_rajneg_d, color = 'm', ls = '--')

        # Linha Vc a Vd para cada manobra
        Vc_Vd = np.linspace(Vec, Ved, int(np.ceil(n_points/2)))
        nc_nd_pos = np.linspace(np.max(n_rajpos_c), np.max(n_rajpos_d), int(np.ceil(n_points/2)))
        nc_nd_neg = np.linspace(np.min(n_rajneg_c), np.min(n_rajneg_d), int(np.ceil(n_points/2)))
        plt.plot(Vc_Vd, nc_nd_pos, color = 'g', ls = '--')
        plt.plot(Vc_Vd, nc_nd_neg, color = 'g', ls = '--')

    # Config plot
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('V [m/s]')
    plt.ylabel('n')
    plt.xlim([0, Vd*1.1])
#     plt.ylim([1.1*n_min, 1.1*n_max])

    plt.plot(V_PHAA, n_max, 'ro')  # PHAA
    plt.text(V_PHAA*1.01, n_max*1.01, 'PHAA', color='r')

    plt.plot(Vd, n_max, 'co')  # PLAA
    plt.text(Vd*1.01, n_max*1.01, 'PLAA', color='c')

    plt.plot(V_NHAA, n_min, 'ro')  # NHAA
    plt.text(V_NHAA*1.01, n_min*1.01, 'NHAA', color='r')

    plt.plot(Vd, n_min, 'co')  # NLAA
    plt.text(Vd*1.01, n_min*1.01, 'NLAA', color='c')

    if plot: plt.show()

    return V_NHAA, V_PHAA

def lift_Schrenk(L, b, lamb, sign = 1, n = 100):
    n_points = n

    if sign == 1:
        x = np.linspace(0, b/2, n_points)
    else:
        x = np.linspace(-b/2, 0, n_points)

    w_elip = 4*L*np.sqrt(1 - (2*sign*x/b)**2)/(np.pi * b)
    w_trap = 2*L*(1 - 2*sign*x*(1 - lamb)/b)/(b*(1 + lamb))

    w_eff = 0.5*(w_elip + w_trap)

    return x, w_eff

def drag_linear(rho, V, S, Cd, b, n = 100):
    D = 0.5*rho*V**2 * S * Cd
    w_drag = np.ones(n)*D*4/b
    return w_drag

def drag_total(rho, V, S, Cd):
    D = 0.5*rho*V**2 * S * Cd
    return D

def shear_force(x, w):
    V = np.zeros(len(x))
    for i in range(len(x)-1, 0, -1):
        V[i-1] = V[i] + w[i]*(x[i] - x[i-1])
    return V

def bending_moment(x, V):
    M = np.zeros(len(x))
    for i in range(len(x)-1, 0, -1):
        M[i-1] = M[i] + V[i]*(x[i] - x[i-1])
    return M

def integrate_deflection(x, beta, cond = [0, 0]):
    # x é a coordenada espacial E [0, b/2]
    # beta é o momento fletor/(E0 * I)

    # Deflexão
    v = np.zeros(len(x))
    v[0] = cond[0]

    # d (Deflexao)/dx
    u = np.zeros(len(x))
    u[0] = cond[1]

    for i in range(1, len(x)):
        u[i] = u[i-1] + beta[i]*(x[i] - x[i-1])
        v[i] = v[i-1] + u[i]*(x[i] - x[i-1])

    return v

def transf_aero2struct(alpha, F):
    f_trans = np.zeros((2, np.shape(F)[1]))
    f_trans[0, :] = F[0,:]*np.cos(alpha) - F[1,:]*np.sin(alpha)
    f_trans[1, :] = F[0,:]*np.sin(alpha) + F[1,:]*np.cos(alpha)
    return f_trans

def plot_3d_forces(x, F_distrib, L_total, D_total, title="Distribuição de Forças Aerodinâmicas"):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Extrai distribuições brutas
    w_drag = F_distrib[0,:]
    w_lift = F_distrib[1,:]

    # Área sob a curva (para normalização)
    drag_integral = np.trapz(w_drag, x)
    lift_integral = np.trapz(w_lift, x)

    # Evita divisão por zero
    drag_integral = drag_integral if drag_integral != 0 else 1
    lift_integral = lift_integral if lift_integral != 0 else 1

    # Escalona distribuições para total desejado
    w_drag_scaled = w_drag * (D_total / drag_integral)
    w_lift_scaled = w_lift * (L_total / lift_integral)

    # Escalas para visualização
    max_force = max(np.max(np.abs(w_lift_scaled)), np.max(np.abs(w_drag_scaled)))
    scale_global = 0.3 * max(x) / max_force if max_force > 0 else 1

    # Linha da asa (eixo x)
    ax.plot(x, np.zeros_like(x), np.zeros_like(x), 'k-', linewidth=2, label='Asa')

    for i in range(len(x)):
        # Vetor de sustentação
        ax.quiver(x[i], 0, 0,
                  0, 0, w_lift_scaled[i] * scale_global,
                  color='blue', alpha=0.7,
                  label='Sustentação' if i == 0 else "",
                  arrow_length_ratio=0.05, linewidth=1.5)

        # Vetor de arrasto
        ax.quiver(x[i], 0, 0,
                  0, -w_drag_scaled[i] * scale_global, 0,
                  color='red', alpha=0.7,
                  label='Arrasto' if i == 0 else "",
                  arrow_length_ratio=0.05, linewidth=1.5)

    # Eixos e título
    ax.set_xlabel('Envergadura [m]')
    ax.set_ylabel('Profundidade [m] (Arrasto)')
    ax.set_zlabel('Altura [m] (Sustentação)')
    # ax.set_title(f"{title}\nL = {L_total:.1f} N | D = {D_total:.1f} N | L/D = {abs(L_total/D_total):.2f}")

    # Limites com padding
    padding = 0.1 * max(x)
    ax.set_xlim([min(x)-padding, max(x)+padding])

    max_drag = np.max(np.abs(w_drag_scaled)) * scale_global
    ax.set_ylim([-max_drag-padding, padding])

    max_lift = np.max(w_lift_scaled) * scale_global
    min_lift = np.min(w_lift_scaled) * scale_global
    ax.set_zlim([min(0, min_lift-padding), max_lift+padding])

    # Ângulo de visão
    ax.view_init(elev=25, azim=-50)

    # Legendas sem duplicação
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')

    plt.tight_layout()
    plt.show()

def det_alpha(cl0, cla, cl):
    return (cl - cl0)/cla

def calc_area_prop(N, E, t, h1, b, e):
    # Propriedades de cada material
    E1, E2, E3 = E
    E0 = E1                     # Definindo E0 = E1
    # Espessura de cada material
    t1, t2, t3 = t

    # Comprimento e altura total da caixa de asa
    l = (N + 2)*e + (N+1)*b
    h = h1 + 2*t2 + 2*t3

    # Área ponderada
    As = (E1/E0)*(N+1)*t1*h1 + 2*(E2/E0)*(N+1)*b*t2 + 2*(E3/E0)*l*t3

    # Centroide ponderado z*
    zs = ((E1/E0)*(N+1)*h1*h*t1/2 + (E2/E0)*(N+1)*b*h*t2 + (E3/E0)*t3*l*h)/As

    # Centroide ponderado y*
    sum_ys = 0
    sum_yss = 0
    for i in range(1, N+2):
        sum_ys = sum_ys + i*e + (2*i - 1)*b/2
        sum_yss = sum_yss + (i*e + (2*i - 1)*b/2)**2

    ys = (((E1/E0)*h1*t1 + 2*(E2/E0)*b*t2)*sum_ys + (E3/E0)*l**2 * t3)/As

    # Momento de inércia ponderado Izz*
    Iz0z0s = (E1/E0)*(N+1)*(h1*t3**3/12 + h**2*t1*h1/4) + \
             (E2/E0)*(N+1)*(t2*b**3/6 + b*t2*( (t3+0.5*t2)**2 + (h-(t3+0.5*t2))**2 )) + \
             (E3/E0)*(1+0)*(t3*l**3/6 + l*t3*( (0.5*t3**2)    + (h - 0.5*t3**2)))
    Izzs = Iz0z0s - ys**2 * As

    # Momento de inércia ponderado Iyy*
    Iy0y0s = (E1/E0)*((N+1)*t1*h1**3/12 + h1*t1*sum_yss) + \
             (E2/E0)*((N+1)*b*t2**3/12 + b*t2*sum_yss) + \
           2*(E3/E0)*(l*t3**2/12 + l**3 * t3 / 4)
    Iyys = Iy0y0s - zs**2 * As

    return [As, zs, ys, Izzs, Iyys]

def eixo_neutro(Mz, My, Izzs, Iyys):    # Considerando que sempre I_yz = 0
    alpha_eixo = np.arctan2(Mz*Iyys, My*Izzs)
    return alpha_eixo

def desenhar_vigas_T(n, t1, t2, t3, b, h1, e, ponto_centroide=None, alpha=None):
    largura_total = n * b + (n + 1) * e
    altura_total = 2 * t3 + 2 * t2 + h1

    fig, ax = plt.subplots(figsize=(12, 6))

    # Revestimentos (superior e inferior)
    ax.add_patch(patches.Rectangle((0, altura_total - t3), largura_total, t3, color='gray'))
    ax.add_patch(patches.Rectangle((0, 0), largura_total, t3, color='gray'))

    for i in range(n):
        x_base = e + i * (b + e)

        # Aba inferior
        ax.add_patch(patches.Rectangle((x_base, t3), b, t2, color='dodgerblue'))

        # Alma
        x_alma = x_base + (b - t1) / 2
        y_alma = t3 + t2
        ax.add_patch(patches.Rectangle((x_alma, y_alma), t1, h1, color='red'))

        # Aba superior
        ax.add_patch(patches.Rectangle((x_base, y_alma + h1), b, t2, color='dodgerblue'))

    # Itens da legenda
    legend_items = [
        patches.Patch(color='gray', label='Revestimento'),
        patches.Patch(color='red', label='Alma'),
        patches.Patch(color='dodgerblue', label='Aba')
    ]

    # Centroide e eixo neutro
    if ponto_centroide is not None:
        x_c, y_c = ponto_centroide
        ax.plot(x_c, y_c, 'ko')  # ponto preto
        legend_items.append(Line2D([0], [0], marker='o', color='black', label='Centroide ponderado', linestyle=''))

        if alpha is not None:
            # Converte o ângulo de graus para radianos
            theta = np.radians(alpha)

            # Define uma linha grande o suficiente para cruzar a imagem
            comprimento = max(largura_total, altura_total) * 2

            # Direção da linha
            dx = np.cos(theta)
            dy = np.sin(theta)

            # Ponto inicial e final da linha
            x_start = x_c - comprimento * dx
            x_end   = x_c + comprimento * dx
            y_start = y_c - comprimento * dy
            y_end   = y_c + comprimento * dy

            # Desenha a linha pontilhada
            ax.plot([x_start, x_end], [y_start, y_end], linestyle='--', color='black', linewidth=1)
            legend_items.append(Line2D([0], [0], linestyle='--', color='black', label='Eixo neutro'))

    # Legenda
    ax.legend(handles=legend_items, loc='center left', bbox_to_anchor=(1.01, 0.5))

    # Ajustes visuais
    ax.set_aspect('equal')
    ax.set_xlim(-0.05 * largura_total, 1.05 * largura_total)
    ax.set_ylim(-0.05 * altura_total, 1.05 * altura_total)
    ax.axis('off')
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    plt.show()

def sigma_x(N_cel, h1, t1, t2, t3, e, b,
                  E1, E2, E3, E0,
                  My, Mz, Iyy_star, Izz_star,
                  resolution=300, point = '', centroide = (0, 0), alpha = 0):
    
    # Determinando a altura total e comprimento total da caixa
    h = 2*t3 + 2*t2 + h1
    l = (N_cel+1)*b + (N_cel+2)*e 

    # Determinando a malha
    y = np.linspace(0, l, resolution)
    z = np.linspace(0, h1 + 2*t2 + 2*t3, resolution)
    Y, Z = np.meshgrid(y, z)

    # Região de existencia do revestimento
    i_rev_inf = (Z>=0) & (Z<=t3)
    i_rev_sup = (Z>= t3+2*t2+h1) & (Z<=2*t3 + 2*t2 + h1)

    # Campo do módulo E
    E_field = np.nan*np.zeros((resolution, resolution))

    # Campo E no revestimento
    E_field[i_rev_inf] = E3 
    E_field[i_rev_sup] = E3

    # Determinando a existencia do material na alma e na aba
    for k in range(resolution):
        for j in range(resolution):
            for i in range(N_cel+1):
                # Almas 
                if ((Y[k,j]>= ((i+1)*e + i*b + (b-t1)/2)) & (Y[k,j]<= (((i+1)*e + i*b + (b+t1)/2))))&((Z[k,j]>=(t3+t2)) & (Z[k,j]<=(t3+t2+h1))):
                    E_field[k,j] = E1 
                
                # Aba superior 
                if ((Y[k,j]>=((i+1)*e + i*b)) & (Y[k,j]<=(i+1)*e + (i+1)*b)) & ((Z[k,j]>=(h-t3-t2)) & (Z[k,j]<= h-t3)):
                    E_field[k,j] = E2
                
                # Aba inferior 
                if ((Y[k,j]>=((i+1)*e + i*b)) & (Y[k,j]<=(i+1)*e + (i+1)*b)) & ((Z[k,j]>=(t3)) & (Z[k,j]<= t3+t2)):
                    E_field[k,j] = E2
                    
    # Campo de tensões
    sigma = E_field * (Z*My/Iyy_star - Y*Mz/Izz_star) / E0 

    # --- Plot ---
    plt.figure(figsize=(10, 4))
    extent = [0, l, 0, 2*t3+2*t2 + h1]
    im = plt.imshow(sigma/1e6, origin='lower', extent=extent,
                    aspect='auto', cmap='plasma')
    

    # Encontra o valor máximo ignorando NaNs
    max_valor = np.nanmax(np.abs(sigma))

    # Encontra os índices do valor máximo (ignorando NaNs)
    indices = np.where(np.abs(sigma) == max_valor)

    # Extrai a primeira ocorrência (se houver múltiplos máximos)
    linha, coluna = indices[0][0], indices[1][0]

    print('Máxima tensão do ponto' + point + ' = ', round(sigma[linha, coluna]/1e6, 3), '[MPa]')

    # Eixo neutro 
    x_c, y_c = centroide
    if alpha is not None:
            # Define uma linha grande o suficiente para cruzar a imagem
            comprimento = max(l, h) * 2

            # Direção da linha
            dx = np.cos(alpha)
            dy = np.sin(alpha)

            # Ponto inicial e final da linha
            x_start = x_c - comprimento * dx
            x_end   = x_c + comprimento * dx
            y_start = y_c - comprimento * dy
            y_end   = y_c + comprimento * dy

            # Desenha a linha pontilhada
            plt.plot([x_start, x_end], [y_start, y_end], linestyle='--', color='black', linewidth=1, label = 'Eixo Neutro')
            

    plt.plot(x_c, y_c, 'ko', label = 'Centroide ponderado') 
    plt.plot(Y[linha, coluna], Z[linha, coluna], 'ro', label = 'Máxima tensão = '+ str(round(sigma[linha, coluna]/1e6, 2)) + 'Mpa')  
    plt.colorbar(im, label=r'$\sigma_x(y,z)$ [MPa]')
    plt.xlabel('y [mm]')
    plt.ylabel('z [mm]')
    plt.title('Campo de Tensões σₓ(y,z) - '+ point)
    plt.xlim([-0.01, l+0.01])
    plt.ylim([-0.01, h+0.01])
    plt.tight_layout()
    plt.legend()
    