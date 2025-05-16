import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import pandas as pd 

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

def plot_esforco_cortante(x, Vz_PHAA, Vz_PLAA, Vz_NHAA, Vz_NLAA, Vy_PHAA, Vy_PLAA, Vy_NHAA, Vy_NLAA):
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

def plot_momento_fletor(x, Mz_PHAA,Mz_PLAA, Mz_NHAA, Mz_NLAA, My_PHAA, My_PLAA, My_NHAA, My_NLAA):
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


def plot_deflexao(x, v_def_PHAA, v_def_PLAA, v_def_NHAA, v_def_NLAA, w_def_PHAA, w_def_PLAA, w_def_NHAA, w_def_NLAA):
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


def plot_torcao(x, mx_PHAA, mx_PLAA, mx_NHAA, mx_NLAA, Tx_PHAA, Tx_PLAA, Tx_NHAA, Tx_NLAA):
    plt.figure(dpi = 100, figsize = (13,4))
    plt.subplot(1,2,1)
    plt.plot(x, mx_PHAA/1000, color = 'r', label = 'PHAA')
    plt.plot(x, mx_PLAA/1000, color = 'b', label = 'PLAA')
    plt.plot(x, mx_NHAA/1000, color = 'k', label = 'NHAA')
    plt.plot(x, mx_NLAA/1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$m_x$ (x) [kN m/ m]')
    plt.xlabel('x[m]')

    plt.subplot(1,2,2)
    plt.plot(x, Tx_PHAA/1000, color = 'r', label = 'PHAA')
    plt.plot(x, Tx_PLAA/1000, color = 'b', label = 'PLAA')
    plt.plot(x, Tx_NHAA/1000, color = 'k', label = 'NHAA')
    plt.plot(x, Tx_NLAA/1000, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$T$ (x) [kNm]')
    plt.xlabel('x[m]')

def plot_angulo_torcao(x, theta_x_PHAA, theta_x_PLAA, theta_x_NHAA, theta_x_NLAA):
    plt.figure(dpi = 100, figsize = (13,4))
    plt.plot(x, theta_x_PHAA*180/np.pi, color = 'r', label = 'PHAA')
    plt.plot(x, theta_x_PLAA*180/np.pi, color = 'b', label = 'PLAA')
    plt.plot(x, theta_x_NHAA*180/np.pi, color = 'k', label = 'NHAA')
    plt.plot(x, theta_x_NLAA*180/np.pi, color = 'g', label = 'NLAA')
    plt.legend()
    plt.grid()
    plt.ylabel(r'$\theta$ (x) [º]')
    plt.xlabel('x[m]')