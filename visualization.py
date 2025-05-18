import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

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
    plt.show()

def plot_subplot(x, u, v, w, z, p, q, r, l,
                 f = 1000, 
                 title = "", 
                 ylabel = ["", ""], 
                 xlabel = ["", ""], 
                 label = ["", "", "", ""]):
    plt.figure(dpi = 100, figsize = (13,4))
    plt.subplot(1,2,1)
    plt.plot(x, u/f, color = 'r', label = label[0])
    plt.plot(x, v/f, color = 'b', label = label[1])
    plt.plot(x, w/f, color = 'k', label = label[2])
    plt.plot(x, z/f, color = 'g', label = label[3])
    plt.legend()
    plt.grid()
    plt.ylabel(ylabel[0])
    plt.xlabel(xlabel[0])

    plt.subplot(1,2,2)
    plt.plot(x, p/f, color = 'r', label = label[0])
    plt.plot(x, q/f, color = 'b', label = label[1])
    plt.plot(x, r/f, color = 'k', label = label[2])
    plt.plot(x, l/f, color = 'g', label = label[3])
    plt.legend()
    plt.grid()
    plt.ylabel(ylabel[1])
    plt.xlabel(xlabel[1])

    plt.suptitle(title)
    # plt.show()

def plot_plot(x, u, v, w, z, 
                f = 1000,
                title = "", 
                ylabel = "", 
                xlabel = "", 
                label = ["", "", "", ""]):
    plt.figure(dpi = 100, figsize = (13,4))
    plt.plot(x, u/f, color = 'r', label = label[0])
    plt.plot(x, v/f, color = 'b', label = label[1])
    plt.plot(x, w/f, color = 'k', label = label[2])
    plt.plot(x, z/f, color = 'g', label = label[3])
    plt.legend()
    plt.grid()
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

def plot_sigma(sigma, l, htot, f = 1e6,
               i_ref_plot = 0, 
               centroide = [], alpha = None, 
               title = ""):
    # Determinando a malha
    resolution = sigma.shape
    y = np.linspace(0, l, resolution[0])
    z = np.linspace(0, htot, resolution[0])
    Y, Z = np.meshgrid(y, z)

    # --- Plot ---
    plt.figure(figsize=(10, 4))
    extent = [0, l, 0, htot]
    im = plt.imshow(sigma[:, :, i_ref_plot]/f, origin='lower', extent=extent,
                    aspect='auto', cmap='plasma')
    

    # Encontra o valor máximo ignorando NaNs
    max_valor = np.nanmax(np.abs(sigma[:, :, i_ref_plot]))

    # Encontra os índices do valor máximo (ignorando NaNs)
    indices = np.where(np.abs(sigma[:, :, i_ref_plot]) == max_valor)

    # Extrai a primeira ocorrência (se houver múltiplos máximos)
    linha, coluna = indices[0][0], indices[1][0]

    # Eixo neutro 
    x_c, y_c = centroide
    if alpha is not None:
            # Define uma linha grande o suficiente para cruzar a imagem
            comprimento = max(l, htot) * 2

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
    plt.plot(Y[linha, coluna], Z[linha, coluna], 'ro', label = 'Máxima tensão = '+ str(round(sigma[linha, coluna, i_ref_plot]/f, 2)) + 'MPa')  
    plt.colorbar(im, label=r'$\sigma_x(y,z)$ [MPa]')
    plt.xlabel('y [mm]')
    plt.ylabel('z [mm]')
    plt.title(title)
    plt.xlim([-0.01, l+0.01])
    plt.ylim([-0.01, htot+0.01])
    plt.tight_layout()
    plt.legend()

def plot_qfield(q, l, htot,
               i_ref_plot = 0,
               f = 1e6, 
               centroide = [], 
               title = ""):
    # Determinando a malha
    resolution = q.shape
    y = np.linspace(0, l, resolution[0])
    z = np.linspace(0, htot, resolution[0])
    Y, Z = np.meshgrid(y, z)

    # --- Plot ---
    plt.figure(figsize=(10, 4))
    extent = [0, l, 0, htot]
    im = plt.imshow(q/f, origin='lower', extent=extent,
                    aspect='auto', cmap='plasma')
    

    # Eixo neutro 
    x_c, y_c = centroide
    plt.text(x_c + 0.2*x_c, y_c, r"$q_1$ = " + str(round(q[0, 0], 3)))
    plt.text(x_c + 0.2*x_c, y_c - 0.2*y_c, r"$q_2$ = " + str(round(q[1, 0], 3)))
    plt.plot(x_c, y_c, 'ko', label = 'Centroide ponderado') 
    # plt.plot(Y[linha, coluna], Z[linha, coluna], 'ro', label = 'Máxima tensão = '+ str(round(sigma[linha, coluna, i_ref_plot]/f, 2)) + 'MPa')  
    plt.colorbar(im, label=r'$q(y,z)$ [kN/m]')
    plt.xlabel('y [mm]')
    plt.ylabel('z [mm]')
    plt.title(title)
    plt.xlim([-0.01, l+0.01])
    plt.ylim([-0.01, htot+0.01])
    plt.tight_layout()
    plt.legend()


