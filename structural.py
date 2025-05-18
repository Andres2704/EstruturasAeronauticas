from func_solver import * 
from fun_plot import * 

class structure():
    def __init__(self, **kwargs):
        self.init_aircraft(**kwargs)
        self.init_section(**kwargs)

    def init_aircraft(self, **kwargs):
        self.Vd = kwargs.get('Vd')
        self.Vc = kwargs.get('Vc')
        self.rho = kwargs.get('rho')
        self.rho0 = kwargs.get('rho0')
        self.factor = kwargs.get('rho0')
        self.Cl_alpha = kwargs.get('CL_alpha')
        self.CMA = kwargs.get('CMA')
        if self.factor == 0:
            self.factor = 1.25

        self.CL_max = kwargs.get('CL_max')
        self.CL_min = kwargs.get('CL_min')
        self.W      = kwargs.get('W')
        self.S      = kwargs.get('S')
        self.n_max  = kwargs.get('n_max')
        self.n_min  = kwargs.get('n_min')
        self.N_points = kwargs.get('npoints')
        self.lamb   = kwargs.get('lamb')
        self.env    = kwargs.get('env')
        self.x_coord = np.linspace(0, self.env/2, self.N_points)

    def init_section(self, **kwargs):
        self.type = kwargs.get('type')

        if self.type.lower() == 'points':
            pass 
        else: 
            self.N_cel = kwargs.get('N_cel')
            self.E1, self.E2, self.E3 = kwargs.get('E')
            self.G1, self.G2, self.G3 = kwargs.get('G')
            self.t1, self.t2, self.t3 = kwargs.get('t')
            self.h1 = kwargs.get('h1')
            self.l = kwargs.get('l')
            self.b_por = kwargs.get('b_por')
            self.b = (self.b_por/100)*self.l/(self.N_cel + 1)
            self.e = (self.l - (self.N_cel + 1)*self.b)/(self.N_cel + 2)
            self.htot = self.h1 + 2*self.t2 + 2*self.t3
            self.i_ref_mat = kwargs.get('i_ref_mat')
            self.E0 = kwargs.get('E')[self.i_ref_mat]
            self.G0 = kwargs.get('G')[self.i_ref_mat]

class aero_struct(structure):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)  # Chama o __init__ da classe pai
        self.vn_load()
        self.area_properties()

    def vn_load(self, **kwargs):
        # ===========================================================================
        # ============================= DIAGRAMA V-n manobra ========================
        # ===========================================================================
        if (self.Vc == 0 and self.Vd == 0):
            print('A velocidade de cruzeiro e mergulho são nulas, favor rever')
            return 0 
        
        self.Ved = self.Vd * np.sqrt(self.rho/self.rho0)
        self.Vec = self.Vc * np.sqrt(self.rho/self.rho0)

        # Calculo de Vs_p e Vs_m para determinar o limite inferior de n = f(v^2)
        self.Vs_p = np.sqrt(2*self.W/(self.factor*self.CL_max*self.S*self.rho0))
        self.Vs_m = np.sqrt(2*self.W/(abs(self.CL_min*self.S*self.rho0)))

        # Determinando o n_max e n_min caso não sejam fornecidos
        if (self.n_max == 0 or self.n_min == 0):
            self.n_max = 2.1 + 24000/(self.W + 10000) # FAR-23 - Aeronave de categoria normal
            self.n_min = -0.4*self.n_max
            if self.n_min > -1: self.n_min = -1

        self.L_max = self.n_max*self.W 
        self.L_min = self.n_min*self.W 
        self.n_ult_pos = 1.5*self.n_max
        self.n_ult_neg = 1.5*self.n_min
        
        # Determinando a velocidade em Clmax Nmax (PHAA)
        self.V_PHAA = np.sqrt(self.n_max)*self.Vs_p
        self.V_NHAA = np.sqrt(abs(self.n_min))*self.Vs_m
       

        # ===========================================================================
        # ============================= DIAGRAMA V-n rajadas ========================
        # ===========================================================================
        if (self.Cl_alpha == 0 or self.CMA == 0):
            print('Valor dCl/dAlpha ou c_barra não fornecidos, curva de rajada não será mostrada')
        else:

            # Voo subsônico
            mu = 2*(self.W/self.S)/(self.rho*9.81*self.CMA*self.Cl_alpha)
            K = 0.88*mu/(5.3 + mu)

            # Rajada na velocidade de cruzeiro Vc
            Ude = 9
            U = K*Ude
            self.Ue_cruise = U*np.sqrt(self.rho/self.rho0)

            # Rajada na velocidade de mergulho Vd
            Ude = 4.5
            U = K*Ude
            self.Ue_vd = U*np.sqrt(self.rho/self.rho0)
          
    def lift_distribution(self, **kwargs):
        L = kwargs.get('L')

        w_lift_elip = 4*L*np.sqrt(1 - (2*self.x_coord/self.env)**2)/(np.pi * self.env)
        w_lift_trap = 2*L*(1 - 2*self.x_coord*(1 - self.lamb)/self.env)/(self.env*(1 + self.lamb))
        w_lift_eff = 0.5*(w_lift_elip + w_lift_trap)
        return w_lift_trap, w_lift_elip, w_lift_eff

    def drag_distribution(self, **kwargs):
        Cd = kwargs.get('Cd')
        V = kwargs.get('V')
        D = 0.5*self.rho*V**2*self.S*Cd
        w_drag = np.ones(self.N_points)* D * 4 / self.env
        return D, w_drag 
  
    def total_loads_struc(self, **kwargs):
        alpha = kwargs.get('AOA')
        F = kwargs.get('loads_aero')
        f_trans = np.zeros((2, self.N_points))
        f_trans[0, :] = F[0,:]*np.cos(alpha) - F[1,:]*np.sin(alpha)
        f_trans[1, :] = F[0,:]*np.sin(alpha) + F[1,:]*np.cos(alpha)
        return f_trans 

    def x_neg_integral(self, **kwargs):
        H = np.zeros(self.N_points)
        V = kwargs.get('V')
        for i in range(self.N_points-1, 0, -1):
            H[i-1] = H[i] + V[i]*(self.x_coord[i] - self.x_coord[i-1])

        return H
    
    def x_pos_2nd_integral(self, **kwargs):
        # Variable
        v = np.zeros(self.N_points)
        v[0] = kwargs.get('cond')

        # Variable derivate 
        u = np.zeros(self.N_points)
        u[0] = kwargs.get('cond_derivate')

        # Function to integrate 
        beta = kwargs.get('beta')
        for i in range(1, self.N_points):
            u[i] = u[i-1] + beta[i]*(self.x_coord[i] - self.x_coord[i-1])
            v[i] = v[i-1] + u[i]*(self.x_coord[i] - self.x_coord[i-1])

        return v  

    def area_properties(self, **kwargs):
        # Área ponderada
        self.As = (self.E1/self.E0)*(self.N_cel+1)*self.t1*self.h1 + \
                 + 2*(self.E2/self.E0)*(self.N_cel+1)*self.b*self.t2 + \
                 + 2*(self.E3/self.E0)*self.l*self.t3
        
        # Centroide ponderado z*
        self.zs =  ((self.E1/self.E0)*(self.N_cel+1)*self.h1*self.htot*self.t1/2 + \
                 + (self.E2/self.E0)*(self.N_cel+1)*self.b*self.htot*self.t2 + \
                 + (self.E3/self.E0)*self.t3*self.l*self.htot)/self.As

        # Centroide ponderado y*
        sum_ys = 0
        sum_yss = 0
        for i in range(1, self.N_cel+2):
            sum_ys = sum_ys + i*self.e + (2*i - 1)*self.b/2
            sum_yss = sum_yss + (i*self.e + (2*i - 1)*self.b/2)**2

        self.ys = (((self.E1/self.E0)*self.h1*self.t1 + \
             + 2*(self.E2/self.E0)*self.b*self.t2)*sum_ys +\
             +   (self.E3/self.E0)*self.l**2 * self.t3)/self.As

        # Momento de inércia ponderado Izz*
        self.Iz0z0s = (self.E1/self.E0)*(self.N_cel+1)*(self.h1*self.t3**3/12 + self.htot**2*self.t1*self.h1/4) + \
                (self.E2/self.E0)*(self.N_cel+1)*(self.t2*self.b**3/6 + self.b*self.t2*( (self.t3+0.5*self.t2)**2 + (self.htot-(self.t3+0.5*self.t2))**2 )) + \
                (self.E3/self.E0)*(1+0)*(self.t3*self.l**3/6 + self.l*self.t3*( (0.5*self.t3**2)    + (self.htot - 0.5*self.t3**2)))
        self.Izzs = self.Iz0z0s - self.ys**2 * self.As

        # Momento de inércia ponderado Iyy*
        self.Iy0y0s = (self.E1/self.E0)*((self.N_cel+1)*self.t1*self.h1**3/12 + self.h1*self.t1*sum_yss) + \
                (self.E2/self.E0)*((self.N_cel+1)*self.b*self.t2**3/12 + self.b*self.t2*sum_yss) + \
            2*(self.E3/self.E0)*(self.l*self.t3**2/12 + self.l**3 * self.t3 / 4)
        self.Iyys = self.Iy0y0s - self.zs**2 * self.As 
        
    def propagate_aero_critic(self, **kwargs):
        '''
        Function for aerodynamic load at four critical points of v-n diagram
        '''
        self.CL_PHAA = kwargs.get('CL_PHAA')
        self.CL_PLAA = 2*self.L_max/(self.rho*self.Vd**2 * self.S) 
        self.CL_NHAA = kwargs.get('CL_NHAA')
        self.CL_NLAA = 2*self.L_min/(self.rho*self.Vd**2 * self.S)
        alpha_PHAA = kwargs.get('alpha_PHAA')
        alpha_PLAA = kwargs.get('alpha_PLAA')
        alpha_NLAA = kwargs.get('alpha_NLAA')
        alpha_NHAA = kwargs.get('alpha_NHAA')
        cd_PHAA = kwargs.get('cd_PHAA')
        cd_PLAA = kwargs.get('cd_PLAA')
        cd_NLAA = kwargs.get('cd_NLAA')
        cd_NHAA = kwargs.get('cd_NHAA')

        _, _, w_sust_pos = self.lift_distribution(L = self.L_max)
        _, _, w_sust_neg = self.lift_distribution(L = self.L_min)
        D_PHAA, w_drag_PHAA = self.drag_distribution(Cd = cd_PHAA, V = self.V_PHAA)
        D_PLAA, w_drag_PLAA = self.drag_distribution(Cd = cd_PLAA, V = self.Vd)
        D_NHAA, w_drag_NHAA = self.drag_distribution(Cd = cd_NHAA, V = self.V_NHAA)
        D_NLAA, w_drag_NLAA = self.drag_distribution(Cd = cd_NLAA, V = self.Vd)

        # Distributed loads in the aerodynamic reference frame
        F_distrib_PHAA = np.concatenate((w_drag_PHAA, w_sust_pos)).reshape(2, len(self.x_coord))
        F_distrib_PLAA = np.concatenate((w_drag_PLAA, w_sust_pos)).reshape(2, len(self.x_coord))
        F_distrib_NHAA = np.concatenate((w_drag_NHAA, w_sust_neg)).reshape(2, len(self.x_coord))
        F_distrib_NLAA = np.concatenate((w_drag_NLAA, w_sust_neg)).reshape(2, len(self.x_coord))

        # Distributed loads on the structural reference of frame 
        self.F_PHAA = self.total_loads_struc(AOA = alpha_PHAA, loads_aero = F_distrib_PHAA)
        self.F_PLAA = self.total_loads_struc(AOA = alpha_PLAA, loads_aero = F_distrib_PLAA)
        self.F_NHAA = self.total_loads_struc(AOA = alpha_NHAA, loads_aero = F_distrib_NHAA)
        self.F_NLAA = self.total_loads_struc(AOA = alpha_NLAA, loads_aero = F_distrib_NLAA)

        # Shear force y and z at each critical point 
        self.Vy_PLAA = self.x_neg_integral(cond = 0, V = -self.F_PLAA[0, :])
        self.Vy_PHAA = self.x_neg_integral(cond = 0, V = -self.F_PHAA[0, :])
        self.Vy_NLAA = self.x_neg_integral(cond = 0, V = -self.F_NLAA[0, :])
        self.Vy_NHAA = self.x_neg_integral(cond = 0, V = -self.F_NHAA[0, :])
        self.Vz_PLAA = self.x_neg_integral(cond = 0, V = -self.F_PLAA[1, :])
        self.Vz_PHAA = self.x_neg_integral(cond = 0, V = -self.F_PHAA[1, :])
        self.Vz_NLAA = self.x_neg_integral(cond = 0, V = -self.F_NLAA[1, :])
        self.Vz_NHAA = self.x_neg_integral(cond = 0, V = -self.F_NHAA[1, :])

        # Bending moment y and z at each critical point 
        self.My_PLAA = self.x_neg_integral(cond = 0, V = -self.Vy_PLAA)
        self.My_PHAA = self.x_neg_integral(cond = 0, V = -self.Vy_PHAA)
        self.My_NLAA = self.x_neg_integral(cond = 0, V = -self.Vy_NLAA)
        self.My_NHAA = self.x_neg_integral(cond = 0, V = -self.Vy_NHAA)
        self.Mz_PLAA = self.x_neg_integral(cond = 0, V = self.Vz_PLAA)
        self.Mz_PHAA = self.x_neg_integral(cond = 0, V = self.Vz_PHAA)
        self.Mz_NLAA = self.x_neg_integral(cond = 0, V = self.Vz_NLAA)
        self.Mz_NHAA = self.x_neg_integral(cond = 0, V = self.Vz_NHAA)

        # Deflections w and v at each critical point 
        self.v_def_NHAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = self.Mz_NHAA/(self.E0*self.Izzs)) 
        self.v_def_NLAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = self.Mz_NLAA/(self.E0*self.Izzs)) 
        self.v_def_PHAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = self.Mz_PHAA/(self.E0*self.Izzs)) 
        self.v_def_PLAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = self.Mz_PLAA/(self.E0*self.Izzs)) 
        self.w_def_NHAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = -self.My_NHAA/(self.E0*self.Iyys)) 
        self.w_def_NLAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = -self.My_NLAA/(self.E0*self.Iyys)) 
        self.w_def_PHAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = -self.My_PHAA/(self.E0*self.Iyys)) 
        self.w_def_PLAA = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = -self.My_PLAA/(self.E0*self.Iyys)) 

        # Neutral axis for each critical point 
        self.alpha_axis_PHAA = self.neutral_axis(Mz = self.Mz_PHAA, My = self.My_PHAA)
        self.alpha_axis_PLAA = self.neutral_axis(Mz = self.Mz_PLAA, My = self.My_PLAA)
        self.alpha_axis_NHAA = self.neutral_axis(Mz = self.Mz_NHAA, My = self.My_NHAA)
        self.alpha_axis_NLAA = self.neutral_axis(Mz = self.Mz_NLAA, My = self.My_NLAA)

        # Tension sigma_x for 
        self.sigma_max_PHAA, self.sigma_PHAA = self.sigma_x(My = self.My_PHAA, Mz = self.Mz_PHAA)
        self.sigma_max_PLAA, self.sigma_PLAA = self.sigma_x(My = self.My_PLAA, Mz = self.Mz_PLAA)
        self.sigma_max_NHAA, self.sigma_NHAA = self.sigma_x(My = self.My_NHAA, Mz = self.Mz_NHAA)
        self.sigma_max_NLAA, self.sigma_NLAA = self.sigma_x(My = self.My_NLAA, Mz = self.Mz_NLAA)

    def propagate_aero(self, **kwargs):
        '''
        Propagate the problem for a given point of analysis
        '''
        AOA     = kwargs.get('AOA')
        CD      = kwargs.get('CD')
        L       = kwargs.get('L')
        V       = kwargs.get('V')

        self.w_sust_trap, self.w_sust_elip, self.w_sust_point = self.lift_distribution(L = L)
        self.D_point, self.w_drag_point = self.drag_distribution(Cd = CD, V = V)
        
        self.F_distrib_point = np.concatenate((self.w_drag_point, self.w_sust_point)).reshape(2, len(self.x_coord))
        
        # Distributed loads on the structural reference of frame 
        self.F_point = self.total_loads_struc(AOA = AOA, loads_aero = self.F_distrib_point)
        
        # Shear force y and z at each critical point 
        self.Vy_point = self.x_neg_integral(cond = 0, V = -self.F_point[0, :])
        self.Vz_point = self.x_neg_integral(cond = 0, V = -self.F_point[1, :])

        # Bending moment y and z at each critical point 
        self.My_point = self.x_neg_integral(cond = 0, V = -self.Vy_point)
        self.Mz_point = self.x_neg_integral(cond = 0, V = self.Vz_point)

        # Deflections w and v at each critical point 
        self.v_def_point = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = self.Mz_point/(self.E0*self.Izzs)) 
        self.w_def_point = self.x_pos_2nd_integral(cond = 0, cond_derivate = 0, beta = -self.My_point/(self.E0*self.Iyys)) 

        # Neutral axis for each critical point 
        self.alpha_axis_point = self.neutral_axis(Mz = self.Mz_point, My = self.My_point)

        # Tension sigma x 
        self.sigma_max_point, self.sigma_point = self.sigma_x(My = self.My_point, Mz = self.Mz_point)

    def neutral_axis(self, Mz, My):
        return np.arctan2(Mz*self.Iyys, My*self.Izzs)
    
    def sigma_x(self, **kwargs):
        My = kwargs.get('My')
        Mz = kwargs.get('Mz')

        # Calculating the mesh
        y = np.linspace(0, self.l, self.N_points)
        z = np.linspace(0, self.htot, self.N_points)
        Y, Z = np.meshgrid(y, z)

        # Region of existence of material 3
        i_rev_inf = (Z>=0) & (Z<=self.t3)
        i_rev_sup = (Z>= self.t3+2*self.t2+self.h1) & (Z<=2*self.t3 + 2*self.t2 + self.h1)

        # Campo do módulo E
        E_field = np.nan*np.zeros((self.N_points, self.N_points))

        # Campo E no revestimento
        E_field[i_rev_inf] = self.E3 
        E_field[i_rev_sup] = self.E3

        # Determinando a existencia do material na alma e na aba
        for k in range(self.N_points):
            for j in range(self.N_points):
                for i in range(self.N_cel+1):
                    # Almas 
                    if ((Y[k,j]>= ((i+1)*self.e + i*self.b + (self.b-self.t1)/2)) & (Y[k,j]<= (((i+1)*self.e + i*self.b + (self.b+self.t1)/2))))&((Z[k,j]>=(self.t3+self.t2)) & (Z[k,j]<=(self.t3+self.t2+self.h1))):
                        E_field[k,j] = self.E1 
                    
                    # Aba superior 
                    if ((Y[k,j]>=((i+1)*self.e + i*self.b)) & (Y[k,j]<=(i+1)*self.e + (i+1)*self.b)) & ((Z[k,j]>=(self.htot-self.t3-self.t2)) & (Z[k,j]<= self.htot-self.t3)):
                        E_field[k,j] = self.E2
                    
                    # Aba inferior 
                    if ((Y[k,j]>=((i+1)*self.e + i*self.b)) & (Y[k,j]<=(i+1)*self.e + (i+1)*self.b)) & ((Z[k,j]>=(self.t3)) & (Z[k,j]<= self.t3+self.t2)):
                        E_field[k,j] = self.E2
                        
        # Campo de tensões
        sigma = np.zeros((self.N_points, self.N_points, self.N_points))
        for i in range(self.N_points):
            sigma[:, :, i] = E_field * (Z*My[i]/self.Iyys - Y*Mz[i]/self.Izzs) / self.E0 

        # Encontra o valor máximo ignorando NaNs
        sigma_max = np.nanmax(np.abs(sigma))

        return sigma_max, sigma
        
class torsion(aero_struct):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.Ami = (self.b+self.e)*(self.h1 + self.t2)
        self.Ci = 1 / (2*self.Ami*self.G0)
        self.betai = (self.G0/self.G1)*(self.h1/self.t1) + (self.G0/self.G2)*(self.t2/self.t1)
        self.alphai = 2*self.betai + 2*(self.G0/(self.G2*self.t2))*(self.e + self.b)
              
    def Matrix(self, **kwargs):
        """
        Construct the matrix for the given N, with coefficients Am, C, alpha, and beta.
        
        Parameters:
        - N: int, size of the matrix (number of rows/columns).
        - Am: list or array of length N, coefficients for the first row.
        - C: list or array of length N, coefficients for the remaining rows.
        - alpha: list or array of length N, coefficients for diagonal terms.
        - beta: list or array of length N-1, coefficients for off-diagonal terms (beta_n is beta[n-1]).
        
        Returns:
        - If N == 1: returns 2*Am[0] as a float.
        - If N > 1: returns the constructed N x N matrix as a numpy array.
        """
        Am = self.Ami*np.ones(self.N_cel)
        C  = self.Ci*np.ones(self.N_cel)
        alpha = self.alphai*np.ones(self.N_cel)
        beta = self.betai*np.ones(self.N_cel)

        if self.N_cel == 1:
            return 2 * Am[0]  # Special case: return a float for N=1
        
        # Initialize the matrix with zeros
        M = np.zeros((self.N_cel, self.N_cel))
        
        # Fill the first row
        M[0, :] = 2 * Am
        
        # Fill the second row (i=1 in zero-based indexing)
        M[1, 0] = C[0] * alpha[0] + C[1] * beta[0]
        M[1, 1] = -(C[1] * alpha[1] + C[0] * beta[0])
        if self.N_cel > 2:
            M[1, 2] = C[1] * beta[1]
        
        # Fill the third row (i=2 in zero-based indexing)
        if self.N_cel > 2:
            M[2, 0] = -C[1] * beta[1]
            M[2, 1] = C[1] * alpha[1] + C[2] * beta[1]
            M[2, 2] = -(C[2] * alpha[2] + C[1] * beta[1])
            if self.N_cel > 3:
                M[2, 3] = C[2] * beta[2]
        
        # Fill rows 4 to N-1 (i=3 to i=N-2 in zero-based indexing)
        for i in range(3, self.N_cel-1):
            M[i, i-1] = -C[i] * beta[i-1]
            M[i, i] = C[i] * alpha[i] + C[i+1] * beta[i]
            M[i, i+1] = -(C[i+1] * alpha[i+1] + C[i] * beta[i])
            if i+2 < self.N_cel:
                M[i, i+2] = C[i+1] * beta[i+1]
        
        # Fill the last row (i=N-1 in zero-based indexing)
        if self.N_cel > 3:
            M[-1, -3] = -C[-2] * beta[-2]
            M[-1, -2] = C[-2] * alpha[-2] + C[-1] * beta[-1]
            M[-1, -1] = -(C[-1] * alpha[-1] + C[-2] * beta[-1])
        
        return M 

    def get_ycp(self, CL):
        if CL>=0:
            return 7.31644534*np.exp(-7.21294987*CL) + 0.06457523
        else:
            return 0.01896724/(CL + 0.05777520) - CL*0.10845673 - 0.18187554
    
    def get_zcp(self, c):
        return 0.0253*c

    def torsion_moment(self, **kwargs):
        type = kwargs.get('type')
        c    = kwargs.get('c')

        if type.lower() == 'fit':
            Fz = kwargs.get('Fz')
            Fy = kwargs.get('Fy')
            CL = kwargs.get('CL')
            zcp = self.get_zcp(c)
            ycp = self.get_ycp(CL)
            return Fz*ycp - Fy*zcp 
        else:
            V = kwargs.get('V')
            Cm = kwargs.get('Cm')
            Ma = 0.5*self.rho*V**2 * c * Cm
            return Ma*np.ones(self.N_points)

    def propagate_torsion_crit(self, **kwargs):
        self.propagate_aero_critic(**kwargs)
        # Torsion distributed moment
        c = kwargs.get('c')
        typ = kwargs.get('type').lower()
        if typ == 'fit':
            self.mx_PHAA = self.torsion_moment(type = typ, Fz = self.F_PHAA[1, :], Fy = self.F_PHAA[0, :], c = c, CL = self.CL_PHAA)
            self.mx_PLAA = self.torsion_moment(type = typ, Fz = self.F_PLAA[1, :], Fy = self.F_PLAA[0, :], c = c, CL = self.CL_PLAA)
            self.mx_NHAA = self.torsion_moment(type = typ, Fz = self.F_NHAA[1, :], Fy = self.F_NHAA[0, :], c = c, CL = self.CL_NHAA)
            self.mx_NLAA = self.torsion_moment(type = typ, Fz = self.F_NLAA[1, :], Fy = self.F_NLAA[0, :], c = c, CL = self.CL_NLAA)
        else:
            Cm_PHAA = kwargs.get('Cm_PHAA')
            Cm_PLAA = kwargs.get('Cm_PLAA')
            Cm_NHAA = kwargs.get('Cm_NHAA')
            Cm_NLAA = kwargs.get('Cm_NLAA')    
            self.mx_PHAA = self.torsion_moment(type = typ, Cm = Cm_PHAA, c=c, V = self.V_PHAA)
            self.mx_PLAA = self.torsion_moment(type = typ, Cm = Cm_PLAA, c=c, V = self.Vd)
            self.mx_NHAA = self.torsion_moment(type = typ, Cm = Cm_NHAA, c=c, V = self.V_NHAA)
            self.mx_NLAA = self.torsion_moment(type = typ, Cm = Cm_NLAA, c=c, V = self.Vd)
            
        # Torsion torque
        self.Tx_PHAA = self.x_neg_integral(V = self.mx_PHAA)
        self.Tx_PLAA = self.x_neg_integral(V = self.mx_PLAA)
        self.Tx_NHAA = self.x_neg_integral(V = self.mx_NHAA)
        self.Tx_NLAA = self.x_neg_integral(V = self.mx_NLAA)

        # Angle and flux 
        self.theta_PHAA, self.q_PHAA = self.flux_angle_calc(self.Tx_PHAA)
        self.theta_PLAA, self.q_PLAA = self.flux_angle_calc(self.Tx_PLAA)
        self.theta_NLAA, self.q_NLAA = self.flux_angle_calc(self.Tx_NLAA)
        self.theta_NHAA, self.q_NHAA = self.flux_angle_calc(self.Tx_NHAA)

    def propagate_torsion(self, **kwargs):
        # Torsion distributed moment
        c = kwargs.get('c')
        typ = kwargs.get('type').lower()
        if typ == 'fit':
            Fy = kwargs.get('Fy')
            Fz = kwargs.get('Fz')
            CL = kwargs.get('CL')
            mx_point = self.torsion_moment(type = typ, Fz = Fz, Fy = Fy, c = c, CL = CL)
        else:
            Cm_point = kwargs.get('Cm')
            V = kwargs.get('V')
            mx_point = self.torsion_moment(type = typ, Cm = Cm_point, c=c, V = V)
            
        # Torsion torque
        Tx_point = self.x_neg_integral(V = mx_point)

        # Angle and flux 
        self.theta_point, self.q_point = self.flux_angle_calc(Tx_point)

    def flux_angle_calc(self, Tx):
        # Shear-stress flux 
        q = np.zeros((self.N_cel, self.N_points))

        # Torsion angle derivative
        d_dtheta = np.zeros(self.N_points)
        
        # Torsion angle
        theta = np.zeros(self.N_points)

        # Solver parameters 
        T = np.zeros(self.N_cel).T

        MA = self.Matrix()

        for i in range(self.N_points):
            if self.N_cel == 1:
                q[0, i] = Tx[i]/MA 
                d_dtheta[i] = Tx[i]*self.alphai/(4*self.Ami**2 * self.G0)
            else: 
                T[0] = Tx[i]
                q[:, i] = np.matmul(np.linalg.inv(MA), T)
                d_dtheta[i] = q[0, i]*self.Ci*self.alphai - q[1, i]*self.betai*self.Ci

        for i in range(1, self.N_points):
            theta[i] = theta[i-1] + d_dtheta[i-1]*(self.x_coord[i] - self.x_coord[i-1])

        return theta, q  