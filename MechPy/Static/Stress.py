# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:11:40 2022

@author: valte
"""
import numpy as np

class Stress():
    def __init__(self):
        return
    
    class Section():
        def __init__(self,section_type,*dims):
            self.type = section_type
            self.dims = list(dims)
            
            def rectangle(b,h):
                return {'A': b*h,
                        'I_x': (b*h**3)/12,
                        'I_y': (b**3*h)/12,
                        'I_xy': 0
                        }

            def circle(D):
                return {'A': np.pi*D**2/4,
                        'I_x': np.pi*D**4/64,
                        'I_y': np.pi*D**4/64,
                        'I_xy': 0,
                        'J': np.pi*D**4/32
                        }

            def hollow_circle(D,d):
                return {'A': np.pi*(D**2-d**2)/4,
                        'I_x': np.pi*(D**4-d**4)/64,
                        'I_y': np.pi*(D**4-d**4)/64,
                        'I_xy': 0,
                        'J': np
                        .pi*(D**4-d**4)/32
                        }
            
            if len(dims) == 1:
                self.dims.append(0)
                
            self.property_dict = {'rectangle': rectangle(self.dims[0],self.dims[1]),
                                  'circle': circle(self.dims[0]),
                                  'hollow_circle': hollow_circle(self.dims[0],self.dims[1])
                                  }
            
        def properties(self):
            match_key = ('rectangle'
                         if self.type in (1,'1','square','r','rectangle')
                         else 'circle'
                         if self.type in (2,'2','circle','round','r','c')
                         else 'hollow_circle'
                         if self.type in (3,'3','hollow_circle','hc')
                         else 'circle'
                         )
            
            return self.property_dict[match_key]
        
    def get_stress(section:Section,
                   pos,
                   normal='x',
                   full_output = False,
                   **loads):
        """
        

        Parameters
        ----------
        section : Section
            DESCRIPTION.
        pos : TYPE
            DESCRIPTION.
        normal : TYPE, optional
            DESCRIPTION. The default is 'x'.
        **loads : keyword arguments
            Cargas internas na seção:
                - N: Normal
                - V_y: Cisalhante no sentido y
                - V_z: Cisalhante no sentido z
                - T: Torque normal
                - Mf_y: Momento fletor no sentido y
                - Mf_z: Momento fletor no sentido y

        Returns
        -------
        S : TYPE
            Tensor das tensões. Nunca considera tensões no plano yz, apenas xy e xz.

        """
        axis_dict = {'x':0,
                     'y':1,
                     'z':2}
        
        local_axis = {'x': {'x_0': 'x',
                            'y_0': 'y',
                            'z_0': 'z'},
                      'y': {'x_0': 'y',
                            'y_0': 'z',
                            'z_0': 'x'},
                      'z': {'x_0': 'z',
                            'y_0': 'y',
                            'z_0': 'x'},
                      }
        
             
        pos_y = pos[0]
        pos_z = pos[1]   
        
        
        x = axis_dict[local_axis[normal]['x_0']] #x indica a orientação *local*
        y = axis_dict[local_axis[normal]['y_0']] #y indica a orientação *local*
        z = axis_dict[local_axis[normal]['z_0']] #z indica a orientação *local*
                
        sect_props = section.properties()
        
        A = sect_props['A']
        I_z = sect_props['I_x'] 
        I_y = sect_props['I_y'] 
        I_zy = sect_props['I_xy']
        J = sect_props['J']
        

        for key,val in loads.items():
            if type(val) not in [tuple,list,np.array]:
                loads[key] = [0,0,0]
                if key in ['N','T']:
                    loads[key][x] = val
                elif key in ['V_y','Mf_y']:
                    loads[key][y] = val
                elif key in ['V_z','Mf_z']:
                    loads[key][z] = val

      

        N = loads.pop('N',[0,0,0])
        V_y = loads.pop('V_y',[0,0,0])
        V_z = loads.pop('V_x',[0,0,0])

        
        T = loads.pop('T',[0,0,0])
        Mf_y = loads.pop('Mf_y',[0,0,0])
        Mf_z = loads.pop('Mf_z',[0,0,0])
        

        sig_tension = N[x]/A
        sig_flex_y = -Mf_z[z]*pos_y/I_y
        sig_flex_z = Mf_y[y]*pos_z/I_z
        

        sig_xx = sig_tension + sig_flex_y + sig_flex_z
        
        
        rho = np.sqrt(pos_z**2+pos_y**2)
        
        theta = (np.arctan(pos_y/pos_z) if pos_z != 0
                 else 0 if pos_y == 0
                 else np.pi/2 if pos_y//pos_y == 1
                 else 3*np.pi/2
                 )
        rho_dir = np.array((0,np.cos(theta),np.sin(theta)))
        T_dir = np.zeros(3)
        T_dir[x] = (T[x]//T[x] if T[x] != 0
                    else 1)
        tau_dir = np.cross(T_dir,rho_dir)
        
        tau_torsion = -T[x]*rho/J
        tau_torsion_xy = tau_dir[y]*tau_torsion
        tau_torsion_xz = tau_dir[z]*tau_torsion
        

        
        if section.type in (2,'2','circle','round','c'):
            R = section.dims[0]/2
            tau_shear_xy = 4*(V_y[y]*(R**2-pos_y**2))/(3*np.pi*R**4)
            tau_shear_xz = 4*(V_z[z]*(R**2-pos_z**2))/(3*np.pi*R**4)
        
        elif section.type in (3,'3','hollow_circle','hc'):
            R = section.dims[0]/2
            tau_shear_xy = (2*V_y[y]/A if pos_y==0
                           else 0 if pos_y==R
                           else V_y[y]/A) 
            tau_shear_xz = (2*V_z[z]/A if pos_z==0
                           else 0 if pos_z==R
                           else V_z[z]/A)
            
        else:
            a,b = section.dims[0:2]
            tau_shear_xy = V_y[y]/(2*I_y)*(b**2-pos_y**2)
            tau_shear_xz = V_z[z]/(2*I_z)*(a**2-pos_z**2)
            
        tau_xy = tau_torsion_xy+tau_shear_xy
        tau_xz = tau_torsion_xz+tau_shear_xz
        
        S = np.zeros((3,3))
        
        S[x,x] = sig_xx 
        S[x,y] = tau_xy
        S[y,x] = tau_xy
        S[x,z] = tau_xz
        S[z,x] = tau_xz
        
        stress_dict = {'sig_normal': sig_xx,
                       'sig_flex': (sig_flex_y,sig_flex_z),
                       'tau_torsion': (tau_torsion_xy,tau_shear_xz),
                       'tau_shear': (tau_shear_xy,tau_shear_xz)
                       }
        
        if full_output == True:
            return S,stress_dict
        
        else:
            return S
    
    def VonMisses(S):
        
        sig_x = S[0,0]
        sig_y = S[1,1]
        sig_z = S[2,2]
        
        tau_xy = S[0,1]
        tau_xz = S[0,2]
        tau_yz = S[1,2]


        sig_vm = (1/np.sqrt(2))*np.sqrt((sig_x-sig_y)**2
                                        +(sig_y-sig_z)**2
                                        +(sig_z-sig_x)**2 
                                        +6*(tau_xy**2+tau_xz**2+tau_yz**2)
                                        )

        return sig_vm        
    