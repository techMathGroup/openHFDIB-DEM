#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
#=FILE DESCRIPTION======================================================
#
# Python program to evaluate results of the settling sphere benchmark
#
# Source of data:
# [correlation] Cheng, N.-S.: Comparison of formulas for drag 
#               coefficient and settling velocity of spherical particles
#               Powder Technology 189 (2009) pp. 395-398.
#               doi:10.1016/j.powtec.2008.07.006
# [experiment]  Brown, P.P. Lawler, D.F.: Sphere drag and settling
#               velocity revisited
#               Journal of Environmental Engineering 129 (2003)
#               pp. 222-231.
#               doi:10.1061/(ASCE)0733-9372(2003)129:3(222)
# Note: experimental data are compiled from a number of previous
#       publications. some data are corrected for wall effects
#
# Code structure:
# - three general classes are prepare to define the "physics" behind the
#   simulations
#   1. particle->contains data on the settling sphere
#   2. fluid->contains data on the fluid in which the test is performed
#   3. environment->contans data common to all the tests
# Note: at the moment, only gravity and settling direction
# - one class for the [correlation] (general)
# - one class for the [experiment] (general)
# - one class for the simulation (OpenHFDIB-DEM specific)
#
# - settings for the correlation and experiment classes are based on the
#   simulations settings
#

#=LICENSE===============================================================
#  settling_sphere_benchmark_eval.py
#
#  Copyright 2025 Martin Isoz     <isozm>
#                 Ondrej Studenik <studenio>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
#=======================================================================
"""

# imports --------------------------------------------------------------
import re
import numpy as np
import os

from scipy.integrate import odeint


# -- custom classes-----------------------------------------------------
# -- A. classes usable irrespective of the solver
class particle:
    """ 
        class for the particle that is simulated, contains material and
        geometry data
    """
    
    def __init__(self,radius,density,spring,damping):
        self.radius     = radius
        self.density    = density
        self.spring     = spring                                        #spring coefficient for linear spring dashpot
        self.damping    = damping                                       #damping coefficient for linear spring dashpot
        self.volume     = 4.0/3.0*np.pi*radius**3.0
        self.mass       = self.volume*self.density
        
        self.cross_sec  = 0.5*np.pi*self.radius**2.0                    #particle cross-section in flow
        
class fluid:
    """
        class for the simulated fluid, contains material properties
    """
    
    def __init__(self,density,dyn_viscosity):
        self.density    = density
        self.mu         = dyn_viscosity
        self.nu         = self.mu/self.density
        
class environment:
    """
        class to store data on the global environment variables (forces)
    """
    
    def __init__(self,gravity,settling_direction):
        self.gravity    = gravity    
        self.settling_direction = settling_direction

class correlation:
    """
        class that contains everything related to the Re-Cd correlation
        used to validate the CFD-DEM results
    """
    
    def __init__(self,particle,fluid,environment):
        self.particle   = particle
        self.fluid      = fluid
        self.environment= environment
        
        self.solved     = False
        
        self.act_force  = []                                            #acting force during test for postprocessing
        self.Re_vals    = []                                            #Reynolds numbers during test for postprocessing
        self.Cd_vals    = []                                            #drag coeffs during test for postprocessing
        
    def Re(self,v):
        """
            for given velocity return Reynolds number
        """
        nu      = self.fluid.nu
        radius  = self.particle.radius
        v       = np.abs(v) + 1e-12
        return 2.0*radius*v/nu
        
    def Cd(self,Re):
        """
            correlation itself, kept in this form for clarity
        """
        return (24./Re)*(1.+0.27*Re)**0.43+0.47*(1.-np.exp(-0.04*Re**0.38))
        
    def acting_force(self,v,pos):
        """
            computation of force acting on the particle. has to be moved
            to a separate method for postprocessings
        """
        
        v = np.asarray(v)
        pos = np.asarray(pos)
        
        spring   = self.particle.spring
        damping  = self.particle.damping
        
        f_contact = np.where(pos <= 0, spring*pos + damping*v, 0.0)     #simplistic contact force
            
        gravity  = self.environment.gravity
        dir_sign = np.sign(self.environment.settling_direction)
        gravity *= dir_sign
        mass_s   = self.particle.mass                                   #mass of the particle
        rho_f    = self.fluid.density
        mass_f   = self.particle.volume*rho_f                           #corresponding mass of the fluid
        
        f_gravity= (mass_s - mass_f)*gravity                            #gravity/buyoancy force
        
        Re_val   = self.Re(v)
        Cd_val   = self.Cd(Re_val)
        cross_sec= self.particle.cross_sec
        f_drag   = cross_sec*rho_f*np.abs(v)*v*Cd_val
        
        return f_gravity-f_drag-f_contact
        
    def model(self,x,t):
        """
            model for the particle in fluid, allows for a very simplistic
            contact between the particle and wall
            
            x = [velocity,position]
        """
        
        v,pos    = x
        
        f_acting = self.acting_force(v,pos)
                
        dv_dx = f_acting/self.particle.mass
        ds_dx = v
        
        return [dv_dx, ds_dx]
        
    def integrator(self,init_cond,t_end,n_points):
        """
            method to integrate the model with given initial condition,
            end time, and number of points to return (uniformly
            distributed)
            
            init_cond = [initial_velocity,initial_position]
        """
        
        t = np.linspace(0.0,t_end,n_points)                             #always integrate from t = 0.0
        
        res_x, out_msg = odeint(self.model, init_cond, t,full_output=True)
                
        self.solved = True
        
        self.vel_vals = res_x[:,0]
        self.pos_vals = res_x[:,1]
        self.time_vals= t
        
    def prep_postprocessings(self):
        """
            just convert the on-the-fly constructed lists into numpy
            arrays for simpler treatment
        """
        
        if self.solved:
            v   = self.vel_vals
            pos = self.pos_vals
            self.act_force = np.array(self.acting_force(v,pos))
            self.Re_vals   = np.array(self.Re(v))
            self.Cd_vals   = np.array(self.Cd(self.Re_vals))
        else:
            raise Exception("The system needs to be solved first")
            
    def save_data(self,file_path):
        """
            save the data in a format unified with the rest of the
            postprocessings
        """
        
        self.prep_postprocessings()
        
        time_vals = self.time_vals[1::]                                 #initial condition has just too big of a drag coeff
        pos_vals  = self.pos_vals[1::]
        vel_vals  = self.vel_vals[1::]
        act_force = self.act_force[1::]
        Re_vals   = self.Re_vals[1::]
        Cd_vals   = self.Cd_vals[1::]
        
        data = np.vstack(
            [
                time_vals,
                pos_vals,
                vel_vals,
                act_force,
                Re_vals,
                Cd_vals,
            ]
        ).T
        
        header = "t\ty\tv\tf\tRe\tCd"
        np.savetxt(file_path,data,header=header,delimiter="\t",comments="")
        
class experiment:
    """
        class that contains everything related to the Re-Cd experimental
        data
    """
    
    def __init__(self,particle,fluid,environment):
        
        self.data = np.array([
            [0.00204546900784068, 11866.0089178156],
            [0.0055575377230813, 4380.30278843168],
            [0.00770746760052448, 3120.31430289127],
            [0.00827150810569307, 2980.36996408625],
            [0.00977086363698181, 2468.39660965411],
            [0.0118563299797709, 1989.90094632461],
            [0.0211837213305344, 1102.98686213207],
            [0.0224133776113131, 1108.79270643122],
            [0.023778216041798, 983.218133609588],
            [0.0238606573128367, 1014.79166386998],
            [0.025498748088403, 965.015497113737],
            [0.0306927259819384, 778.874648730925],
            [0.0327371778149632, 730.762311604458],
            [0.0347354182177486, 707.100300030685],
            [0.0349178112077496, 682.209367873707],
            [0.0361359864049185, 651.521441323929],
            [0.0372436972555025, 629.891908461823],
            [0.0390892303209355, 617.272576442247],
            [0.0423705730866492, 558.346975374594],
            [0.0486066039674844, 492.806545920061],
            [0.0499166789753353, 479.693198247336],
            [0.0530987810402053, 461.774880661631],
            [0.0566357070301909, 424.271401865121],
            [0.0594421703437372, 404.302533431214],
            [0.0665433658919302, 356.111351074039],
            [0.0693747956277823, 349.074765493573],
            [0.070264285982522, 335.947791794671],
            [0.0738948461466191, 337.820729790391],
            [0.0822792721534839, 291.287772630967],
            [0.08732309155114, 275.34331815751],
            [0.0879439582088186, 286.592074938936],
            [0.0896994333596537, 268.441541799379],
            [0.0920733292202874, 258.020482674603],
            [0.0951789200006259, 251.387660569206],
            [0.0979804468466893, 243.121244918419],
            [0.104882025843813, 228.978336158334],
            [0.105065281261169, 235.36770486247],
            [0.111483424425933, 227.732755742445],
            [0.112597012723336, 217.900640730698],
            [0.119533753193868, 202.668209731928],
            [0.123583419129259, 194.608225653723],
            [0.129916610763742, 187.433872123621],
            [0.140822451768427, 169.541203638363],
            [0.161070967098876, 152.97278062092],
            [0.166796549106829, 143.544454383215],
            [0.177906918083646, 136.373800259738],
            [0.179150854561301, 141.467263344513],
            [0.193533402406766, 124.868432315364],
            [0.193889213103715, 130.320376994377],
            [0.207909081493811, 119.798236362804],
            [0.212426512021475, 113.141568682569],
            [0.2229528609288, 111.417721369797],
            [0.258691357722499, 93.9754212090243],
            [0.266006778628873, 96.9666130258414],
            [0.272728859964002, 93.3944680231639],
            [0.279404414014782, 88.3016442251428],
            [0.290560630451821, 83.3687631748461],
            [0.311106836544361, 79.5423543726814],
            [0.316588002166498, 75.4081997510268],
            [0.334386973724359, 71.2710367226875],
            [0.359685184072818, 68.87724961098],
            [0.367500392620728, 65.4265996396933],
            [0.394972195079861, 65.3877361801436],
            [0.398350187943434, 62.0465819124161],
            [0.44491282695152, 55.0411888806209],
            [0.468036221306616, 52.1133706297119],
            [0.494207506086721, 49.2445584320696],
            [0.52819014385276, 49.4324801596181],
            [0.529611248051735, 46.2981740042282],
            [0.54461931066493, 44.9045655725577],
            [0.57902660968478, 45.3708468347602],
            [0.596075293669277, 42.6109199129009],
            [0.64018797956533, 38.9654702497964],
            [0.642650020632084, 40.7983563816721],
            [0.691623899248309, 37.2618465951626],
            [0.723299969900596, 37.8084235683141],
            [0.725246020269589, 34.9814030487596],
            [0.778977813573051, 32.6860210688732],
            [0.784017336402516, 34.7965430259514],
            [0.831200699065172, 31.0413333032064],
            [0.854410741575355, 29.0370399033862],
            [0.921170675136768, 29.3950486647797],
            [0.921170675136768, 27.758704631835],
            [0.993146939102581, 27.538945204503],
            [1.01472594983613, 25.9914093492833],
            [1.08231715465062, 24.6138632692021],
            [1.1544106302201, 23.6022314107512],
            [1.17317219417076, 22.2086573436135],
            [1.27165405376924, 22.1253359952495],
            [1.29492794834535, 20.3221364070051],
            [1.3673367204939, 20.9692435251814],
            [1.38396962994623, 18.7857943995332],
            [1.48610536984704, 18.2242262595885],
            [1.52657686503572, 17.3723220867778],
            [1.6195361150257, 17.3861488350745],
            [1.65028513647563, 16.2162343137309],
            [1.76178834436439, 15.3354069374148],
            [1.80930061000866, 15.923977784272],
            [1.84541659029827, 14.6738782769229],
            [1.88975044502019, 15.2343714133322],
            [1.94419729010411, 14.7378610006818],
            [2.09128648875364, 13.2010836242715],
            [2.19997663940691, 12.5915783272521],
            [2.38465352793061, 12.6365763889369],
            [2.56408126465704, 12.1247650190306],
            [3.02473056758774, 10.5574611579657],
            [3.05543104285586, 10.2214167496229],
            [3.21330536084292, 10.1002690968208],
            [3.34545820110019, 9.57415197270263],
            [3.6937243052164, 8.9866942268261],
            [3.91088169014498, 8.66681559669222],
            [3.93070133358267, 8.34512546947053],
            [4.06857250718552, 8.1170728925863],
            [4.14978647664066, 8.3856780627555],
            [4.19252654751918, 7.88312040843625],
            [4.69338226464326, 7.00697416766672],
            [4.95768546495015, 7.23265045410368],
            [5.118216921532, 6.92720141622767],
            [5.36554366971215, 7.11250600927409],
            [5.48118959297627, 6.69230585180858],
            [5.69330967400627, 6.55276630334357],
            [5.97733483936208, 6.30179340038182],
            [6.27352909240652, 5.95342233148521],
            [7.13712768894764, 5.5450479304592],
            [7.67414434690492, 5.13080379944845],
            [9.38745809410786, 4.40910795339828],
            [9.5914279710694, 4.17800277237127],
            [10.0993133102223, 4.41456416581126],
            [10.2696597587308, 4.24324055109419],
            [10.6797125039858, 4.01917688954708],
            [11.2996401999996, 3.97538202590604],
            [12.2340916855743, 3.58225000261246],
            [12.5519860048572, 3.80676273436524],
            [13.1226774664425, 3.62931540691943],
            [14.1988030613728, 3.46767448708558],
            [14.5073134253495, 3.2859149228946],
            [15.9117057233791, 3.33634017708667],
            [15.9235457622563, 3.16639328730628],
            [18.2785197247219, 2.79211067339707],
            [19.3889213103715, 2.83068641250691],
            [19.7067236353767, 2.73517487306702],
            [21.0759475273067, 2.76533068676122],
            [23.1636432145553, 2.48166554020204],
            [23.5304540729889, 2.38385815459712],
            [24.6964562219817, 2.41842182818623],
            [29.4883532975375, 2.24247217125596],
            [29.8871881525564, 2.07975152074457],
            [31.2000704352565, 1.99118346158458],
            [33.0111480062284, 1.97298754274766],
            [38.7860064873938, 1.66154696483151],
            [40.0569569195009, 1.7923051758498],
            [42.2684282948393, 1.66447833123015],
            [43.5998483342044, 1.75527847149688],
            [44.4823105803735, 1.66587452545857],
            [45.4488183746123, 1.74517793920367],
            [50.1993584145248, 1.55235966355617],
            [54.2673329784576, 1.50089442932488],
            [55.7452094131318, 1.4491147062004],
            [67.3482290506772, 1.34854789189168],
            [68.4355048152733, 1.2909073098378],
            [71.4623433494449, 1.32628605116715],
            [74.3133145643206, 1.26575244337626],
            [78.8852321873609, 1.28299367471285],
            [79.1214922179568, 1.21486472109778],
            [80.4598834504296, 1.16312305930164],
            [84.3162485666749, 1.15742051301084],
            [88.25690623659, 1.10639420223742],
            [94.135725662361, 1.01924906502964],
            [101.673041952807, 1.0424740174505],
            [108.584356698612, 1.00468021615896],
            [125.828918310473, 0.986095747690538],
            [130.652331877926, 0.913612331868951],
            [140.357199211905, 0.888926762246088],
            [152.685497609445, 0.883153434496971],
            [160.251506228761, 0.881750975784564],
            [163.440408996158, 0.825462615479964],
            [169.553661787917, 0.789751275038507],
            [177.389144250003, 0.821983648648174],
            [182.311334729479, 0.797372784807196],
            [188.225349226617, 0.769763670237879],
            [192.894002161198, 0.807822496701825],
            [193.066842064128, 0.735357161240833],
            [207.407865005634, 0.718671873186272],
            [224.818704047503, 0.705326544120135],
            [245.168799482026, 0.695673040113509],
            [262.027067318352, 0.703579734553748],
            [269.162918792578, 0.658883939098202],
            [286.321589947353, 0.647072565304574],
            [311.470695360891, 0.584858106081946],
            [314.669084977483, 0.641440119307019],
            [322.748938375163, 0.628463646297964],
            [345.571253269281, 0.616980470125877],
            [364.649610526438, 0.629251992420427],
            [367.177764018658, 0.597716888046654],
            [395.260097141223, 0.579124754801287],
            [428.440178961227, 0.576567895705743],
            [461.60645112804, 0.566822540699384],
            [490.37382879683, 0.546259538605236],
            [524.696711964303, 0.545024894817701],
            [549.430124901383, 0.526043091156658],
            [588.281584447436, 0.52260992061576],
            [598.646090842168, 0.506288674337833],
            [615.775885691037, 0.582202798257149],
            [616.16861128124, 0.525408829772642],
            [649.771746823622, 0.513708664961811],
            [680.170866071406, 0.516626737414872],
            [704.587190867712, 0.507175330425234],
            [736.12707191152, 0.499815649428302],
            [741.206605919459, 0.473392520480244],
            [807.75622734324, 0.465879651222809],
            [864.874515852325, 0.473079144758565],
            [950.520651577153, 0.464370650331555],
            [1010.10611960252, 0.469372342742773],
            [1039.84737801082, 0.449716520470051],
            [1091.37482723148, 0.448404683901014],
            [1094.31118900893, 0.423014329945459],
            [1145.45561079652, 0.419405393306165],
            [1164.07163713241, 0.481360318925198],
            [1164.07163713241, 0.375476065162891],
            [1221.75475817494, 0.453540256336232],
            [1324.31487778585, 0.43500725674847],
            [1413.23567586323, 0.437084644790506],
            [1440.63629968169, 0.43243793868952],
            [1470.21880807082, 0.436992751923354],
            [1541.00777463235, 0.426626272306968],
            [1686.60318101675, 0.432139458792007],
            [1828.18496968898, 0.426806792325686],
            [1900.94677476563, 0.422373131397945],
            [1993.66825882525, 0.415395116263262],
            [2148.00145942013, 0.405933455638877],
            [2328.31529501478, 0.398793702975756],
            [2409.76253990185, 0.388765211538419],
            [2427.98245977418, 0.426163124574572],
            [2443.69006122622, 0.350133752634645],
            [2448.0712900488, 0.370797644359418],
            [2606.46500105914, 0.386357675669461],
            [2965.26456725111, 0.412865738506946],
            [3266.42063814349, 0.378428353888336],
            [3598.16250567504, 0.394791916501171],
            [4093.47671307848, 0.378478115477859],
            [4296.32033987362, 0.417577573076171],
            [4446.05390211485, 0.389021903850068],
            [4485.04879523159, 0.372626328886771],
            [4921.62845479172, 0.37803543639156],
            [5213.31517262662, 0.370267489448483],
            [5283.82612249424, 0.383628528207615],
            [5560.57568533195, 0.378869945984844],
            [6987.25956459117, 0.377990157678849],
            [7313.82055921511, 0.383391012294209],
            [8371.07808014398, 0.395758886413992],
            [8803.61222154089, 0.39204614642201],
            [10157.7780319194, 0.407380049517587],
            [10868.7544506367, 0.400461456511887],
            [11381.0911697539, 0.389838346153973],
            [11673.099063229, 0.459308020999499],
            [12251.5349298261, 0.410854974112053],
            [13173.3733261112, 0.408514060667224],
            [14707.5454441506, 0.413290413305074],
            [15747.5495883495, 0.420743363938538],
            [16420.3873705606, 0.408035584183897],
            [17467.1587029171, 0.459853193488931],
            [17942.8463921943, 0.424470831513726],
            [19108.6633924153, 0.425541626629214],
            [20032.4717274609, 0.417217259302482],
            [21192.4056102305, 0.452467102845885],
            [21195.3011611963, 0.43482016105877],
            [21689.106154233, 0.401658617610004],
            [22486.9901921166, 0.452324420023159],
            [22974.5392674827, 0.435720636239032],
            [23347.9231193117, 0.405108145517519],
            [24504.8803308977, 0.453958755340234],
            [24504.8803308977, 0.407167591294298],
            [27334.2282239219, 0.433670281579829],
            [27432.3297512082, 0.404222115222988],
            [27948.0770140197, 0.451801634515514],
            [29735.1347958905, 0.454903524287998],
            [30709.5037975016, 0.420636398412134],
            [33287.4109933518, 0.425338020008995],
            [33828.4003724818, 0.45580422224941],
            [34735.4182177486, 0.429703572910985],
            [36300.514344129, 0.421828009158884],
            [37593.0269579777, 0.436594772695552],
            [39005.6526414914, 0.427439621393764],
            [39110.5979862648, 0.459413064180541],
            [41715.7668393268, 0.435910092230656],
            [42507.7954933606, 0.46074130236028],
            [45952.4735987646, 0.467971064798171],
            [49809.9508520118, 0.461580987007695],
            [50566.4727900349, 0.488696302052939],
            [51580.542806304, 0.43493514389388],
            [52701.2803759784, 0.470545969081585],
            [56098.5661733304, 0.432935061303671],
            [56933.0967400626, 0.450094920994635],
            [60441.250182061, 0.453190362557792],
            [65120.108636702, 0.440954763945154],
            [69039.1363925624, 0.43789158285929],
            [71506.9897147395, 0.469633847124598],
            [74705.4720777156, 0.440528000422674],
            [75744.9445550405, 0.486080718603803],
            [80790.338399466, 0.449188501197169],
            [86222.8388269085, 0.48731108287779],
            [87395.9818388584, 0.463586511115189],
            [89475.0407037211, 0.448844958292144],
            [93316.2077995966, 0.463817425726435],
            [96688.6541818735, 0.452231300207722],
            [104190.147207809, 0.464354427317498],
            [112130.07634271, 0.457353796449381],
            [120459.108817016, 0.458874242997533],
            [122856.168255049, 0.408559603497381],
            [128944.047212168, 0.406284024585887],
            [131039.655069979, 0.459824276542911],
            [142039.769952138, 0.457794133859348],
            [153963.288725708, 0.459889342228734],
            [166887.726467196, 0.458621892708051],
            [180897.105250902, 0.472687153841876],
            [189861.080736239, 0.431643181192369],
            [196082.500378411, 0.474376378588238],
            [205798.956172758, 0.434920225485522],
            [207464.545788608, 0.473918285046219],
        ])
        
        self.particle   = particle
        self.fluid      = fluid
        self.environment= environment
        
        self.Re_vals    = self.data[:,0]
        self.Cd_vals    = self.data[:,1]
        
    def vel_from_Re(self,Re):
        """ 
            for the given particle and fluid, return velocities requied
            to achieve the given Re
        """
        nu      = self.fluid.nu
        radius  = self.particle.radius
                
        return nu*Re/(2.0*radius)
        
    def create_time(self,time_vals):
        """
            experimental data are not time-specific, this is to create
            time for experiment based on specific case
        """
        self.time_vals = time_vals
        
    def acting_force(self):
        """
            for the given particle, fluid and environment, computes
            force acting on the particle for the compete range of
            Re available in the data
        """
        
        dir_sign = np.sign(self.environment.settling_direction)
        
        gravity  = self.environment.gravity
        gravity *= dir_sign
        mass_s   = self.particle.mass                                   #mass of the particle
        rho_f    = self.fluid.density
        mass_f   = self.particle.volume*rho_f                           #corresponding mass of the fluid
        
        f_gravity= (mass_s - mass_f)*gravity                            #gravity/buyoancy force
        
            
        f_drag = []
        vel_vals = []
        for Re_val,Cd_val in zip(self.Re_vals,self.Cd_vals):
            v   = dir_sign*self.vel_from_Re(Re_val)
            cross_sec= self.particle.cross_sec
            f_drag.append(cross_sec*rho_f*np.abs(v)*v*Cd_val)
            vel_vals.append(v)
            
        f_drag = np.array(f_drag)
        
        self.act_force = f_gravity*np.ones(f_drag.shape)-f_drag
        self.vel_vals  = np.array(vel_vals)
        
    def filter_exp_data_based_on_velocity(self,vel_filtr):
        
        time_vals   = self.time_vals
        vel_vals    = self.vel_vals
        act_force   = self.act_force
        Re_vals     = self.Re_vals
        Cd_vals     = self.Cd_vals
        
        # -- from the experimental data, select the closest velocities
        sel_exp_ind = np.array([np.abs(vel_vals - val).argmin() for val in vel_filtr])
        
        # -- filter out the data
        vel_vals    = vel_vals[sel_exp_ind]
        act_force   = act_force[sel_exp_ind]
        Re_vals     = Re_vals[sel_exp_ind]
        Cd_vals     = Cd_vals[sel_exp_ind]
        
        # -- prepare data for comparison
        aux_t       = [time_vals[0]]
        aux_pos     = [0.0]
        aux_vel     = [vel_vals[0]]
        aux_force   = [act_force[0]]
        aux_Re      = [Re_vals[0]]
        aux_Cd      = [Cd_vals[0]]
        
        vel_max     = np.max(np.abs(vel_vals))
        
        ind_init = 1
        ind_end  = 1
        for ind in range(1,len(time_vals)):
            t_init = time_vals[ind_init]
            t_end  = time_vals[ind_end]
            if np.abs(vel_vals[ind] - vel_vals[ind-1]) > 1e-4*vel_max or ind==len(time_vals)-1:
                aux_t.append(0.5*(time_vals[ind_init] + time_vals[ind_end]))
                aux_pos.append(0.0)
                aux_vel.append(0.5*(vel_vals[ind_init] + vel_vals[ind_end]))
                aux_force.append(0.5*(act_force[ind_init] + act_force[ind_end]))
                aux_Re.append(0.5*(Re_vals[ind_init] + Re_vals[ind_end]))
                aux_Cd.append(0.5*(Cd_vals[ind_init] + Cd_vals[ind_end]))
                ind_init = ind+1
                ind_end  = ind+1
            else:
                ind_end = ind+1
                
        self.time_vals = np.array(aux_t)
        self.pos_vals  = np.array(aux_pos)                              #just a placeholder
        self.vel_vals  = np.array(aux_vel)
        self.act_force = np.array(aux_force)
        self.Re_vals   = np.array(aux_Re)
        self.Cd_vals   = np.array(aux_Cd)
        
    def save_data(self,file_path,time_vals,vel_filtr):
        """
            save the data in a format unified with the rest of the
            postprocessings
        """
        
        self.acting_force()
        self.create_time(time_vals)
        self.filter_exp_data_based_on_velocity(vel_filtr)
        
        time_vals = self.time_vals[1::]                                 #initial condition has just too big of a drag coeff
        pos_vals  = self.pos_vals[1::]
        vel_vals  = self.vel_vals[1::]
        act_force = self.act_force[1::]
        Re_vals   = self.Re_vals[1::]
        Cd_vals   = self.Cd_vals[1::]
        
        data = np.vstack(
            [
                time_vals,
                pos_vals,
                vel_vals,
                act_force,
                Re_vals,
                Cd_vals,
            ]
        ).T
        
        header = "t\ty\tv\tf\tRe\tCd"
        np.savetxt(file_path,data,header=header,delimiter="\t",comments="")
    
# -- B. OpenHFDIB-DEM specific classes
class OpenHFDIB_DEM_data_reader():
    """
        class to read the case data from OpenHFDIB-DEM case directory
    """
    
    def __init__(self,case_path,environment):
        # -- global settings
        self.environment = environment
        
        # -- data files
        self.case_path = case_path
        self.controldict = self.case_path + "/" + "system/controlDict"
        self.hfdibdemdict = self.case_path + "/" + "constant/HFDIBDEMDict"
        self.transportproperties = self.case_path + "/" +  "constant/transportProperties"
        
        self.data_available = False
        
        # -- strings to look for in the log file
        self.sought_dict = {                                            #what to look for in the log
            '-- body 0 linear velocity      : (': 'velocity',
            '-- body 0 CoM                  : (': 'position',
            '-- body 0 Acting Force    : (':'force',
        }
        self.sought_strings = list(self.sought_dict.keys())
        self.time_string    = 'Time = '                                 #how are different (CFD) time levels differentiated
        
    def read_dict_value(self,file_path,var_nm):
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('//') or line.startswith('/*') or line.startswith('*') or line.startswith('//~'):
                    continue
                if var_nm in line:
                    parts = line.split()
                    if len(parts) >= 2 and parts[0] == var_nm:
                        return float(parts[-1].rstrip(';'))             #the value is the last item, assume the line ends with ";"
        return None                                                     #if not found
    
    def read_dict_vector(self,file_path,var_nm):
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('//') or line.startswith('/*') or line.startswith('*') or line.startswith('//~'):
                    continue
                if var_nm in line:
                    parts = line.split("(")[-1]
                    parts = parts.split(");")[0]
                    parts = parts.split(" ")
                    if len(parts) >= 1:
                        return np.array([float(val) for val in parts])
        return None
        
    def read_data_from_log(self,file_name="log.pimpleHFDIBFoam"):
        """ 
            method to read simulation data from a given log file
        """
        
        sought_strings = self.sought_strings
        time_string    = self.time_string
        
        with open(self.case_path + "/" + file_name, 'r') as file:
            data = file.readlines()
            file.close()
    
        outer_list = [[] for i in range(len(sought_strings)+1)]
        for i in range(len(data)):
            f_ind = data[i].find(time_string)
    
            if f_ind==0:
                outer_list[-1].append(
                    [
                        float(string_to_value) for string_to_value in re.findall(
                            r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?",
                            data[i][f_ind+len(time_string)::]
                        )
                    ]
                )
                for pos in range(len(sought_strings)):
                    string = sought_strings[pos]
                    if len(outer_list[-1]) > 0:
                        for k in range(1,2000):
                            if i-k > 0:
                                f_ind = data[i-k].find(string)    
                                if f_ind==0:
                                    outer_list[pos].append(
                                        [
                                            float(string_to_value) for string_to_value in re.findall(
                                                r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?", 
                                                data[i - k][f_ind+len(string)::]
                                            )
                                        ]
                                    )
                                    break
                            else:
                                break        
    
        sim_time  = outer_list[-1]
        sim_time  = [val[0] for val in sim_time]
        data_lists= outer_list[:-1]
        
        data_to_export = [[] for i in range(len(sought_strings))]    
    
        for pos in range(len(data_lists)):
            data_list = data_lists[pos]
            if(len(data_list) > 0):
                for data in data_list:
                    if(len(data) == 3):
                        data_to_export[pos].append(np.array(data))
                    else:
                        data_to_export[pos].append(data)
                        
        self.simulation_data = np.array(data_to_export)
        self.time_vals       = np.array(sim_time[1::])                  #I need to skip the zero time
        self.data_available  = True
        
        
    def create_particle(self):
        """
            method to read data on particle and return a particle class
        """
        
        density = self.read_dict_value(
            self.hfdibdemdict, 
            "rho"
        )
        radius = self.read_dict_value(
            self.hfdibdemdict, 
            "radius"
        )
        
        return particle(radius,density,1.0e6, 100.0)                    #spring and damping are not used in this test - hardcoded
        
    def create_fluid(self):
        """
            method to read data on fluid and return a fluid class
        """
        
        density = self.read_dict_value(
            self.transportproperties, 
            "rho"
        )
        kin_viscosity = self.read_dict_value(
            self.transportproperties, 
            "nu"
        )
        
        return fluid(density,kin_viscosity*density)
        
    def read_case_data(self,axis=0):
        """
            read the particle initial position and simulation end time
            
            Note: axis is the direction along which the particle falls
                  => axis = 0 for x or 1 for y or 2 for z
        """
        
        startPos = self.read_dict_vector(
            self.hfdibdemdict,
            "startPosition"
        )
        s0     = startPos[np.abs(axis)-1]
        
        t_end  = self.read_dict_value(
            self.controldict, 
            "endTime"
        )
        
        return s0,t_end
        
    def prep_postprocessings(self):
        
        if self.data_available:
            sim_data  = self.simulation_data
            
            dir_sign = np.sign(self.environment.settling_direction)
            
            vel_vals = dir_sign*np.linalg.norm(np.array(sim_data[0]),axis=1)
            pos_vals = np.linalg.norm(np.array(sim_data[1]),axis=1)    #just a specific test sign convention - beware
            act_force= dir_sign*np.linalg.norm(np.array(sim_data[2]),axis=1)
            
            aux_correlation = correlation(                               #auxiliary correlation class to compute Re and Cd
                self.create_particle(),
                self.create_fluid(),
                self.environment
            )
            
            self.vel_vals = vel_vals
            self.pos_vals = pos_vals
            self.act_force= act_force
            self.Re_vals   = np.array(aux_correlation.Re(vel_vals))
            self.Cd_vals   = np.array(aux_correlation.Cd(self.Re_vals))
            
        else:
            raise Exception("You need to load the data first")
        
    def save_data(self,file_path):
        """
            save the data in a format unified with the rest of the
            postprocessings
        """
        
        self.prep_postprocessings()
        
        time_vals = self.time_vals
        pos_vals  = self.pos_vals
        vel_vals  = self.vel_vals
        act_force = self.act_force
        Re_vals   = self.Re_vals
        Cd_vals   = self.Cd_vals
        
        data = np.vstack(
            [
                time_vals,
                pos_vals,
                vel_vals,
                act_force,
                Re_vals,
                Cd_vals,
            ]
        ).T
        
        header = "t\ty\tv\tf\tRe\tCd"
        np.savetxt(file_path,data,header=header,delimiter="\t",comments="")
    
        
# -- custom functions---------------------------------------------------        
        
def list_folders(directory):
    return [name for name in os.listdir(directory)
        if os.path.isdir(os.path.join(directory, name))]
        
def create_folder(folder_name):
    if(not os.path.exists(folder_name)):
        os.makedirs(folder_name)
    else:
        os.system('rm -r ' + folder_name)
        os.makedirs(folder_name)
    
# -- main script--------------------------------------------------------
def main():
    
    case_dir    = "../20_testCases/"
    results_dir = "../30_exported_results/"
    
    case_id_str = "settling"                                            #process only these directories
    skip_id_str = "_clean"                                              #skip these
        
    # -- parameters that are assumed constant
    g = 9.81                                                            #direction is driven by settling_direction
    settling_direction = -2                                             #1 - x axis; 2 - y axis; 3 - z axis; + along, - against
    bench_env = environment(g,settling_direction)                       #benchmark environment is assumed constant over the cases

    # -- list the folders in the case_directory (available data files)
    data_folders = list_folders(case_dir)    
    
    for folder in data_folders:
        
        # -- skip irrelevant folders
        if folder.find(case_id_str) < 0: continue
        if folder.find(skip_id_str) >= 0: continue
        
        print("Working on the case: " + folder)
        
        # -- prepare folder names
        case_sub_dir    = case_dir + folder 
        results_sub_dir = results_dir + folder
        
        # -- create directory to store the results in
        create_folder(results_sub_dir)
        
        # -- read the case properties
        data_reader    = OpenHFDIB_DEM_data_reader(case_sub_dir,bench_env)
        bench_particle = data_reader.create_particle()
        bench_fluid    = data_reader.create_fluid()
        s0,t_end       = data_reader.read_case_data(axis=settling_direction)
        
        # -- prepare and export data from the correlation
        bench_corr     = correlation(bench_particle,bench_fluid,bench_env)
        bench_corr.integrator([0.0,s0],t_end,1000)
        bench_corr.save_data(results_sub_dir + '/corr.dat')
        
        # -- prepare and export experimental data
        # Note: to extract relevant experimental data, I need times and
        #       velocities from the correlation
        bench_exp      = experiment(bench_particle,bench_fluid,bench_env)
        bench_exp.save_data(
            results_sub_dir + '/exp.dat',
            bench_corr.time_vals,bench_corr.vel_vals
        )
        
        # -- prepare and export the simulation data
        data_reader.read_data_from_log()
        data_reader.save_data(results_sub_dir + '/sim.dat')
            

if __name__ == "__main__":
    main()
