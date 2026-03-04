import sys
import numpy as np
from collections import deque
import openseespy.opensees as ops
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QSlider, QFrame)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Scipy Check
try:
    from scipy.linalg import solve_continuous_are
except ImportError:
    print("ERROR: 'scipy' library is required.")
    sys.exit()

# ---------------------------------------------------------
# GLOBAL PARAMETERS
# ---------------------------------------------------------
DT = 0.01   
MAX_FORCE = 2000.0   
HISTORY_WINDOW = 20.0  
MAX_INTEGRAL = 10.0  

ACADEMIC_STYLE = """
QMainWindow {
    background-color: #F5F7FA;
}
QWidget {
    color: #2C3E50;
    font-family: 'Segoe UI', sans-serif;
}
QFrame#HUD {
    background-color: #FFFFFF;
    border: 1px solid #D1D9E6;
    border-radius: 8px;
    padding: 10px;
}
QLabel#HUDLabel {
    color: #5E6C84;
    font-size: 11px;
    font-weight: bold;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}
QLabel#HUDValue {
    color: #0052CC;
    font-size: 16px;
    font-family: 'Consolas', monospace;
    font-weight: bold;
}
QPushButton {
    background-color: #0052CC;
    color: white;
    border-radius: 4px;
    padding: 8px;
    font-weight: bold;
    font-size: 12px;
    border: none;
}
QPushButton:hover {
    background-color: #0747A6;
}
QSlider::groove:horizontal {
    border: 1px solid #BBB;
    height: 4px;
    background: #E1E4E8;
    border-radius: 2px;
}
QSlider::handle:horizontal {
    background: #0052CC;
    width: 14px;
    height: 14px;
    margin: -5px 0;
    border-radius: 7px;
}
"""

class OpenSeesControlApp(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("Inverted Pendulum in OpenSeesPy- LQI Control Simulation")
        self.setGeometry(100, 100, 1280, 900) 
        self.setStyleSheet(ACADEMIC_STYLE)
        
        # --- Simulation States ---
        self.shock_steps_remaining = 0
        self.current_time = 0.0
        self.shock_magnitude_N = 5.0
        self.set_point = 0.0 
        self.integral_error = 0.0
        self.initial_disturbance = 0.02 
        
        self.is_recovering_from_shock = False
        self.shock_recovery_timer = 0
        
        # --- Data Buffers ---
        max_points = int(HISTORY_WINDOW / DT)
        self.time_data = deque(maxlen=max_points)
        self.actual_pos_data = deque(maxlen=max_points)
        self.setpoint_data = deque(maxlen=max_points)
        self.shear_force_data = deque(maxlen=max_points)
        self.bearing_disp_data = deque(maxlen=max_points)
        
        # --- Physical Properties ---
        self.L = 3.0
        self.m_cart = 10.0
        self.m_pend = 10.0
        self.g = 9.81
        
        # --- LQI Weights (Initial)---
        self.q_pos = 1200.0
        self.q_vel = 1800.0    
        self.q_theta = 2500.0
        self.q_int = 5000.0    
        self.R_val = 0.1       
        
        self.K_LQI = np.zeros(5) 
        
        # --- INIT ORDER ---
        self.build_opensees_model()
        self.init_ui()               
        self.calculate_lqi_gains()   
        
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start(int(DT * 1000))

    def calculate_lqi_gains(self):
        """Calculates LQI gains dynamically."""
        A = np.array([
            [0, 1, 0, 0, 0],
            [0, -0.1/self.m_cart, -(self.m_pend * self.g) / self.m_cart, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0.1/(self.m_cart * self.L), ((self.m_cart + self.m_pend) * self.g) / (self.m_cart * self.L), 0, 0],
            [1, 0, 0, 0, 0]   
        ])
        B = np.array([[0], [1.0/self.m_cart], [0], [-1.0/(self.m_cart*self.L)], [0]])
        
        Q = np.diag([self.q_pos, self.q_vel, self.q_theta, 100.0, self.q_int]) 
        R = np.array([[self.R_val]]) 
        
        try:
            P = solve_continuous_are(A, B, Q, R)
            self.K_LQI = np.linalg.inv(R).dot(B.T).dot(P).flatten()
            
            if hasattr(self, 'lbl_k_gains_val'):
                k_str = f"[{self.K_LQI[0]:.0f}, {self.K_LQI[1]:.0f}, {self.K_LQI[2]:.0f}, {self.K_LQI[3]:.0f}, {self.K_LQI[4]:.0f}]"
                self.lbl_k_gains_val.setText(k_str)
        except Exception as e:
            print(f"LQI Calc Error: {e}")

    def build_opensees_model(self):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        
        x1, y1 = self.L * self.initial_disturbance, self.L * np.sqrt(1 - self.initial_disturbance**2)
        
        ops.node(0, 0.0, 0.0); ops.fix(0, 1, 1, 1) 
        ops.node(1, 0.0, 0.0) 
        ops.node(2, 0.0, 0.0) 
        ops.node(3, x1, y1) 
        
        ops.fix(1, 0, 0, 1) 
        ops.fix(2, 0, 0, 0)
        ops.equalDOF(1, 2, 1, 2)
        
        ops.geomTransf('Corotational', 1)
        ops.element('elasticBeamColumn', 1, 2, 3, 100.0, 2.0e7, 0.1, 1)
        ops.mass(1, self.m_cart, 1e-6, 1e-6)
        ops.mass(3, self.m_pend, self.m_pend, 1e-6)        
        
        ops.frictionModel('Coulomb', 1, 0.05)
        ops.uniaxialMaterial('Elastic', 1, 1.0e9) 
        ops.uniaxialMaterial('Elastic', 2, 1.0e-6) 
        ops.element('singleFPBearing', 2, 0, 1, 1, 1.0e4, 1.0e6, '-P', 1, '-Mz', 2, '-orient', 0, 1, 0, 1, 0, 0)
        
        ops.fix(3, 1, 0, 0)
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        ops.load(1, 0.0, -self.m_cart * self.g, 0.0)
        ops.load(3, 0.0, -self.m_pend * self.g, 0.0)
        
        ops.system('BandGeneral'); ops.numberer('RCM'); ops.constraints('Transformation')
        ops.test('NormDispIncr', 1.0e-6, 100, 0); ops.algorithm('Newton')
        ops.integrator('LoadControl', 0.1); ops.analysis('Static')
        ops.analyze(10)
        
        ops.loadConst('-time', 0.0) 
        ops.remove('sp', 3, 1)
        ops.wipeAnalysis()
        
        ops.wipeAnalysis()
        ops.timeSeries('Constant', 2)
        ops.timeSeries('Constant', 3)
        
        ops.constraints('Transformation') 
        ops.numberer('RCM')            
        ops.system('BandGeneral')
        ops.test('EnergyIncr', 1.0e-6, 1250, 0) 
        ops.algorithm('KrylovNewton')        
        ops.integrator('TRBDF2') 
        ops.analysis('Transient')      

    def init_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        
        left_panel = QWidget(); left_layout = QVBoxLayout(left_panel)
        
        try:
            plt.style.use('seaborn-v0_8-whitegrid')
        except:
            available = plt.style.available
            if 'seaborn-whitegrid' in available:
                plt.style.use('seaborn-whitegrid')
        else:
            plt.style.use('bmh') 
        self.fig = Figure(facecolor='#F5F7FA')  
        self.canvas = FigureCanvas(self.fig)
        
        self.fig.subplots_adjust(top=0.95, bottom=0.08, left=0.10, right=0.95, hspace=0.3)
        gs = self.fig.add_gridspec(4, 1) 
        
        self.ax_anim = self.fig.add_subplot(gs[0:2, 0])
        self.ax_anim.set_facecolor('#FFFFFF') 
        self.ax_anim.set_xlim(-4.0, 4.0); self.ax_anim.set_ylim(-1.0, 3.5); self.ax_anim.set_aspect('equal')
        self.ax_anim.set_xlabel("Horizontal Position (m)", fontsize=10, color='#333')
        self.ax_anim.set_ylabel("Height (m)", fontsize=10, color='#333')
        self.ax_anim.grid(True, linestyle=':', alpha=0.6) 
        self.ax_anim.annotate('Frictional Surface', 
                              xy=(0.0, -0.42),            
                              xytext=(-2.5, -0.98),       
                              arrowprops=dict(facecolor='#34495E', shrink=0.05, width=1.5, headwidth=8),
                              fontsize=10, color='#34495E', fontweight='bold', ha='center')
        # --------------------------------------------------
        

        ground = Rectangle((-5.0, -0.42), 10.0, 0.15, facecolor='#2C3E50', edgecolor='none', zorder=1)
        self.ax_anim.add_patch(ground)

        self.line_pendulum, = self.ax_anim.plot([], [], '-', lw=4, color='#34495E', zorder=2) 
        self.cart_marker, = self.ax_anim.plot([], [], 's', markersize=30, color='#95A5A6', zorder=3) 
        self.top_marker, = self.ax_anim.plot([], [], 'o', markersize=14, color='#C0392B', zorder=4) 
        self.quiver_force = self.ax_anim.quiver([], [], [], [], color='#E74C3C', scale=5000, width=0.015, zorder=5)

        self.ax_graph = self.fig.add_subplot(gs[2, 0])
        self.ax_graph.set_facecolor('#FFFFFF')
        self.ax_graph.set_ylabel("Position (m)", fontsize=10, color='#333')
        self.ax_graph.set_xlabel("Time (sec)", fontsize=10, color='#333')
        self.ax_graph.grid(True, linestyle=':', alpha=0.6)
        
        self.line_actual, = self.ax_graph.plot([], [], '-', color='#C0392B', lw=2.5, label='Actual') 
        self.line_setpoint, = self.ax_graph.plot([], [], '--', color='#2980B9', lw=2, label='Target') 
        
        self.ax_graph.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.9)
        

        self.ax_hyst = self.fig.add_subplot(gs[3, 0])
        self.ax_hyst.set_facecolor('#FFFFFF')
        self.ax_hyst.set_xlabel("Displacement (m)", fontsize=10, color='#333')
        self.ax_hyst.set_ylabel("Shear Force (N)", fontsize=10, color='#333')
        self.ax_hyst.grid(True, linestyle=':', alpha=0.6)
        self.line_hyst, = self.ax_hyst.plot([], [], '-', color='#8E44AD', lw=2) 
        
        left_layout.addWidget(self.canvas)
        main_layout.addWidget(left_panel, stretch=3) 

    
        right_panel = QWidget(); right_layout = QVBoxLayout(right_panel); right_panel.setFixedWidth(400)
        
        hud_container = QFrame(); hud_container.setObjectName("HUD")
        hud_grid = QVBoxLayout(hud_container)
        self.lbl_time_val = self.create_hud_row(hud_grid, "TIME", "0.00 s")
        self.lbl_pos_val = self.create_hud_row(hud_grid, "POSITION", "0.000 m")
        self.lbl_theta_val = self.create_hud_row(hud_grid, "ANGLE", "0.000 rad")
        self.lbl_force_val = self.create_hud_row(hud_grid, "CART FORCE", "0.0 N")
        self.lbl_shear_val = self.create_hud_row(hud_grid, "SHEAR FORCE", "0.0 N")
        self.lbl_k_gains_val = self.create_hud_row(hud_grid, "LQI GAINS (K)", "[0, 0, 0, 0, 0]")
        self.lbl_k_gains_val.setStyleSheet("font-size: 11px; color: #555;")
        right_layout.addWidget(hud_container)

        # --- CONTROLS ---
        self.create_slider(right_layout, 'Setpoint', 'TARGET (m)', -200, 200, 0, is_setpoint=True)
        
        self.create_slider(right_layout, 'q_pos', 'Q-POS (Accuracy)', 10, 10000, self.q_pos)
        self.create_slider(right_layout, 'q_vel', 'Q-VEL (Damping)', 10, 5000, self.q_vel)
        self.create_slider(right_layout, 'q_theta', 'Q-THETA (Balance)', 100, 25000, self.q_theta)
        
        self.create_slider(right_layout, 'q_int', 'Q-INT (Steady Error Correction)', 10, 10000, self.q_int)
        
        self.create_slider(right_layout, 'R_val', 'R (Control Cost)', 1, 500, self.R_val*100, scale=100.0)
        
        btn_shock = QPushButton("APPLY SHOCK IMPULSE"); btn_shock.clicked.connect(self.trigger_shock)
        right_layout.addWidget(btn_shock)
        
        btn_reset = QPushButton("REBOOT SYSTEM"); btn_reset.setStyleSheet("background-color:#7F8C8D;"); btn_reset.clicked.connect(self.reset_simulation)
        right_layout.addWidget(btn_reset)

        right_layout.addStretch()
        attribution = QLabel(
            "Written by\n"
            "Huseyin Cilsalar, PhD\n"
            "Yozgat Bozok University\n"
            "Dept. of Civil Engineering\n"
            "huseyin.cilsalar@bozok.edu.tr"
        )
        attribution.setFont(QFont("Segoe UI", 11, QFont.Bold)); attribution.setStyleSheet("color: #7F8C8D;"); attribution.setAlignment(Qt.AlignCenter)
        right_layout.addWidget(attribution)
        
        main_layout.addWidget(right_panel, stretch=1)

    def create_hud_row(self, layout, title, init_val):
        row = QWidget(); rl = QHBoxLayout(row); rl.setContentsMargins(0,2,0,2)
        lt = QLabel(title); lt.setObjectName("HUDLabel")
        lv = QLabel(init_val); lv.setObjectName("HUDValue"); lv.setAlignment(Qt.AlignRight)
        rl.addWidget(lt); rl.addWidget(lv); layout.addWidget(row)
        return lv

    def create_slider(self, layout, name, label, min_v, max_v, init_v, is_setpoint=False, scale=1.0):
        c = QWidget(); v = QVBoxLayout(c)
        l_row = QHBoxLayout(); ln = QLabel(label); lv = QLabel(f"{init_v/scale if not is_setpoint else init_v/100.0:.2f}"); lv.setAlignment(Qt.AlignRight)
        l_row.addWidget(ln); l_row.addWidget(lv)
        
        s = QSlider(Qt.Horizontal); s.setMinimum(min_v); s.setMaximum(max_v); s.setValue(int(init_v))
        
        if is_setpoint: 
            s.valueChanged.connect(lambda val, l=lv: self.update_setpoint(val, l))
        else: 
            s.valueChanged.connect(lambda val, n=name, l=lv, sc=scale: self.update_lqi_weight(val, n, l, sc))
            
        v.addLayout(l_row); v.addWidget(s); layout.addWidget(c)

    def update_lqi_weight(self, val, name, label, scale): 
        real_val = val / scale
        setattr(self, name, float(real_val))
        label.setText(f"{real_val:.2f}") 
        self.calculate_lqi_gains()       

    def update_setpoint(self, val, label): 
        self.set_point = val/100.0; label.setText(f"{self.set_point:.2f}")

    def trigger_shock(self): 
        self.shock_steps_remaining = 10; self.integral_error = 0.0; self.is_recovering_from_shock = True; self.shock_recovery_timer = 1.0
    def reset_simulation(self): 
        self.current_time = 0.0; self.integral_error = 0.0; self.build_opensees_model()

    def update_simulation(self):
        try: 
            cx = ops.nodeCoord(1, 1) + ops.nodeDisp(1, 1); cv = ops.nodeVel(1, 1)
            tx = ops.nodeCoord(3, 1) + ops.nodeDisp(3, 1); ty = ops.nodeCoord(3, 2) + ops.nodeDisp(3, 2); tv = ops.nodeVel(3, 1)
            sf = ops.eleResponse(2, 'localForce')[1]
        except: return
            
        self.current_time += DT
        theta = np.arctan2(tx - cx, ty)
        tv_approx = (tv - cv) / (self.L * np.cos(theta)) if abs(np.cos(theta)) > 0.01 else 0.0
        err = cx - self.set_point
        

        if not self.is_recovering_from_shock:
            self.integral_error += (cx - self.set_point) * DT 
            self.integral_error = np.clip(self.integral_error, -MAX_INTEGRAL, MAX_INTEGRAL)
        else:
            self.shock_recovery_timer -= DT
            if self.shock_recovery_timer <= 0: self.is_recovering_from_shock = False

        Z = np.array([err, cv, theta, tv_approx, self.integral_error if not self.is_recovering_from_shock else 0.0])
        u = np.clip(float(-np.dot(self.K_LQI, Z)), -MAX_FORCE, MAX_FORCE)
        
        ops.remove('loadPattern', 2); ops.pattern('Plain', 2, 2); ops.load(1, u, 0.0, 0.0)
        if self.shock_steps_remaining > 0:
            ops.remove('loadPattern', 3); ops.pattern('Plain', 3, 3); ops.load(3, self.shock_magnitude_N, 0.0, 0.0); self.shock_steps_remaining -= 1

        if ops.analyze(1, DT) != 0: self.reset_simulation(); return
        
        # HUD & Plots
        self.lbl_time_val.setText(f"{self.current_time:.2f} s"); self.lbl_pos_val.setText(f"{cx:+.3f} m")
        self.lbl_theta_val.setText(f"{theta:+.3f} rad"); self.lbl_force_val.setText(f"{u:+.1f} N"); self.lbl_shear_val.setText(f"{sf:+.1f} N")
        
        self.time_data.append(self.current_time); self.actual_pos_data.append(cx); self.setpoint_data.append(self.set_point)
        self.bearing_disp_data.append(cx); self.shear_force_data.append(sf)
        
        self.line_pendulum.set_data([cx, tx], [0.0, ty]); self.cart_marker.set_data([cx], [0.0]); self.top_marker.set_data([tx], [ty])
        self.line_actual.set_data(self.time_data, self.actual_pos_data); self.line_setpoint.set_data(self.time_data, self.setpoint_data)
        self.quiver_force.set_offsets([cx, 0.0]) 

        self.quiver_force.set_UVC([u*2250/40], [0])     
        
        
        
        self.ax_graph.set_xlim(max(0, self.current_time - HISTORY_WINDOW), max(HISTORY_WINDOW, self.current_time))
        self.ax_graph.set_ylim(-0.75, 0.75) 
        
        self.line_hyst.set_data(self.bearing_disp_data, self.shear_force_data)
        if len(self.bearing_disp_data) > 1:
            self.ax_hyst.set_xlim(min(self.bearing_disp_data)-0.1, max(self.bearing_disp_data)+0.1)
            self.ax_hyst.set_ylim(min(self.shear_force_data)-10, max(self.shear_force_data)+10)
        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv); ex = OpenSeesControlApp(); ex.show(); sys.exit(app.exec_())
