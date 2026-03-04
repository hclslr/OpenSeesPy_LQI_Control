# OpenSeesPy_LQI_Control
This repository demonstrates a high-fidelity structural control simulation framework. It bridges the gap between Structural Engineering (OpenSeesPy) and Control Theory (LQI) to manage chaotic dynamics in real-time.
Project OverviewThis project simulates a Single Inverted Pendulum mounted on a frictional cart system. While the pendulum represents a simplified structural mass, the base interaction is modeled as a Flat Slider Bearing with Coulomb friction—a common element in seismic isolation design.

Key Highlights:

Physics Engine: Utilizes OpenSeesPy for nonlinear transient analysis.Control Algorithm: Implements a Linear Quadratic Integral (LQI) controller to ensure zero steady-state error despite non-linear friction.

Real-Time Monitoring: A custom PyQt5 dashboard for dynamic parameter tuning and data visualization.Structural Insight: Real-time Force-Displacement (Hysteresis) plotting to monitor energy dissipation. Installation & SetupTo ensure portability and prevent common Matplotlib version errors, follow these steps:Clone the repository:Bashgit clone https://github.com/hclslr/OpenSees_LQI_Control.git
cd OpenSees_LQI_Control

Install dependencies:Bashpip install -r requirements.txt

Technical Details

Control StrategyThe system manages five distinct states to achieve stability:Cart Position & Velocity: ($x, \dot{x}$)Pendulum Angle & Angular Velocity: ($\theta, \dot{\theta}$)Integral of Position Error: This state is crucial for overcoming the "stick-slip" behavior of the Coulomb friction model used in the OpenSees simulation.

Structural ModelingThe model is built with a singleFPBearing element in OpenSees, representing a frictional surface with a defined friction coefficient. The LQI controller dynamically calculates the required force ($u$) to counteract both gravity and frictional resistance.

About the AuthorHuseyin Cilsalar, PhD Department of Civil Engineering, Yozgat Bozok University Expertise in Structural Dynamics, Seismic Dampers, and Mechatronics. Contact: huseyin.cilsalar@bozok.edu.tr

LicenseThis project is licensed under the MIT License.
