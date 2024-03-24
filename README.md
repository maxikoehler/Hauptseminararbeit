# 1 Paper

## 1.1 Supervisors and general conditions

Seminar paper at the Chair of Electrical Energy Systems (EES) at FAU Erlangen-Nuremberg. Seminar Electrical Energy Supply (Elektrische Energieversorgung).<br>
Supervisor: M.Sc. Ilya Burlakin<br>
Seminar leader: Dr. Gert Mehlmann<br><br>
Date of submission and presentation: 28th March 2024

## 1.2 Titel

> _"Critical clearing time of synchronous generators"_

## 1.3 Task definition

The critical clearing time (CCT) is an essential parameter in power system stability analysis.
For example, in the case of synchronous generators, the CCT determines the maximum fault-clearing time a generator can withstand without losing synchronism. This seminar will introduce the concept of CCT computing.
We will discuss the factors influencing CCT, such as generator characteristics, system parameters, and fault type, and explore the methods used to calculate CCT in practical power system analysis.

Set framework:
- Swing equation of synchronous generators
- Solving the Swing equation with the help of Python
- Equal-area criterion (Derivation of the equations)
- Simulation of a fault
- Comparison between analytical and (numerical) simulation results

## 1.4 Abstract

The goal of this Student Research Paper is the determination of the critical clearing time (CCT), looking at a synchronous generator (SG) in a simplified single machine infinite bus (SMIB) model. For this, a simple three-phase fault scenario is applied, and the generator swing equation is solved in the time domain with a Python-integrated solver. A function is implemented to calculate the CCT numerically, the result is compared to the analytical solution. Two additional interruption scenarios are constructed, simulated and evaluated. These three cases illustrate the transient stability of a generator against an infinite bus bar (IBB), especially in a visual context with selected plots. The numerical algorithm shows satisfying results compared to the analytical. Limitations rely on the complexity of the considered electrical network, also in the interrupting scenario, and the missing possible machine interaction. Further, the damping in the system is neglected completely.

# 2 Python Code

This section shall be used as a mini-documentation of the included python code. Some general usage and preface, followed by some information about the single (created and used) modules.

## 2.1 General

Firstly, the code is historically grown, and due to the short time frame not cleaned up, organized or fully commented. It is devided into a main module (`smib_model.py`), containing all relevant functionalities for generating a time domain solution and determing the critical clearing time and the critical power angle. The other modules are looking into one specific assessment using the simulation method from the main module (faults and comparisons). If one is interested in the result of the fault 1 parameter set, running only the module `fault1.py` in the same directory as the main module shall be satisfying.

## 2.2 `smib_model.py`

### 2.2.1 Purpose and structure

### 2.2.2 Functions

## 2.3 `fault1.py`, `fault2.py`, `fault3.py`

### 2.3.1 Purpose and structure

### 2.3.2 Functions

## 2.4 `parameter_comparison.py` and `comparison_alg-vs-nonalg.py`

### 2.4.1 Purpose and structure

### 2.4.2 Functions
