VASP K-point visualizer tool

Description:
The purpose of the VASP k-point visualizer tool is to create an interactive
three-dimensional reciprocal space plot which displays the k-points in the first Brillouin
zone or paths within it. The tool also displays information about transitions within the
electronic band structure at each k-point depending on the user-inputted transition
energy value within a desired range of error, in units of eV. The plot will alter the size
and color of the k-point if there is a probability of a transition occurring, and display
information about the occupied and unoccupied bands in which the transition occurs.
The Γ point on the plot is located at (0.5, 0.5, 0.5). The tool requires the OUTCAR and
Transmatrix files from VASP as inputs.

Usage:
The tool consists of 3 functions, which should be run in the following order:

Initialize(outcarfile, transmatrix file):
This function only needs to be run once to initialize data for every k-point.

RunTransition(outcarfile,ev,deltaE):
This function selects the k-points that satisfy the user-inputted transition energy within
the desired accuracy range and stores data about each possible transition in the
variable coordlist3.

Plot(coordlist3):
This function produces an interactive 3D plot that displays the location of every k-point,
with the k-points that satisfy the transition energy having size and color corresponding
to their transition probability.

There are 4 total user inputs to the tool:
outcarfile: OUTCAR file from VASP
transmatrixfile: Transmatrix file from VASP
ev: Desired transition energy between the occupied and unoccupied band
deltaE: Desired accuracy range, such that the tool will search for transition energies that
satisfy ≥ev-deltaE and ≤ev+deltaE.

Example Output:

![Example output](https://user-images.githubusercontent.com/79879109/271623528-793ebf20-6c1d-44bf-9adf-b643d437c8e8.png)

Installation:
The tool requires Python 3, as well as the following libraries:

-numpy

-pandas

-ipywidgets

-Plotly

Contributors:
Alp Kurbay, Brian Robinson, Andre Schleife (UIUC)

Contact:
akurbay2@illinois.edu
