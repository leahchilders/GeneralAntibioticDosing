#  General Analytic Approach to Predicting the Best Antibiotic Dosing Regimen
The code used to generate the numerical results in [this paper] is available in this repository in `input_drug_data.ipynb`. 

The code is written in a Jupyter Notebook using Python 3.12.8 and uses the SciPy library for numerical integration. The code is designed to be user-friendly and can be used to analyze any antibiotic data with a selection of PD models. The user can choose a PD model and input their own PK and PD parameters. The code outputs, among other things, the value of $\tilde{D}$, an analytic and numerical analysis of the concavity of the Hill function as performed in this paper, and plots of the dose response curve, the treatment regimens, and the bacteria population curves. The code to reproduce Figure 1 is also provided in a separate Jupyter Notebook.

### Requirements
Any modern version of Python + NumPy, SciPy, and Matplotlib.
