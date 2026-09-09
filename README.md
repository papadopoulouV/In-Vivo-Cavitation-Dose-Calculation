# In-Vivo-Cavitation-Dose-Calculation
These scripts can be used to determine the total cavitation dose across multiple treatments. Example voltage and time data for one animal across multiple treatment days (pre-processed for ease of using with this script) can be used with this script.

# Requirements
The code is written in MATLAB and requires you add the natsortfiles MATLAB function to your MATPATH: https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

# Processing Steps
1) Download and save the following MATLAB scripts from GitHub and natsortfiles to your MATPATH.
	GITHUB scripts:
		cavitationAnalysis_Step1
		cavitationAnalysis_Step2
		pcdFFT4
	MATLAB  function:
		https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort

2) Save the three example data files with voltage and time data for each treatment day (D2, D3, D4) to a folder.
	This example dataset contains data from all treatment days for two animals (one from each treatment group)

3) Run the cavitationAnalysis_Step1 script and select the matfile for one treatment day.
	You will be prompted to select a location to save the cavitation data (for use in cavitationAnalysis_Step2).
	Repeat for all treatment days.

3) Run the cavitationAnalysis_Step2 script. You will be prompted to select the export files generated in (2) for D2, D3, D4.
	Compare MATLAB figures to supplemental figures for animals 2 and 16.
	T_summary table displayed in the Command window should match the result below.



| PressureAnimal | harmonicAUCnet_TotalAUC | harmonicAUCnet_TotalAUC_SD | broadbandnet_TotalAUC | broadbandnet_TotalAUC_SD | totalcavnet_TotalAUC | totalcavnet_TotalAUC_SD |
| --- | --- | --- | --- | --- | --- | --- |
| 1700_2 | 1.0639e+06 | 16304 | 9.9026e+05 | 15312 | 2.1209e+06 | 27798 |
| 700_16 | 36655 | 1801 | 11471 | 970.19 | 47159 | 2739.3 |

# License
The codes are licensed under GPL-2.0 license.

# Citation
For any utilization of the code content of this repository, the following paper needs to get cited by the user:
  
> VanTreeck KE, Liu JD, Sherman K, Tabashsum Z, Angeles-Solano M, Lifschin ZJ, Dayton PA, Rowe SE, Papadopoulou V. (2026) Ultrasound Cavitation From Low-Boiling Point Perfluorocarbon Phase Change Contrast Agents Correlates With Bacterial Burden Reduction in Antibiotic Treatment of Biofilm-Infected Murine Wounds. Ultrasound in Medicine and Biology.
