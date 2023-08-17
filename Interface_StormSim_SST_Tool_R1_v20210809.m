%% HEADER
%{
LICENSING:
    This code is part of StormSim software suite developed by the U.S. Army
    Engineer Research and Development Center Coastal and Hydraulics
    Laboratory (hereinafter “ERDC-CHL”). This material is distributed in
    accordance with DoD Instruction 5230.24. Recipient agrees to abide by
    all notices, and distribution and license markings. The controlling DOD
    office is the U.S. Army Engineer Research and Development Center
    (hereinafter, "ERDC"). This material shall be handled and maintained in
    accordance with For Official Use Only, Export Control, and AR 380-19
    requirements. ERDC-CHL retains all right, title and interest in
    StormSim and any portion thereof and in all copies, modifications and
    derivative works of StormSim and any portions thereof including,
    without limitation, all rights to patent, copyright, trade secret,
    trademark and other proprietary or intellectual property rights.
    Recipient has no rights, by license or otherwise, to use, disclose or
    disseminate StormSim, in whole or in part.

DISCLAIMER:
    STORMSIM IS PROVIDED “AS IS” BY ERDC-CHL AND THE RESPECTIVE COPYRIGHT
    HOLDERS. ERDC-CHL MAKES NO OTHER WARRANTIES WHATSOEVER EITHER EXPRESS
    OR IMPLIED WITH RESPECT TO STORMSIM OR ANYTHING PROVIDED BY ERDC-CHL,
    AND EXPRESSLY DISCLAIMS ALL WARRANTIES OF ANY KIND, EITHER EXPRESSED OR
    IMPLIED, INCLUDING WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY,
    NON-INFRINGEMENT, FITNESS FOR A PARTICULAR PURPOSE, FREEDOM FROM BUGS,
    CORRECTNESS, ACCURACY, RELIABILITY, AND RESULTS, AND REGARDING THE USE
    AND RESULTS OF THE USE, AND THAT THE ASSOCIATED SOFTWARE’S USE WILL BE
    UNINTERRUPTED. ERDC-CHL DISCLAIMS ALL WARRANTIES AND LIABILITIES
    REGARDING THIRD PARTY SOFTWARE, IF PRESENT IN STORMSIM, AND DISTRIBUTES
    IT “AS IS.” RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
    ERDC-CHL, THE UNITED STATES GOVERNMENT AND ITS CONTRACTORS AND
    SUBCONTRACTORS, AND SHALL INDEMNIFY AND HOLD HARMLESS ERDC-CHL, THE
    UNITED STATES GOVERNMENT AND ITS CONTRACTORS AND SUBCONTRACTORS FOR ANY
    LIABILITIES, DEMANDS, DAMAGES.

SOFTWARE NAME:
    Interface-StormSim-SST Tool

DESCRIPTION:
    Script that serves as an interface to the StormSim-SST Tool. Input 
    arguments must be specified in this script. Output arrays can be viewed 
    in this script. 

    The StormSim-SST tool is used to estimate the hazard curve
    of either a peaks-over-threshold (POT) or a raw time series dataset. A 
    description of the tool's features is included in the Quick Start Guide.

INPUT/OUTPUT ARGUMENTS:
    A description of the input and output arguments is included in the 
    Quick Start Guide.

AUTHORS:
	Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

Point of Contact: 
	Efrain.Ramos-Santiago@usace.army.mil   
	U.S. Army Engineer R&D Center Coastal & Hydraulics Laboratory                  
	Vicksburg, MS                  

HISTORY OF REVISIONS:
20201015-ERS: alpha version 0.1
20201026-ERS: alpha version 0.2: added capability to incorporate skew tides and sea level change (without steric adjustment).
20201027-ERS: alpha version 0.2.1: added types of application
20201215-ERS: alpha version 0.2.2: minor corrections.
20210113-ERS: alpha version 0.2.3: corrected rules for using of the parallel computing toolbox.
20210311-ERS: alpha version 0.3: major corrections.
20210420-ERS: alpha version 0.4: expanded capability and updated documentation.
20210809-ERS: alpha version 0.5: minor corrections and updated documentation.

**************  ALPHA VERSION  **  FOR INTERNAL TESTING ONLY  *************
%}
close all;clear;clc;

%% Input data files
load .\INPUT\InputData.mat InputData

%% Transfer of input data
% Transfer one (or multiple) dataset(s) into the structure array
% NOTE: use a for-loop for entering more than one dataset 
i=1;
input_data(i).time_values = InputData(:,1); % Time values as serial date numbers
input_data(i).data_values = InputData(:,2); % Dataset 
input_data(i).data_values2 = []; 

%% Dataset(s) information
% NOTE: one row per dataset or gauge station
staID = {'8720218','Mayport (Bar Pilots Dock)'};

%% General settings
ExecMode = 'Regular';
DataType = 'Timeseries';
use_AEP = 1;
GPD_TH_crit = 0;
prc = []; 
HC_tbl_rsp_y = []; 
apply_GPD_to_SS = 1;
apply_Parallel = 0;
path_out = []; 
flag_value = []; 
tLag = 48; 
lambda = 12;
Nyrs = []; 

%% Plot settings
yaxis_Label = {'Still Water Level (m, MSL)'}; 
yaxis_Limits = [0 5]; 

%% Skew tides settings (for ADCIRC-simulated data with skew tides)
ind_Skew = 0; 
gprMdl(1).mdl = []; 
SLC = []; 

%% USER: DO NOT CHANGE ANYTHING BELOW THIS LINE 
[HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,Removed_datasets,Check_datasets,SST_output] = StormSim_SST_Tool_R1_v20210809(input_data,flag_value,tLag,lambda,Nyrs,path_out,staID,yaxis_Label,yaxis_Limits,prc,use_AEP,GPD_TH_crit,SLC,ind_Skew,gprMdl,DataType,ExecMode,HC_tbl_rsp_y,apply_GPD_to_SS,apply_Parallel);
%% End of Script