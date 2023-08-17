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
    Interface-StormSim-JPM

DESCRIPTION:
    Main script to run the StormSim-JPM Tool. Input arguments must be 
    specified in this script. For more details, refer to the provided Quick 
    Start Guide.

INPUT/OUTPUT ARGUMENTS:
    Refer to the Quick Start Guide.

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

CONTRIBUTORS:
    Alexandros A. Taflanidis, PhD (AAT)
    Victor M. Gonzalez, PE (VMG)

POC: 
	Efrain.Ramos-Santiago@usace.army.mil   
	U.S. Army Engineer Research & Development Center 
	Coastal & Hydraulics Laboratory                  
	Vicksburg, MS                                   

HISTORY OF REVISIONS:
20210405-ERS: alpha v0.1
20210601-ERS: alpha v0.2
20210809-ERS: alpha v0.3
20210831-ERS: alpha v0.3: minor corrections.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
close all;clear;clc;

%% Input data files
load .\INPUT\Response_per_savepoint.mat Resp
load .\INPUT\ProbMass_per_Event.mat ProbMass 
load .\INPUT\VG_info.mat vg_id

%% General settings
vg_ColNum = [1 100 400];
prc = []; 
integrate_Method = 'PCHA Standard';  
ind_aep = 0;
apply_Parallel=0;
path_out = [];
HC_tbl_rsp_y=[];

%% Uncertainty settings
U_a = 0.2; 
U_r = 0.15; 
U_tide_app = 0;
U_tide = []; 
U_tide_type = [];
uncert_treatment = 'combined'; 
SLC = 0; 

%% Plot settings
plot_results = 1;
yaxis_label = 'Storm Surge (m, MSL)';
yaxis_limits = [];

%% USER: DO NOT CHANGE ANYTHING BELOW THIS LINE 
[JPM_output,HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,Removed_vg] = StormSim_JPM_Tool_R1_v20210831(Resp,ProbMass,vg_id,vg_ColNum,U_a,U_r,U_tide,U_tide_app,U_tide_type,uncert_treatment,prc,integrate_Method,path_out,yaxis_label,yaxis_limits,SLC,plot_results,ind_aep,apply_Parallel,HC_tbl_rsp_y);
%% END of Main script