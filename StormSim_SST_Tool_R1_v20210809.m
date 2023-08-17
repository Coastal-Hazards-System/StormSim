%% StormSim_SST_Tool.m
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
    StormSim-SST-Tool (Statistics)

GENERAL DESCRIPTION:
   This script receives and processes the input arguments to the StormSim's
   Stochastic Simulation Technique (SST) Tool, executes the diferent methodologies,
   and conveys the output arguments to the workspace window in addition to
   save them in a .mat file.

   The tool can be used to estimate the hazard curve (HC) in terms annual
   exceedance probability (AEP) or annual exceedance frequency (AEF) from one
   or many datasets consecutively (either raw time series or POT samples).

   *** NEED more info about the methods and logic inside the SST (i.e., MRL  ***
   *** method, bootstrap process, combination of emp with GPD fit)           ***

   This script enables the communication between the Peaks-Over-Threshold
   (POT) methodology (StormSim_POT.m) and the SST methodology. The SST is
   executed with one of the following three scripts:
     1) StormSim_SST_Fit.m: when GPD_TH_crit = 0, executes the SST by evaluating
        each of the GPD thresholds identified by the MRL method and produces
        a sample plot of the bootstrap process.
  
     2) StormSim_SST_Fit_Simple.m: when GPD_TH_crit = 1 or 2, executes the SST
        using the GPD threshold identified with the MRL method criterion specified
        by the user. No sample plot of the bootstrap process is produced.

     3) StormSim_SST_Fit_SimplePar.m: this is a parallelized version of the
        script StormSim_SST_Fit_Simple.m for faster execution. Will only be
        used when the Parallel Computing Toolbox is installed in the User's
        MATLAB software.

DESCRIPTION: MODES OF EXECUTION

   The tool offers two modes or types of execution for the user to select:
   'Regular' and 'Fast'. This modes are specified in ExecMode
   (see Section INPUT ARGUMENTS) and have the following implications:
    1) 'Regular': under this mode, the tool will return to the User all
       available outputs and plots for each input dataset (see Section OUTPUT
       ARGUMENTS). This provides the opportunity to test different setting
       values and see how the HC changes, so the User can decide which values
       are best per dataset.
    2) 'Fast': this mode limits the output information to increase
       execution time and decrease output size for a straight forward application
       of several datasets. The user is adviced to use this mode after selecting
       the best settings under the 'Regular' mode.

   The tool can also plot the results depending on the execution type selected
   by the User, as follows:
    1) 'Regular': when ExecMode = 'Regular', the tool will produce histogram
       plots of the GPD shape parameter values, plots of the HC for each GPD
       threshold identified by the MRL method, and a stacked plot of the
       MRL method results.

    2) 'Fast': when ExecMode = 'Fast', no plots
       will be produced for faster evaluation of the input datasets.

   The tool will store all available outputs on a .mat file named StormSim_SST_output.mat.
   This file will be saved inside the folder directory specified in path_out
   (see Section INPUT ARGUMENTS). For a detailed list and description of the
   outputs, see the Section OUTPUT ARGUMENTS. Some of the results are organized
   in a structure array with a different set of fields depending on the mode
   of execution selected by the User, as follows:
    1) 'Regular': when ExecMode = 'Regular', the output structure
       array will contain the station information (staID); the record length
      (RL); the POT sample (POT); the outputs of the MRL method (MRL_out);
       the HC as plotted (HC_plt) for the values in HC_plt_x; the summarized
       HC (HC_tbl) for the values in HC_tlb_x; the summarized HC (HC_tbl_rsp_x)
       for the values in HC_tbl_rsp_y; the empirical HC (HC_emp); and an error
       message when the tool failed to evaluate a dataset (ME).

    2) 'Fast': when ExecMode = 'Fast', the output
       structure array will only contain the station information (staID);
       the HC as plotted (HC_plt) for the values in HC_plt_x; the summarized
       HC (HC_tbl) for the values in HC_tlb_x; the summarized HC (HC_tbl_rsp_x)
       for the values in HC_tbl_rsp_y; and an error message when the tool
       failed to evaluate a dataset (ME).
      
DESCRIPTION: TYPES OF APPLICATIONS

   Current design of the StormSim-SST Tool enables the evaluation of
   general and case-specific applications. A general application considers the
   use of raw time series datasets or POT datasets of any type of parameter.
   Refer to Section INPUT ARGUMENTS for a description of each input. Settings
   for general applications are as follows:

    G1) For raw time series dataset(s) with regular setting:
       DataType = 'Timeseries' and ExecMode = 'Regular'. This is
       intended for a detailed evaluation of a dataset, when there is no
       knowledge of the input arguments values.
        
    G2) For raw time series dataset(s) with mass-production setting:
       DataType = 'Timeseries' and ExecMode = 'Fast'.
       This is intended for the evaluation of hundreds of datasets, after
       the input arguments values have been defined.

    G3) For POT dataset(s) with regular setting: DataType = 'POT' and
       ExecMode = 'Regular'. Same as #1 but using POT dataset.

    G4) For POT dataset(s) with mass-production setting, DataType = 'POT'
       and ExecMode = 'Fast'. Same as #2 but using POT dataset.

   Available case-specific applications are described below:
    C1) POT dataset of storm surge with skew tides added.
       The following must be specified for this application:
       - set input argument ind_Skew = 1.
       - storm surge POT dataset values cannot include tides.
       - storm surge values may/may not include sea level change.
         > If SLC is included, also specify that value in the input argument SLC.
         > If SLC is not included, set input argument SLC = 0.
       - a second POT dataset of storm surge with tides must be provided
         for replacement of the empirical distribution. Those need to be
         entered in the input argument input_data.data_values2.
       - a Gaussian regression model (GRM) per dataset. Those need to be
         entered in the input argument gprMdl.mdl.
    C2) Additional cases will be added in future versions.

INPUT ARGUMENTS:
  - DataType: indicator for specifying the type of input dataset. Current
      options are:
      > 'POT' for a Peaks-Over-Threshold sample
      > 'Timeseries' for a raw time series data set
  - ind_Skew: indicator for computing/adding skew tides to the storm surge.
      This applies when the input is a storm surge POT dataset with/without SLC.
      Use as follows:
      > 1 for the tool to compute and add the skew tides
      > 0 otherwise
      Example: ind_Skew = 1;
  - input_data: raw time series datasets or POT samples to be evaluated;
      specified as a structure array with fields time_values, data_values
      and data_values2. Each record of input_data must correspond to an input
      dataset. Use as follows:
      > When DataType = 'Timeseries':
        input_data.time_values: must have the timestamps as serial date
          numbers (see 'datenum' in the MATLAB Documentation)
        input_data.data_values: must have the time series data values
      > When DataType = 'POT':
          input_data.time_values can be empty []
          input_data.data_values must have the POT values
          input_data.data_values2 can be empty []
      > When ind_Skew = 1: (evaluating POT storm surge with skew tides)
          input_data.time_values can be empty []
          input_data.data_values must have the POT storm surge values without tides
          input_data.data_values2 must have the POT storm surge values with tides
  - ExecMode: mode or type of execution. Available options are 'Regular'
      or 'Fast'; specified as a character vector.
      Example: ExecMode = 'Regular'.
  - flag_value: any flag value to search and remove from the input dataset;
      specified as a scalar. Leave empty [] otherwise. Example: flag_value = -999;
  - tLag = inter-event time in hours; specified as a positive scalar or
      vector. Use as follows:
      > When DataType = 'POT': leave empty [].
      > When DataType = 'Timeseries': Enter one value per dataset or one
        value for all datasets. Examples: tLag = [48 24]; tLag = 48;
  - lambda = mean annual rate of events (events/year); specified as a positive
      scalar or vector.  Use as follows:
      > When DataType = 'POT': leave empty [].
      > When DataType = 'Timeseries': Enter one value per dataset or one
        value for all datasets. Examples: lambda = [12 6 3]; lambda = 12;
  - Nyrs: record length in years; specified as a positive scalar or vector.
      Use as follows:
      > When DataType = 'POT': enter one value per dataset or one value
        for all datasets. Examples: Nyrs = [50 15 30]; Nyrs = 50;
      > When DataType = 'Timeseries': leave empty [] since the tool will
        compute the effective duration.
  - path_out: path to output folder; specified as a character vector. Leave
      empty [] to apply default: '.\SST_output\'
  - staID: gauge station information for each input dataset; specified as a cell array with format:
      Col(01): gauge station ID number; as a character vector. Example: '8770570'
      Col(02): station name (OPTIONAL); as a character vector. Example: 'Sabine Pass North TX'
  - yaxis_Label: parameter name/units/datum for label of the plot y-axis;
      specified as a cell array of character vectors. Available options for this input are as follows:
      > Enter one label that applies to all input datasets. Example : yaxis_Label = {'Still Water Level (m, MSL)'};
      > Enter one label per input dataset. Example for three : yaxis_Label = [{'Still Water Level (m, MSL)'},{'Significant wave height (m)'},{'Peak wave period (s)'}];
      > When ExecMode = 'Fast', leave empty [].
      > When ind_Skew = 1 and the input data is of storm surge, set
        yaxis_Label = 'Still water level' with associated units and datum.
        Example: 'Still Water Level (m, MSL)'
  - yaxis_Limits: lower and upper limits for the plot y-axis; specified as a
      two-column numerical array. Available options for this input are as follows:
      > Enter a two-value vector that applies to all input datasets. Example: yaxis_Limits = [0 10];
      > Enter as a two-column matrix with one row per input dataset. Example for three input datasets: yaxis_Limits = [[0 10];[2 10];[0 20]];
      > When ExecMode = 'Fast', leave empty [].
      > Leave empty [] to apply default values.
  - prc: percentage values for computing the percentiles; specified as a
      scalar or vector of positive values. Leave empty [] to apply default
      values 2.28%, 15.87%, 84.13%, 97.72%. User can enter 1 to 4 values.
      Example: prc = [2 16 84 98];
  - use_AEP: indicator for expressing the hazard as AEF or AEP. Use 1 for
      AEP, 0 for AEF. Example: use_AEP = 1;
  - GPD_TH_crit: indicator for specifying the GPD threshold option of the
      Mean Residual Life (MRL) selection process; specified as a scalar. Use as follows:
      > GPD_TH_crit = 0: to evaluate all thresholds identified by the MRL method (one per criterion)
      > GPD_TH_crit = 1: to only evaluate the MRL threshold selected by the lambda criterion
      > GPD_TH_crit = 2: to only evaluate the MRL threshold selected by the minimum error criterion
  - SLC: magnitude of the sea level change implicit in the storm surge input
      dataset; specified as a positive scalar.
      > When ind_Skew = 1: Must have same units as the input POT storm surge dataset. Example: SLC = 1.8;
      > When ind_Skew = 0: Leave empty [].
  - gprMdl: Gaussian process regression (GPR) model created with the MATLAB
function 'fitrgp'; specified as a structure array with field 'mdl'.
      > When ind_Skew = 1: Field 'mdl' must contain one GPR object per
        input dataset included in input_data. Refer to MATLAB Documentation
        for 'fitrgp'. Train the GPR with storm surge as a predictor and skew
        tides as the response.
      > When ind_Skew = 0: Leave empty [].
  - HC_tbl_rsp_y: response values used to summarize the HC; specified as a
      numerical vector of positive values. Set as empty [] to apply default
      values (0 to 20). Example: HC_tbl_rsp_y = [0:05:10];
  - apply_GPD_to_SS: indicates if the hazard of small POT samples should evaluated
      with either the empirical or the GPD. A small sample has a sample size < 20 events
      and a record length < 20 years. Use 1 for GPD, 0 for empirical. Example: apply_GPD_to_SS = 1;
  - DELETED K_par_restriction: options to restrict the range of values of GPD shape
      parameter from the bootstrap process. Available options are:
        1 = range -.7 to .1
        2 = range -.8 to .2
        3 = apply fillsoutlier() to remove outliers from bootstrap results;
            TH_otlr cannot be [].
  - DELETED TH_otlr: 'ThresholdFactor' for the fillsoutlier function used to modify
      the GPD shape parameter collection resulting from the bootstrap process.
      Specified as a nonnegative scalar or vector. Refer to MATLAB Documentation
      for 'fillsoutlier', 'ThresholdFactor' option. Enter one value per dataset or one
      value for all datasets. Leave empty [] to apply default value: 10.
      Examples: TH_otlr = [1 2 3]; TH_otlr = 3;
  - apply_Parallel: Enter one (1) to run tool in parallel; enter zero (0)
      otherwise. The tool will execute in parallel ONLY when GPD_TH_crit =
      0, and the Parallel Computing Toolbox is installed.


OUTPUT ARGUMENTS:
  - HC_plt_x: predefined vector of probabilities used to plot the HC. The type is:
     > AEP when use_AEP = 1
     > AEF when use_AEP = 0

  - HC_tbl_x: predefined vector of probabilities used to summarize the HC. The type is:
     > AEP when use_AEP = 1
     > AEF when use_AEP = 0

  - HC_tbl_rsp_y: vector of response values used to summarize the HC.

  - Removed_datasets: Message with a list of the input datasets not evaluated due to one of the following reasons:
     > Dataset resulted empty after removing flag, zeros, NaN and/or Inf values.
     > Dataset has less than 3 unique values

  - Check_datasets:
     > When ind_Skew = 1: Message with a list of the input datasets to which the skew tides were not applied.
       For example: 'These stations were evaluated without applying skew tides since: ind_Skew = 1
            but entered invalid values for either the second POT dataset (with tides) or the GPR model'.
     > When ind_Skew = 0: will be empty.

  - SST_output: structure array containing the output data of the StormSim-SST
     Tool, with the following fields:
     > staID: identifier of the gauge station as specified by the user
     > RL: record length in years
     > POT:
         When DataType = 'Timeseries', this is the POT sample computed
           with the tool; as a matrix array with format:
           col(01): time of POT values in seriel date number format
           col(02): POT values
           col(03): data time range used for selection of POT value: lower bound
           col(04): data time range used for selection of POT value: upper bound
         When DataType = 'POT', this is the same input POT dataset.
     > MRL_output: output of Mean Residual Life (MRL) function; as a structure array with fields:
         Summary = summary of the threshold selection results, as a table array with format:
            Threshold: list of threshold values evaluated
            MeanExcess: mean excess
            Weight: weights
            WMSE: weighted mean square error (MSE)
            GPD_Shape: GPD shape parameter
            GPD_Scale: GPD scale parameter
            Events: number of events above each threshold
            Rate: sample intensity (mean annual rate of events using sample of events above threshold)

         Selection = selected threshold with other parameters, as a table array with format:
            Criterion: criterion applied by MRL method for selecting the threshold
            Threshold: selected threshold value
            id_Summary: location (row ID) of selected threshold in the Summary field
            Events: number of events above the selected threshold
            Rate: sample intensity (mean annual rate of events using sample of events above threshold)

         pd_TH_wOut = GPD threshold parameter values used in the bootstrap process
         pd_k_wOut = initial GPD shape parameter values obtained in the bootstrap process
         pd_sigma = GPD scale parameter values used in the bootstrap process
         pd_k_mod = modified GPD shape parameter values used in the bootstrap process
         eMsg = Status message indicating when the MRL methodology was not able
            to objectively determine a theshold for the GPD. In this case, the
            GPD threshold is set to 0.99 times the minimum value of the bootstrap sample.
            Otherwise, eMsg will be empty.

     > HC_plt: full HC; as a structure array with two fields: 'out' and 'MRL_Crit'. It may contain
         up to two records (one per MRL GPD threshold) with the following format:

         out = numerical array with the following 5 rows
            row(01): mean values
            row(02): values of 2% percentile or 1st percentage of input "prc"
            row(03): values of 16% percentile or 2nd percentage of input "prc"
            row(04): values of 84% percentile or 3rd percentage of input "prc"
            row(05): values of 98% percentile or 4th percentage of input "prc"

         MRL_Crit = character vector indicating the MRL threshold
            selection criterion used for the HC_plt.out record.

     > HC_tbl: summarized HC; as a structure array with two fields: 'out' and 'MRL_Crit'. It may contain
         up to two records (one per MRL GPD threshold) with the following format:

         out = numerical array with the following 5 rows
            row(01): mean values
            row(02): values of 2% percentile or 1st percentage of input "prc"
            row(03): values of 16% percentile or 2nd percentage of input "prc"
            row(04): values of 84% percentile or 3rd percentage of input "prc"
            row(05): values of 98% percentile or 4th percentage of input "prc"

         MRL_Crit = character vector indicating the MRL threshold
            selection criterion used for the HC_tbl.out record.

     > HC_tbl_rsp_x: hazard values interpolated from the HC, that correspond
         to the responses in HC_tbl_rsp_y; as a structure array with two fields: 'out' and 'MRL_Crit'. It may contain up
         to two records (one per MRL GPD threshold) with the following format:

         out = numerical array with the following 5 rows
            row(01): mean values
            row(02): values of 2% percentile or 1st percentage of input "prc"
            row(03): values of 16% percentile or 2nd percentage of input "prc"
            row(04): values of 84% percentile or 3rd percentage of input "prc"
            row(05): values of 98% percentile or 4th percentage of input "prc"

         MRL_Crit = character vector indicating the MRL threshold
            selection criterion used for the HC_tbl_rsp_x.out record.

     > HC_emp: empirical HC; as a table array with column headings as follows:
         Response: Response vector sorted in descending order
         Rank: Rank or Weibull's plotting position
         CCDF: Complementary cumulative distribution function (CCDF)
         Hazard: hazard as AEF or AEP
         ARI: annual recurrence interval (ARI)

     > ME = error message when the tool failed to evaluate a dataset.

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

CONTRIBUTORS:
    ERDC-CHL Coastal Hazards Group

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: Alpha v0.1: Updated documentation.
20201026-ERS: Added capability to include skew surge. Updated documentation and inputs.
20201031-ERS: Alpha v0.2: Reformated logic to reduce execution time. Added
    option to execute the tool in Regular mode or Fast mode.
20201202-ERS: Alpha v0.2: Reviewed logic, added extra layer to run simple
    mode w/o PCT, updated documentation.
20201203-ERS: Alpha v0.2: Updated documentation.
20201208-ERS: Alpha v0.2: Finalized draft documentation for this version.
20201218-ERS: Alpha v0.2: Corrected the input order in the plot function
    call when hasPCT=0 and GPD_TH_crit~=0. Also moved the checkpoint for TH_otlr
    to verify first if it's empty.
20201221-ERS: Alpha v0.2: Second input dataset is now being pre-processed
    for removal of invalid POT values.
20210113-ERS: Alpha v0.2: Changed hasPCT.m to assume pool is
    inactive when PCT is not found. In addition, when ind_Skew = 1, the
    second input dataset will be pre-processed.
20210303-ERS: alpha version 0.2: minor corrections in StormSim_SST_Fit.m.
20210311-ERS: alpha version 0.3: corrections in StormSim_SST_Fit.m,
    StormSim_SST_FitSimple.m and StormSim_SST_FitSimplePar.m: adjustment of
    hazards to be monotonic; creation of hazard tables.
20210325-ERS: alpha v0.4: modified outputs and plots to account for new predefined AEF table values. Also
    organized MRL script and output.
20210405-ERS: alpha v0.4: updated documentation
20210406-ERS: alpha v0.4: modified preprocessing of inputs datasets.
    Modified the MRL, the fit and plotting scripts to account for the
    computation of the default GPD threshold.
20210407-ERS: alpha v0.4: Modified the fit functions to account for small
    bootstrap samples.
20210412-ERS: alpha v0.4: Modified the fit functions to discard bootstrap
    samples with spurious values.
20210419-ERS: alpha v0.4: corrected evaluation of input TH_otlr.
20210420-ERS: alpha v0.4: expanded use of inputs yaxis_Label and yaxis_Limits
    (see description). Took out HC_tbl_rsp_y as another user input.
20210421-ERS: alpha v0.5: corrections for yaxis_Label and yaxis_Limits.
20210503-ERS: alpha v0.5: corrections applied to StormSim_MRL.m. Limits applied
    for the GPD shape parameter. Option provided to evaluate small POT samples
    with either GPD or empirical distribution.
20210505-ERS: alpha v0.5: the tool now will not apply skew tides when
    ind_Skew = 1 but the GPR model and/or 2nd dataset (with tides) have
    invalid formats. Added output Check_datasets to display a message and
    the list of stations, so the user can double check them. Check points
    expanded in all 3 fit scripts.
20210507-ERS: alpha v0.5: added input to modify the GPD shape parameter.
20210511-ERS: alpha v0.5: removed the options to modify the GPD shape
    param since a unique range is now enforced by default. TH_otlr was
    removed. Removed LICENSING and DISCLAIMER comments from internal
    functions.
20210512-ERS: alpha v0.5: re-enabled the patch of 20210412 to ignore bad
    random samples. Also had a succesful run while testing with time series data for NRC Pilot Study.
20210513-ERS: alpha v0.5: cleaned script and updated documentation.
20210517-ERS: alpha v0.5: corrected MRL script to ignore, error due to Inf
    weights. Found while doing Nantucket.
20210520-ERS: alpha v0.5: corrected error in POT script
20210525-ERS: alpha v0.5: corrected bug in preprocessing
20210526-ERS: alpha v0.5: now changing to NaN all values < AEF=10^-4 in the
    HC table summary for HC_tbl_rsp_y.
20210527-ERS: alpha v0.5: minor correction
20210601-ERS: alpha v0.5: added input to control activation of PCT.
20210713-ERS: alpha v0.5: to ensure a GPD threshold is found, reduced the kernel smoothing bandwidth in MRL script from 1/4 to 1/7.
20210719-ERS: alpha v0.5: added warning message when mean HC > 1.75*emp HC.
20210727-ERS: alpha v0.5: minor correction.
20210809-ERS: alpha v0.5: revised.


***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function [HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,Removed_datasets,Check_datasets,SST_output] = StormSim_SST_Tool_R1_v20210809(input_data,flag_value,tLag,lambda,Nyrs,path_out,staID,yaxis_Label,yaxis_Limits,prc,use_AEP,GPD_TH_crit,SLC,ind_Skew,gprMdl,DataType,ExecMode,HC_tbl_rsp_y,apply_GPD_to_SS,apply_Parallel)
%% Other settings
clc;
disp(['***********************************************************' newline...
    '***         StormSim-SST Tool Alpha Version 0.5         ***' newline...
    '***                Release 1 - 20210809                 ***' newline...
    '***                 FOR  TESTING  ONLY                  ***' newline...
    '***********************************************************' newline...
    '**********           NOTES TO THE USER           **********' newline...
    '*** 1) Refer to Quick Start Guide for a description     ***' newline...
    '***    of settings, inputs and outputs.                 ***' newline...
    '*** 2) Please report any error using the Feedback Form. ***' newline...
    '*** 3) For questions or help contact:                   ***' newline...
    '***    Efrain.Ramos-Santiago@usace.army.mil             ***' newline...
    '***********************************************************'])

disp([newline '*** Step 1: Processing input arguments '])

% Turn off all warnings
warning('off','all');

% Check apply_GPD_to_SS
if sum(apply_GPD_to_SS==[0 1])~=1||isempty(apply_GPD_to_SS)||isnan(apply_GPD_to_SS)
    error('Input apply_GPD_to_SS must be 0 or 1')
end

% Check use_AEP
if sum(use_AEP==[0 1])~=1||isempty(use_AEP)||isnan(use_AEP)
    error('Input use_AEP must be 0 or 1')
end

% Set up probabilities for HC summary
if use_AEP %Select AEPs
    HC_tbl_x = 1./[2 5 10 20 50 100 200 500 1000 2000 5000 1e4 2e4 5e4 1e5 2e5 5e5 1e6];
else %Select AEFs
    %     HC_tbl_x = 1./[0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000];
    HC_tbl_x = 1./[0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000];
end

% Set up responses for HC summary
if isempty(HC_tbl_rsp_y)
    HC_tbl_rsp_y =(0.01:0.01:20)';
elseif isrow(HC_tbl_rsp_y)
    HC_tbl_rsp_y =sort(HC_tbl_rsp_y)';
    if sum(HC_tbl_rsp_y<=0)~=0
        error('Input HC_tbl_rsp_y must be specified as a numerical vector of positive values')
    end
end

% Set up AEFs for full HC; in log10 scale (for plotting), from 10^1 to 10^-6
d=1/90; v=10.^(1:-d:0)'; HC_plt_x=v; x=10;
for i=1:6, HC_plt_x=[HC_plt_x; v(2:end)/x]; x=x*10; end %#ok<AGROW>
HC_plt_x=flipud(HC_plt_x);

% Check path to output folder
if isempty(path_out)
    path_out='SST_output';mkdir(path_out);
elseif ~exist(path_out,'dir')
    mkdir(path_out);
end
path_out=['.\',path_out,'\'];

% Check size and value of input arguments
sz = length(input_data);
if isempty(input_data)||~isstruct(input_data)
    error('Input input_data cannot be an empty or matrix array. Must be a structure array.')
end
if isempty(staID)||size(staID,1)~=sz||~iscellstr(staID)||isstring(staID)
    error('Input staID cannot be empty. Provide one ID number per dataset.')
end

% Evaluate parameters for POT or Time series input data
DataType = find(strcmp(DataType,{'Timeseries','POT'}));
switch DataType
    case 1 %Time series
        if sum(tLag<=0)~=0||isempty(tLag)
            error('When DataType = ''Timeseries'': Input tLag cannot be an empty array and must have positive values only')
        end
        if and(length(tLag)==1,sz>1)
            tLag=repmat(tLag,1,sz);
        end
        if sum(lambda<=0)~=0||isempty(lambda)
            error('When DataType is set to ''Timeseries'': Input lambda cannot be an empty array and must have positive values only')
        end
        if length(lambda)==1 && sz>1
            lambda=repmat(lambda,1,sz);
        end
        Nyrs=NaN(1,sz);
    case 2 %POT
        if sum(Nyrs<=0)~=0||isempty(Nyrs)
            error('When DataType is set to ''POT'': Input "Nyrs" cannot be empty and must have positive values only')
        end
        if length(Nyrs)==1 && sz>1
            Nyrs=repmat(Nyrs,1,sz);
        end
    otherwise
        error('Unrecognized value for DataType. Available options are ''Timeseries'' or ''POT'' ')
end

% Check if Parallel Computing Toolbox exists, start pool session and Turn off warnings
[id_PCT,id_act]=hasPCT;
if GPD_TH_crit~=0 && id_PCT && apply_Parallel
    if id_act
        parpool;
    end
    parfevalOnAll(gcp,@warning,0,'off','all');
end

% Other settings
if sum(GPD_TH_crit<0||GPD_TH_crit>2)~=0
    error('Input GPD_TH_crit has wrong values. Available options are 0, 1, or 2')
end

if strcmp(ExecMode,'Regular')
    if sum(strcmp(yaxis_Label,''))~=0||isempty(yaxis_Label)||~iscellstr(yaxis_Label)||isstring(yaxis_Label)
        error('Error found in input yaxis_Label. Refer to the Quick Start Guide for instructions on how to set up this input.');
    end
    if length(yaxis_Label)~=sz
        yaxis_Label = repmat(yaxis_Label,sz,1);
    end
else
    yaxis_Label=cell(sz,1);
end

if strcmp(ExecMode,'Regular')
    if ~isempty(yaxis_Limits)
        if isvector(yaxis_Limits)
            n=length(yaxis_Limits);
            if n~=2
                error(['When input yaxis_Limits is entered as a numeical vector, it must' newline...
                    'have two values.'])
            end
            yaxis_Limits=sort(yaxis_Limits);
            if iscolumn(yaxis_Limits)
                yaxis_Limits=yaxis_Limits';
            end
            yaxis_Limits = repmat(yaxis_Limits,sz,1);
        elseif ismatrix(yaxis_Limits)
            [m,n]=size(yaxis_Limits);
            if or(m~=sz,n~=2)
                error(['When input yaxis_Limits is entered as a numerical matrix, it must' newline...
                    'have two columns and one row per input dataset in input_data'])
            end
        else
            error('Error found in input yaxis_Limits. Refer to the Quick Start Guide for instructions on how to set up this input.');
        end
    else
        yaxis_Limits=double.empty(sz,0);
    end
else
    yaxis_Limits=double.empty(sz,0);
end

if isempty(prc)
    prc=[2.28 15.87 84.13 97.72]';
else
    if length(prc)>4||sum(isnan(prc))>0||sum(isinf(prc))>0||sum(prc<0)~=0
        error('Input prc can have 1 to 4 percentages in the interval [0,100].');
    end
    prc=sort(prc,'ascend'); prc=prc(:);
end
if isempty(SLC)
    SLC=0;
else
    if length(SLC)~=1||SLC<0||isnan(SLC)||isinf(SLC)
        error('Input SLC must be a positive scalar')
    end
end
if sum(ind_Skew==[0 1])~=1||isempty(ind_Skew)||isnan(ind_Skew)||isinf(ind_Skew)
    error('Input ind_Skew must be 0 or 1')
end
if ind_Skew==1
    if isempty(gprMdl)
        error('When ind_Skew = 1: Input gprMdl.mdl cannot be an empty structure')
    end
    if length(gprMdl)~=sz
        error(['Input gprMdl must be same size as input input_data:' newline 'Provide one GPR model per dataset included in input_data'])
    end
else
    for i=1:sz
        gprMdl(i).mdl=[];
    end
end

% Compatibility check
a=version('-release'); if a(end)=='a',a(end)='0';else; a(end)='1';end
a = str2double(a)<20200;

% Types of SST execution
ExecMode = find(strcmp(ExecMode,{'Regular','Fast'}));

% Pre-allocate output structure array
switch ExecMode
    case 1 %Regular
        SST_output = struct('staID','','RL',double.empty,'POT',double.empty,'MRL_output',double.empty,...
            'HC_plt',double.empty,'HC_tbl',double.empty,'HC_tbl_rsp_x',double.empty,'HC_emp',double.empty,'Warning','','ME',cell(1));
    case 2 %Fast
        SST_output = struct('staID','','HC_plt',double.empty,'HC_tbl',double.empty,'HC_tbl_rsp_x',double.empty,'Warning','','ME',cell(1));
    otherwise
        error('Unrecognized value for ExecMode. Available options are ''Regular'' or ''Fast''')
end

disp('*** Step 2: Verifying input datasets')

%% Data preprocessing: remove NaN, inf values; compute record length
procData(sz).dat = [];
procData(sz).dat2 = [];
procData(sz).id = [];
procData(sz).id2 = [];
sz = 1:sz;

switch DataType
    case 1 %Timeseries
        for i=sz %dataset loop
            data_time = input_data(i).time_values;
            data_values = input_data(i).data_values;
            
            % Remove flag values?
            if ~isempty(flag_value)
                data_time(data_values==flag_value)=[];
                data_values(data_values==flag_value)=[];
            end
            
            % Merge data and remove NaN, Inf values
            data_values=[data_time(:) data_values(:)];
            data_values(isinf(data_values(:,2))|isnan(data_values(:,2)),:)=[];
            
            procData(i).dat = data_values;
            
            % Compute record length: Effective duration method
            dt=[]; [dt(:,1),dt(:,2),dt(:,3)]=ymd(datetime(data_time,'ConvertFrom','datenum'));
            dt=unique(dt,'rows'); Nyrs(i)=size(dt,1)/365.25;
            
            % Can dataset be evaluated? Mark empty datasets to delete later
            if isempty(data_values)
                procData(i).id = i;
            end
        end
        
    case 2 %POT
        for i=sz %dataset loop
            
            %take datasets
            tt = input_data(i).time_values;
            data_values = input_data(i).data_values;
            
            if isempty(tt)
                tt = NaN(length(data_values),1);
            end
            
            if ind_Skew
                data_values2 = input_data(i).data_values2;
            end
            
            % Remove flag values?
            if ~isempty(flag_value)
                id_1 = data_values==flag_value;
                data_values(id_1)=[];
                tt(id_1)=[];
                if ind_Skew && ~isempty(data_values)
                    data_values2(id_1)=[];
                    data_values(data_values2==flag_value)=[];
                    tt(data_values2==flag_value)=[];
                    data_values2(data_values2==flag_value)=[];
                end
            end
            
            % Remove NaN, Inf and nonpositive values
            data_values = data_values(:);
            id_1 = isinf(data_values)|isnan(data_values)|data_values<=0;
            data_values(id_1,:)=[];
            tt(id_1,:)=[];
            
            if ind_Skew && ~isempty(data_values) %&& ~isempty(data_values2)
                data_values2 = data_values2(:);
                data_values2(id_1)=[];
                id_1 = isinf(data_values2)|isnan(data_values2)|data_values2<=0;
                data_values(id_1)=[];
                data_values2(id_1)=[];
                tt(id_1)=[];
            end
            
            % Sample only has the same value?
            [~,ia,~] = unique(data_values,'stable');
            ia = length(ia)<=3;
            
            % Merge data and store
            data_values = [tt data_values]; %#ok<AGROW>
            procData(i).dat = data_values;
            
            if ind_Skew && ~isempty(data_values) %&& ~isempty(data_values2)
                procData(i).dat2 = data_values2;
            end
            
            % Can dataset be evaluated? Mark empty datasets to delete later
            if isempty(data_values) || ia
                procData(i).id = i;
            end
            
            % If skew tides are needed, is there a valis model provided?
            if ind_Skew && (isempty(data_values2) || ~isobject(gprMdl(i).mdl))
                procData(i).id2 = i;
            else
                procData(i).id2 = 0;
            end
        end
end


%% Select dataset to evaluate
Removed_datasets = '';
id = [procData.id];
if ~isempty(id)
    sz(id)=[];
    Removed_datasets = {'Cannot evaluate these stations for any of the following reasons:';...
        '- dataset was empty after removal of NaN/Inf/flag values (and nonpositive values when input dataset is a POT)';...
        '- dataset consists of 3 or less unique values repeated many times'};
    
    for i=id
        Removed_datasets = [Removed_datasets; staID{i}]; %#ok<AGROW>
    end
end
j=0;


%% Store message that ind_Skew = 1 but invalid POT_samp2 and/or GPR model object
Check_datasets = '';
if ind_Skew
    id2 = [procData.id2];
    if ~isempty(id)
        id2(id)=[];
    end
    id2(id2==0)=[];
    if ~isempty(id2)
        Check_datasets = {'These stations were evaluated without applying skew tides since:';...
            '- ind_Skew = 1 but entered invalid values for either the second POT dataset (with tides) or the GPR model'};
        
        for i=id2
            Check_datasets = [Check_datasets; staID{i}]; %#ok<AGROW>
        end
    end
end


%% Perform SST
switch DataType
    case 1 %Timeseries
        if GPD_TH_crit==0 % Evaluate all MRL thresholds
            for i=sz %dataset loop
                disp(['*** Step 3: Performing SST for station ',staID{i,1}])
                j=j+1;
                try
                    % Execute StormSim-POT
                    [POT_samp,~] = StormSim_POT(procData(i).dat(:,1),procData(i).dat(:,2),tLag(i),lambda(i),Nyrs(i));
                    
                    % Execute StormSim-SST-Fit
                    [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output] = StormSim_SST_Fit(POT_samp(:,2),Nyrs(i),HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,procData(i).dat2,SLC,gprMdl(i).mdl,staID(i,1),yaxis_Label{i},path_out,yaxis_Limits(i,:),apply_GPD_to_SS);
                    
                    % Gather the output
                    switch ExecMode
                        case 1 %Regular
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).RL = Nyrs(i);
                            SST_output(j).POT = POT_samp;
                            SST_output(j).MRL_output = MRL_output;
                            SST_output(j).HC_plt = HC_plt;
                            SST_output(j).HC_tbl = HC_tbl;
                            SST_output(j).HC_tbl_rsp_x = HC_tbl_rsp_x;
                            SST_output(j).HC_emp = HC_emp;
                            
                            % Plot results
                            StormSim_SST_Plot(HC_plt,HC_emp,MRL_output,prc,staID(i,:),yaxis_Label{i},path_out,yaxis_Limits(i,:),use_AEP,GPD_TH_crit,a,HC_plt_x2)
                            
                        case 2 %Fast
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                    end
                catch ME
                    SST_output(j).staID = staID{i,1};
                    SST_output(j).ME = ME;
                end
            end
            
        elseif id_PCT && GPD_TH_crit~=0 && apply_Parallel % PCT available and evaluate only one MRL threshold
            for i=sz %dataset loop
                disp(['*** Step 3: Performing SST for station ',staID{i,1}])
                j=j+1;
                try
                    % Execute StormSim-POT
                    [POT_samp,~] = StormSim_POT(procData(i).dat(:,1),procData(i).dat(:,2),tLag(i),lambda(i),Nyrs(i));
                    
                    % Execute StormSim-SST-Fit-Simple
                    [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output,str1] = StormSim_SST_Fit_SimplePar(POT_samp(:,2),Nyrs(i),HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,procData(i).dat2,SLC,gprMdl(i).mdl,apply_GPD_to_SS);
                    
                    % Gather the output
                    switch ExecMode
                        case 1 %Regular
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).RL = Nyrs(i);
                            SST_output(j).POT=POT_samp;
                            SST_output(j).MRL_output=MRL_output;
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).HC_emp=HC_emp;
                            SST_output(j).Warning=str1;
                            
                            % Plot results
                            StormSim_SST_Plot_Simple(HC_emp,HC_plt,HC_plt_x2,MRL_output,prc,use_AEP,staID(i,:),yaxis_Label{i},path_out,yaxis_Limits(i,:),a,GPD_TH_crit)
                            
                        case 2 %Fast
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).Warning=str1;
                    end
                catch ME
                    SST_output(j).staID = staID{i,1};
                    SST_output(j).ME = ME;
                end
            end
            
        else % PCT not available and evaluate only one MRL threshold
            for i=sz %dataset loop
                disp(['*** Step 3: Performing SST for station ',staID{i,1}])
                j=j+1;
                try
                    % Execute StormSim-POT
                    [POT_samp,~] = StormSim_POT(procData(i).dat(:,1),procData(i).dat(:,2),tLag(i),lambda(i),Nyrs(i));
                    
                    % Execute StormSim-SST-Fit-Simple
                    [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output,str1] = StormSim_SST_Fit_Simple(POT_samp(:,2),Nyrs(i),HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,procData(i).dat2,SLC,gprMdl(i).mdl,apply_GPD_to_SS);
                    
                    % Gather the output
                    switch ExecMode
                        case 1 %Regular
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).RL = Nyrs(i);
                            SST_output(j).POT=POT_samp;
                            SST_output(j).MRL_output=MRL_output;
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).HC_emp=HC_emp;
                            SST_output(j).Warning=str1;
                            
                            % Plot results
                            StormSim_SST_Plot_Simple(HC_emp,HC_plt,HC_plt_x2,MRL_output,prc,use_AEP,staID(i,:),yaxis_Label{i},path_out,yaxis_Limits(i,:),a,GPD_TH_crit)
                            
                        case 2 %Fast
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).Warning=str1;
                    end
                catch ME
                    SST_output(j).staID = staID{i,1};
                    SST_output(j).ME = ME;
                end
            end
        end
        
    case 2 %POT
        if GPD_TH_crit==0
            for i=sz %dataset loop
                disp(['*** Step 3: Performing SST for station ',staID{i,1}])
                j=j+1;
                try
                    % Execute StormSim-SST-Fit
                    [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output] = StormSim_SST_Fit(procData(i).dat(:,2),Nyrs(i),HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,procData(i).dat2,SLC,gprMdl(i).mdl,staID(i,1),yaxis_Label{i},path_out,yaxis_Limits(i,:),apply_GPD_to_SS);
                    
                    % Gather the output
                    switch ExecMode
                        case 1 %Regular
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).RL = Nyrs(i);
                            SST_output(j).POT=procData(i).dat;
                            SST_output(j).MRL_output=MRL_output;
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).HC_emp=HC_emp;
                            
                            % Plot results
                            StormSim_SST_Plot(HC_plt,HC_emp,MRL_output,prc,staID(i,:),yaxis_Label{i},path_out,yaxis_Limits(i,:),use_AEP,GPD_TH_crit,a,HC_plt_x2)
                            
                        case 2 %Fast
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                    end
                catch ME
                    SST_output(j).staID = staID{i,1};
                    SST_output(j).ME = ME;
                end
            end
            
        elseif id_PCT && GPD_TH_crit~=0 && apply_Parallel % PCT available and evaluate only one MRL threshold
            for i=sz %dataset loop
                disp(['*** Step 3: Performing SST for station ',staID{i,1}])
                j=j+1;
                try
                    % Execute StormSim-SST-Fit-Simple
                    [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output,str1] = StormSim_SST_Fit_SimplePar(procData(i).dat(:,2),Nyrs(i),HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,procData(i).dat2,SLC,gprMdl(i).mdl,apply_GPD_to_SS);
                    
                    % Gather the output
                    switch ExecMode
                        case 1 %Regular
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).RL = Nyrs(i);
                            SST_output(j).POT=procData(i).dat;
                            SST_output(j).MRL_output=MRL_output;
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).HC_emp=HC_emp;
                            SST_output(j).Warning=str1;
                            
                            % Plot results
                            StormSim_SST_Plot_Simple(HC_emp,HC_plt,HC_plt_x2,MRL_output,prc,use_AEP,staID(i,:),yaxis_Label{i},path_out,yaxis_Limits(i,:),a,GPD_TH_crit)
                            
                        case 2 %Fast
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).Warning=str1;
                    end
                catch ME
                    SST_output(j).staID = staID{i,1};
                    SST_output(j).ME = ME;
                end
            end
            
        else % PCT not available and evaluate only one MRL threshold
            for i=sz %dataset loop
                disp(['*** Step 3: Performing SST for station ',staID{i,1}])
                j=j+1;
                try
                    % Execute StormSim-SST-Fit-Simple
                    [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output,str1] = StormSim_SST_Fit_Simple(procData(i).dat(:,2),Nyrs(i),HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,procData(i).dat2,SLC,gprMdl(i).mdl,apply_GPD_to_SS);
                    
                    % Gather the output
                    switch ExecMode
                        case 1 %Regular
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).RL = Nyrs(i);
                            SST_output(j).POT=procData(i).dat;
                            SST_output(j).MRL_output=MRL_output;
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).HC_emp=HC_emp;
                            SST_output(j).Warning=str1;
                            
                            % Plot results
                            StormSim_SST_Plot_Simple(HC_emp,HC_plt,HC_plt_x2,MRL_output,prc,use_AEP,staID(i,:),yaxis_Label{i},path_out,yaxis_Limits(i,:),a,GPD_TH_crit)
                            
                        case 2 %Fast
                            SST_output(j).staID = staID{i,1};
                            SST_output(j).HC_plt=HC_plt;
                            SST_output(j).HC_tbl=HC_tbl;
                            SST_output(j).HC_tbl_rsp_x=HC_tbl_rsp_x;
                            SST_output(j).Warning=str1;
                    end
                catch ME
                    SST_output(j).staID = staID{i,1};
                    SST_output(j).ME = ME;
                end
            end
        end
end

if use_AEP %Convert to AEP
    HC_plt_x = aef2aep(HC_plt_x);
end

%% Save the output
disp(['*** Step 4: Saving results here: ',path_out])
save([path_out,'StormSim_SST_output.mat'],'SST_output','HC_tbl_x','HC_plt_x','HC_tbl_rsp_y','Removed_datasets','Check_datasets','-v7.3')

[~,id_act]=hasPCT;if id_act==0,delete(gcp);end
disp('*** Evaluation finished.' )
disp(['****** Remember to check outputs Check_datasets and Removed_datasets.' newline])
disp('*** StormSim-SST Tool terminated.')

end

%%
function [in] = aef2aep(in)
in = (exp(in)-1)./exp(in);
end

%% ecdf_boot.m
%{
DESCRIPTION:
   This script performs the bootstrap simulation.

INPUT ARGUMENTS:
  - empHC = empirical HC; as a matrix with format:
     col(01): Response vector sorted in descending order
     col(02): Rank or Weibull's plotting position
     col(03): Complementary cumulative distribution function (CCDF)
     col(04): hazard as AEF or AEP
     col(05): annual recurrence interval (ARI)
  - Nsim: number of simulations or bootstrap iterations
 
OUTPUT ARGUMENTS:
  - boot: results from bootstrap process

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)

CONTRIBUTORS:
    Efrain Ramos-Santiago (ERS)

HISTORY OF REVISIONS:
20201015-ERS: revised. Updated documentation.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function boot = ecdf_boot(empHC,Nsim)
%     rng('default');
Nstrm = size(empHC,1);
dlt = abs(diff(empHC)); dlt = [dlt;dlt(end)];
boot=zeros(Nsim,Nstrm);
for i = 1:Nsim
    [x,idx] = datasample(empHC,Nstrm);
    y = x + randn(Nstrm,1).*dlt(idx);
    boot(i,:) = sort(y,'descend')';
end
end

%% hasPCT.m
%{
By: E. Ramos
Description: Function to determine if MATLAB has the Parallel Computing
   Toolbox (PCT) and if a pool is active.
   Inputs: None.
   Outputs: id_PCT = 1 when the PCT is found.
            id_PCT = 0 otherwise.

            id_act = 1 when the parallel pool is not active.
            id_act = 0 otherwise

History of revisions:
20201201-ERS: created.
20210113-ERS: now initially assuming id_act=1, to avoid evaluating gcp()
    when id_PCT=0.
%}
function [id_PCT,id_act]=hasPCT()
ver2=ver;
ver2={ver2.Name};
id_PCT=sum(ver2=="Parallel Computing Toolbox");
id_act=1;%initially assume pool is inactive
if id_PCT
    id_act=isempty(gcp('nocreate'));
end
end

%% KernReg_LocalMean.m
%{
Developed by: E. Ramos-Santiago

Description:
   This function is an adaptation from the original script, to only apply the
    local mean or Nadaraya-Watson estimator (non-parametric regression model) to uni- or multi-variate data,
   following Wand & Jones (1994; Ch5, Sec2). Multivariate regression
   possible to parameters with up to 3 predictors, with order p set to 0
   or 1, only. Normal kernel function is used. Same kernel and bandwidth
   applied to all predictors.

Input:
     x = predictor matrix (one predictor per column)
     y = response parameter (column vector)
     h = bandwidth

Output:
     Y = model prediction (column vector)
     X = predictor grid (same as x)

References:
  Wand, M. P., and Jones, M. C. (1994). Kernel smoothing. Chapman &
    Hall/CRC Monographs on Statistics & Applied Probability, Chapman
    and Hall/CRC press.

History of revisions:
20171206-ERS: developed first draft of function.
20180228-ERS: made revision.
20190109-ERS: expanded algorithm to incorporate multivariate regression using p=0/p=1, only.
20210505-ERS: adapted.
%}
function [X,Y] = KernReg_LocalMean(x,y,H)
[n,k] = size(x); %sample size and dimensionality
X = x;

% setup regression parameters
NN=size(X,1);Y=zeros(NN,1);
Kh =@(H,t,k)H^(-k)*(2*pi)^(k/2)*exp(-0.5*t); %d-variate Normal kernel

%do regression
for i = 1:NN %grid point loop
    t = sum(bsxfun(@minus, X(i,:), x).^2,2)./(H^2);
    Xx = ones(n,1);
    Wx = diag(Kh(H,t,k));
    Y(i) = (Xx'*Wx*Xx)\(Xx'*Wx)*y;
end
end

%% Monotonic_adjustment.m
%{
By: E. Ramos-Santiago
Description: Script to adjust the trend of hazard curves with jumps. The
  StormSim-SST tool can produce non-monotonic curves when the GPD threshold
  parameter returned by the MRL selection method is too low. This causes
  incomplete random samples and the potential to have jumps in the mean
  curve and CLs.
History of revisions:
20210310-ERS: created function to patch the SST tool.
%}
function [y_in] = Monotonic_adjustment(x_in,y_in)

% Change inputs to column vector
y_in=y_in(:);x_in=log(x_in(:));

% Take numbers
idx3=~isnan(y_in);y=y_in(idx3);x=x_in(idx3);

% Compute slope
dx=diff(x);dy=diff(y);s=dy./dx;

% Identify positive slopes
id_pos=find(s>0);

% Adjust y values with positive slopes
if ~isempty(id_pos)
    for i=1:length(id_pos)
        y(id_pos(i)+1)=mean(s(id_pos(i)-2:id_pos(i)-1))*dx(id_pos(i))+y(id_pos(i));
    end
    
    % Return values to original position
    y_in(idx3)=y;
end
y_in=y_in';
end

%% StormSim_MRL.m
%{
SOFTWARE NAME:
    StormSim-SST-Fit (Statistics)

DESCRIPTION:
   This script applies the mean residual life (MRL) methodology to
   objectively select the parameters of the Generalized Pareto Distribution
function.

INPUT ARGUMENTS:
  - POT_samp: empirical distribution of the Peaks-Over-Threshold sample, as
      computed in StormSim_SST_Fit.m
  - Nyrs: record length in years of the input time series data set;
      specified as a positive scalar
     
AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: revised. Updated documentation.
20201215-ERS: added break in for loop to avoid evaluating thresholds
    returning excesses of same magnitude.
20210325-ERS: organized the threshold computation; the 3rd threshold not included
    in the output anymore; output organized into table arrays stored in a structure array.
20210406-ERS: identified error when input sample size <10. When this
    happens, an empty array is returned and the Default GPD threshold computed
    in the three fit scripts.
20210412-ERS: now selecting a minima when no inflexion point exists; avoiding the script to crash.
20210429-ERS: now removing noise from the min WMSE through kernel
    regression. Also corrected the minimum WMSE criterion based on Langousis. Removed patch
    applied on 20210412.
20210430-ERS: script will stop and return empty arrays when no minima is
    found by WRMS criterion.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function [MRL_output] = StormSim_MRL(PEAKS,Nyrs)

%% Step 1: Use sorted POT dataset as initial values of threshold parameter.
th = flipud(PEAKS); %threshold parameter; sorted in ascending order
N = size(th,1);

% preallocate output
mrl=NaN(1,8);mrl = array2table(mrl,'VariableNames',{'Threshold','MeanExcess','Weight','WMSE','GPD_Shape','GPD_Scale','Events','Rate'});
MRL_output.Summary=mrl;
mrl=table;mrl.Criterion(1)={''};mrl.Threshold(1)=NaN;mrl.id_Summary(1)=NaN;mrl.Events(1)=NaN;mrl.Rate(1)=NaN;
MRL_output.Selection=mrl;

if N>20
    
    %% Step 2: Estimate mean values of excesses
    mrl=NaN(N-10,8); %pre-allocation for speed
    % col(01): thresholds, u
    % col(02): mean excesses, e(u)
    % col(03): weights
    % col(04): weighted mean square errors
    % col(05): GPD shape parameter
    % col(06): GPD scale parameter
    % col(07): Number of events above threshold
    % col(08): Sample intensity or annual rate of events
    
    % Note: Looping up to N-10 as max threshold ensures a min of 10 excesses
    %   e(u) to calculate the conditional mean e(u).
    
    for i = 1:N-10 %threshold loop
        mrl(i,1) = th(i);                %Store threshold
        u = PEAKS(PEAKS(:,1)>th(i),1);   %Take sample of values above threshold
        mrl(i,2) = mean(u-th(i),'omitnan');     %Compute mean excess: difference between sample and threshold
        mrl(i,3) = (N-i)/var(u-th(i),'omitnan');%Compute weights: difference in rank normalized by variance of e(u)
    end
    
    % Note: Weights are used to acount for the increase in the estimation
    %   variance of e(u) with increasing threshold u. Under the assumption of
    %   independence of e(u).
    
    id = isinf(mrl(:,3))|isnan(mrl(:,3));%|(mrl(:,3)<=0);
    %     mrl(id,:) = [];
    
    x = mrl(:,1); %threshold, u
    y = mrl(:,2); %mean excess, e(u)
    w = mrl(:,3); %weights
    
    %% Step 3: Fit a linear model to (u,e(u)) using the weighted least squares
    % method and compute parameters of interest
    
    % Note: Looping up to N-20 as max threshold ensures a min of 10 conditional
    %   means e(u) for the linear regression.
    
    for j = 1:N-20
        
        %Fit GPD and take shape/scale parameters
        u = PEAKS(PEAKS(:,1)>th(j),1);
        
        try
            %Do linear regression
            mdl = fitlm(x(j:end),y(j:end),'Weights',w(j:end),'Exclude',id(j:end));
            %Note: if an error arises with the lack of data, enable the 'Exclude' option.
            
            %Compute weighted mean square error (WMSR) using fit residuals
            mrl(j,4) = mean(w(j:end).*(mdl.Residuals.Raw).^2,'omitnan');
            
            %Fit GPD and take parameters
            pd = fitdist(u,'gp','theta',th(j));
            mrl(j,5) = pd.k; %shape
            mrl(j,6) = pd.sigma;%scale
        catch
        end
    end
    
    %Estimate sample intensity for each threshold considered
    for k = 1:size(mrl,1)
        mrl(k,7) = sum(mrl(:,1)>mrl(k,1)); %no. of events
    end
    mrl(:,8)=mrl(:,7)./Nyrs; %lambda
    
    
    %% Step 4: Threshold selection: Minimum weighted MSE criterion
    mrl(isinf(mrl(:,4))|isnan(mrl(:,4))|mrl(:,4)==0,:)=[]; %delete values with WMSR=0
    
    % remove noise
    [~,~,H] = ksdensity(mrl(:,1)); H=H/7; %no difference found for values < 6
    [~,mrl(:,4)] = KernReg_LocalMean(mrl(:,1),mrl(:,4),H);
    
    %Identify local minima
    TH_id = islocalmin(mrl(:,4),'FlatSelection','first');
    
    % If no min found, return empty arrays
    if sum(TH_id)==0
        %Check: When no inflexion point exists, select the minimum available value
        % [~,TH_id] = min(mrl(:,4),[],'omitnan');
        return;
    end
    TH_id = find(TH_id);
    mrl2=mrl(TH_id,:); %take data of minima
    
    [~,I]=min(mrl2(:,1),[],'omitnan'); %select minimum threshold
    
    % Take results
    TH_selected = table;
    TH_selected.Criterion = {'CritWMSE'};
    TH_selected.Threshold = mrl2(I,1);
    TH_selected.id_Summary = TH_id(I);
    TH_selected.Events = mrl2(I,7);
    TH_selected.Rate = mrl2(I,8);
    
    
    %% Step 5: Threshold selection: Sample intensity (lambda) criterion
    
    % Select the minimum with nearly 2 events/year
    aux = mrl2(:,8)-2;
    aux(aux<-1) = NaN;
    [C,I] = min(abs(aux),[],'omitnan');
    
    if isempty(I)||isnan(C)||isempty(C)
        [~,I] = max(mrl2(:,8),[],'omitnan');
    end
    
    % Take results
    TH_selected.Criterion(2) = {'CritSI'};
    TH_selected.Threshold(2) = mrl2(I,1);
    TH_selected.id_Summary(2) = TH_id(I);
    TH_selected.Events(2) = mrl2(I,7);
    TH_selected.Rate(2) = mrl2(I,8);
    
    %% Store MRL output
    mrl = array2table(mrl,'VariableNames',{'Threshold','MeanExcess','Weight',...
        'WMSE','GPD_Shape','GPD_Scale','Events','Rate'});
    
    MRL_output.Summary=mrl;
    MRL_output.Selection=TH_selected;
end
end

%% StormSim_POT.m
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
    StormSim-POT (Statistics)

DESCRIPTION:
   This script generates the Peaks-Over-Threshold (POT) sample from a raw
   time series data set.

INPUT ARGUMENTS:
    dt = timestamps of response data; as a datenum vector
    Resp = response data; specified as a vector. Cannot be a PDS/POT sample.
    tLag = inter-event time in hours; specified as a scalar
    lambda = mean annual rate of events; specified as a scalar
    Nyrs = record length in years; specified as a scalar

OUTPUT ARGUMENTS:
   POTout = response POT sample; as a matrix with format:
     col(01): timestamps of POT values
     col(02): POT values
     col(03): data time range used for selection of POT value: lower bound
     col(04): data time range used for selection of POT value: upper bound
   Threshold = selected threshold for identification of excesses.

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

CONTRIBUTORS:
    ERDC-CHL Coastal Hazards Group

HISTORY OF REVISIONS:
20200903-ERS: revised.
20200914-ERS: revised.

***************  ALPHA  VERSION  ***************  FOR TESTING  ************
%}
function [POTout,Threshold] = StormSim_POT(dt,Resp,tLag,lambda,Nyrs)
%% A few check points

% Required input
if nargin<5
    error('Missing input arguments: Received less than 5 input arguments.')
elseif nargin>5
    error('Too many input arguments: Function only requires 5 input arguments.')
end

% Reformat some inputs
if isrow(dt),dt=dt';end
if isrow(Resp),Resp=Resp';end

%% Set default and compute additional parameters
thMult = 1; %Default = 1; If threshold multiplier is set too low (e.g., ~0)
% an almost continuous time series will be created in Resp(index), with
% dLag <=48 hr and very few peaks will be identified

Nstm = 100000;
Resp_avg = mean(Resp); %arithmetic mean of response
Resp_std = std(Resp); %standard deviation of response
tLag = tLag/24;  %Convert inter-event time from hours to days

%% Determine the threshold that yields a required amount of peak response events
while Nstm > (Nyrs*lambda)
    
    % Compute threshold and identify response values above it
    Threshold = Resp_avg + Resp_std*thMult;
    index = find(Resp >= Threshold);
    
    % Take corresponding time values
    dt2 = dt(index);
    
    % Compute sample inter-event times
    dt2 = [1;diff(dt2)<=tLag;0];
    
    % Identify inter-event times longer than tLag
    id = find(dt2==0);
    
    % Compute amount of resulting events and increase the multiplier
    Nstm = length(id);
    thMult = thMult + 0.01;
end

%% Classify the response values

% NOTE:
%   stm_col has data indices; Resp_peak has response values; each column
%   has response events inside a time segment defined by the inter-event
%   time

stm_col = NaN(length(dt2),length(id));
stm_col(1:id(1)-1,1) = index(1:id(1)-1);
Resp_peak = NaN(length(dt2),length(id));
Resp_peak(1:id(1)-1,1) = Resp(index(1:id(1)-1));

for k = 2:length(id)
    id2 = id(k-1):id(k)-1;
    stm_col(id2-(id(k-1)-1),k) = index(id2);
    Resp_peak(id2-(id(k-1)-1),k) = Resp(index(id2));
end

%% Extract peak response value for each event
[Pk_val,I] = max(Resp_peak,[],1,'omitnan');

%% Locate events inside arrays
ParINDEX = NaN(length(I),3);
for i=1:Nstm
    ParINDEX(i,1) = stm_col(I(i),i);
end
ParINDEX(:,2) = min(stm_col,[],1,'omitnan')';
ParINDEX(:,3) = max(stm_col,[],1,'omitnan')';

%% Extract timestamps of peak events and store results for output
POTout = dt(ParINDEX(:,1),1); %timestamp of peak event
POTout(:,2) = Pk_val'; %peak event

% Data time range used for selection of POT value:
POTout(:,3) = dt(ParINDEX(:,2),1); %lower bound
POTout(:,4) = dt(ParINDEX(:,3),1); %upper bound

%% Remove NaN / Inf / negative values
POTout(isnan(POTout(:,2))|isinf(POTout(:,2))|POTout(:,2)<=0,:)=[];

end

%% StormSim_SST_Fit.m
%{
SOFTWARE NAME:
    StormSim-SST-Fit (Statistics)

DESCRIPTION:
   This script performs the Stochastic Simulation Technique (SST) to
   generate the hazard curve (HC) of a Peaks-Over-Threshold (POT) sample.
   The HC can be expressed in terms of annual exceedance frequencies (AEF)
   or equivalent annual exceedance probabilities (AEP).

INPUT ARGUMENTS:
  - POT_samp: POT values (no timestamps) or second column from the first
      output argument of StormSim_POT.m; specified as a column vector.
  - Nyrs: record length in years; specified as a positive scalar.
  - HC_plt_x: predefined AEF values used to plot the HC; specified as
      a vector of positive values. The script will convert it to AEP when
      use_AEP = 1.
  - HC_tbl_x: predefined AEP values used to summarize the HC; specified as
      a vector of positive values. The script will convert it to AEF when
      use_AEP = 0.
  - HC_tbl_rsp_y: predefined response values used to summarize the HC;
      specified as a vector.
  - prc: percentage values for computing the percentiles; specified as a
      scalar or vector of positive values. Leave empty [] to apply default
      values 2.28%, 15.87%, 84.13%, 97.72%. User can enter 1 to 4 values.
      Example: prc = [2 16 84 98];
  - use_AEP: indicator for expressing the hazard as AEF or AEP. Use 1 for
      AEP, 0 for AEF. Example: use_AEP = 1;
  - GPD_TH_crit: indicator for specifying the GPD threshold option of the
      Mean Residual Life (MRL) selection process; use as follows:
       0 to evaluate all thresholds identified by the MRL method (up to 3)
       1 to only evaluate the MRL threshold selected by the lambda criterion
       2 to only evaluate the MRL threshold selected by the minimum weight criterion
      Example: GPD_TH_crit = 0;
  - TH_otlr: 'ThresholdFactor' for the fillsoutlier function used to modify
      the GPD shape parameter collection resulting from the bootstrap process.
      Specified as a nonnegative scalar. Leave empty [] to apply default value: 3.
      Refer to MATLAB Documentation for 'fillsoutlier', 'ThresholdFactor' option.
      Example: TH_otlr = 3;
  - ind_Skew: indicator for computing/adding skew tides to the storm surge.
      This applies when the input is a storm surge dataset with/without SLC.
      Use as follows:
       1 for the tool to compute and add the skew tides
       0 otherwise
      Example: ind_Skew = 1;
  - POT_samp2: POT values (no timestamps) with tides, to replace POT_samp
      without tides when evaluating storm surge with skew tides; specified
      as a column vector. Cannot be empty when ind_Skew = 1.
  - SLC: magnitude of the sea level change implicit in the storm surge input
      dataset; specified as a positive scalar. Must have same units as the
      input dataset. Example: SLC = 1.8;
  - gprMdl: Gaussian process regression (GPR) model created with the MATLAB
function 'fitrgp'; specified as an object. Refer to MATLAB
      Documentation for 'fitrgp'. Train the GPR with storm surge as a
      predictor and skew tides as the response. Leave empty when ind_Skew = 0.
  - staID: gauge station information; specified as a cell array with format:
      Col(01): gauge station ID number; as a character vector. Example: '8770570'
      Col(02): station name (OPTIONAL); as a character vector. Example: 'Sabine Pass North TX'
  - yaxis_Label: parameter name/units/datum for label of the plot y-axis;
      specified as a character vector.
      > When ind_Skew = 1: Set yaxis_Label = 'Water level' with associated units
        and datum. Example: 'Still Water Level (m, MSL)'
  - path_out: path to output folder; specified as a character vector. Leave
      empty [] to apply default: '.\SST_output\'
  - yaxis_Limits: lower and upper limits for the plot y-axis; specified as a
      vector. Leave empty [] otherwise. Example: yaxis_Limits = [0 10];
  - apply_GPD_to_SS: indicates if the hazard of small POT samples should evaluated
      with either the empirical or the GPD. A small sample has a sample size < 20 events
      and a record length < 20 years. Use 1 for GPD, 0 for empirical. Example: apply_GPD_to_SS = 1;

OUTPUT ARGUMENTS:
  - HC_emp: empirical HC; as a table array with column headings as follows:
     Response: Response vector sorted in descending order
     Rank: Rank or Weibull's plotting position
     CCDF: Complementary cumulative distribution function (CCDF)
     Hazard: hazard as AEF or AEP
     ARI: annual recurrence interval (ARI)
  - HC_plt: full HC; as a structure array with two fields: 'out' and 'MRL_Crit'. It may contain
     up to two records (one per MRL GPD threshold) with the following format:
     
     out = numerical array with the following 5 rows
        row(01): mean values
        row(02): values of 2% percentile or 1st percentage of input "prc"
        row(03): values of 16% percentile or 2nd percentage of input "prc"
        row(04): values of 84% percentile or 3rd percentage of input "prc"
        row(05): values of 98% percentile or 4th percentage of input "prc"

     MRL_Crit = character vector indicating the MRL threshold
        selection criterion used for the HC_plt.out record.

  - HC_plt_x2: same as HC_plt_x but converted to either AEP or AEF.
  - HC_tbl: summarized HC; as a structure array with two fields: 'out' and 'MRL_Crit'. It may contain
     up to two records (one per MRL GPD threshold) with the following format:

     out = numerical array with the following 5 rows
        row(01): mean values
        row(02): values of 2% percentile or 1st percentage of input "prc"
        row(03): values of 16% percentile or 2nd percentage of input "prc"
        row(04): values of 84% percentile or 3rd percentage of input "prc"
        row(05): values of 98% percentile or 4th percentage of input "prc"

     MRL_Crit = character vector indicating the MRL threshold
        selection criterion used for the HC_tbl.out record.

  - HC_tbl_rsp_x: hazard values interpolated from the HC, that correspond
     to the responses in HC_tbl_rsp_y; as a structure array with two fields: 'out' and 'MRL_Crit'. It may contain up
     to two records (one per MRL GPD threshold) with the following format:

     out = numerical array with the following 5 rows
        row(01): mean values
        row(02): values of 2% percentile or 1st percentage of input "prc"
        row(03): values of 16% percentile or 2nd percentage of input "prc"
        row(04): values of 84% percentile or 3rd percentage of input "prc"
        row(05): values of 98% percentile or 4th percentage of input "prc"

     MRL_Crit = character vector indicating the MRL threshold
        selection criterion used for the HC_tbl_rsp_x.out record.

  - MRL_output: output of Mean Residual Life (MRL) function; as a structure array with fields:
     Summary = summary of the threshold selection results, as a table array with format:
        Threshold: list of threshold values evaluated
        MeanExcess: mean excess
        Weight: weights
        WMSE: weighted mean square error (MSE)
        GPD_Shape: GPD shape parameter
        GPD_Scale: GPD scale parameter
        Events: number of events above each threshold
        Rate: sample intensity (mean annual rate of events using sample of events above threshold)

     Selection = selected threshold with other parameters, as a table array with format:
        Criterion: criterion applied by MRL method for selecting the threshold
        Threshold: selected threshold value
        id_Summary: location (row ID) of selected threshold in the Summary field
        Events: number of events above the selected threshold
        Rate: sample intensity (mean annual rate of events using sample of events above threshold)

     pd_TH_wOut = GPD threshold parameter values used in the bootstrap process
     pd_k_wOut = initial GPD shape parameter values obtained in the bootstrap process
     pd_sigma = GPD scale parameter values used in the bootstrap process
     pd_k_mod = modified GPD shape parameter values used in the bootstrap process
     eMsg = Status message indicating when the MRL methodology was not able
        to objectively determine a theshold for the GPD. In this case, the
        GPD threshold is set to 0.99 times the minimum value of the bootstrap sample.
        Otherwise, eMsg will be empty.

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

CONTRIBUTORS:
    Victor M. Gonzalez-Nieves, PE

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: revised. Updated documentation.
20201026-ERS: Added capability to include skew surge. Updated documentation and inputs.
20201125-ERS: Reviewed script. Updated documentation.
20201201-ERS: Reviewed script. Updated documentation.
20201222-ERS: Applied preprocessing rules of 1st POT sample to 2nd input POT sample.
20210303-ERS: Corrected bootstrap check plot; x values now being converted
    to AEP or AEF, depending on the case specified by user.
20210311-ERS: alpha version 3: created and applied a function to adjust the
    hazard mean and CLs, to make sure they are monotonic. Also corrected
    the interpolation scheme to compute the hazard tables.
20210324-ERS: alpha version 4: removed conversion of hazard table values
    since will be input in the correct definition.
20210325-ERS: alpha version 4: adjusted the script to read new format of
    MRL output. Bootstrap sample plot will now apply different x-axis
    labels based on value of use_AEP. Plot filename identifies by MRL
    criteria. Output structure arrays now include the MRL criteria. Update
    description of output. Negatives are changed to NaN in both hazard tables.
20210406-ERS: alpha v0.4: now storing the default GPD threshold inside MRL_output.
20210407-ERS: alpha v0.4: only bootstrap samples with >1 value are
    evaluated to avoid error of fitdist().
20210412-ERS: alpha v0.4: now discarding random samples from bootstrap
    with spurious values (>100) when skew tides are used (peak data from ADCIRC).
20210413-ERS: alpha v0.4: removed the patch of 20210413.
20210429-ERS: alpha v0.4: applied GPD Shape Parameter limits, resulting in
    representative fits of samples with less than 20 events.
20210503-ERS: alpha v0.4: added input argument apply_GPD_to_SS as an option to evaluate
    small POT samples with either GPD or empirical distribution.
20210505-ERS: alpha v0.4: added correction for the evaluation of skew tides
20210513-ERS: alpha v0.4: re-enabled the patch of 20210412. A try/catch
    statement was applied to ignore bad random samples. The preallocation of
    the arrays will make those NaN.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function [HC_emp,HC_plt,HC_plt_x2,HC_tbl,HC_tbl_rsp_x,MRL_output] = StormSim_SST_Fit(POT_samp,Nyrs,HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,POT_samp2,SLC,gprMdl,staID,yaxis_Label,path_out,yaxis_Limits,apply_GPD_to_SS)

%% Bootstrap Input Parameters
Nsim = 1e3; %Number of simulations (no less than 10,000)

%% Develop empirical CDF (using output of POT function)
POT_samp = sort(POT_samp,'descend'); %Sort POT sample in descending order
HC_emp = POT_samp; %Response vector sorted in descending order; positive values only

% Weibull Plotting Position, P = m/(n+1)
% where P = Exceedance probability
%       m = rank of descending response values, with largest equal to 1
%       n = number of response values
Nstrm_hist = length(HC_emp); % Weibull's "n"
HC_emp(:,2) = (1:Nstrm_hist)'; %Rank of descending response values
HC_emp(:,3) = HC_emp(:,2)/(Nstrm_hist+1); %Webull's "P"; NEEDS Lambda Correction

% Lambda Correction - Required for Partial Duration Series (PDS)
Lambda_hist = Nstrm_hist/Nyrs; % Lambda = sample intensity = events/year
HC_emp(:,4) = HC_emp(:,3)*Lambda_hist;

% Compute Annual Recurrence Interval (ARI)
HC_emp(:,5) = 1./HC_emp(:,4);
ecdf_y = HC_emp(:,1);


%% Perform bootstrap
% rng('default');

if ind_Skew && ~isempty(POT_samp2) && isobject(gprMdl)
    
    ecdf_y = ecdf_y - SLC;
    boot = ecdf_boot(ecdf_y,Nsim)';
    
    % Compute and add skew tides to bootstrap sample or surge
    for i=1:Nsim
        [skew_tide_mean,skew_tide_sd] = predict(gprMdl,boot(:,i));
        skew_tide_pred = skew_tide_mean + randn(length(skew_tide_mean),1).*skew_tide_sd;
        boot(:,i) = boot(:,i) + skew_tide_pred;
        boot(:,i) = sort(boot(:,i),'descend');
    end
    
    % Add the SLC amount
    boot = boot + SLC;
    
    %Substitute the empirical dist with user supplied surge + tides + SLC (WL)
    ecdf_y = sort(POT_samp2,'descend'); %Sort POT sample in descending order
    
    % Resizing
    szH = size(HC_emp,1); szy = length(ecdf_y);
    if szH>szy
        HC_emp = HC_emp(1:szy,:);
    elseif szH<szy
        ecdf_y = ecdf_y(1:szH);
    end
    HC_emp(:,1)=ecdf_y;
    boot = boot(1:length(ecdf_y),:);
else
    boot = ecdf_boot(ecdf_y,Nsim)';
end
boot(boot<0)=NaN;


%% Apply the GPD when empirical POT sample size is >20 and RL >20 yrs. Otherwise, compute HC using empirical.
if (length(ecdf_y)>=20 && Nyrs>=20) || apply_GPD_to_SS %apply GPD
    
    
    %% Apply "Mean Residual Life" Automated Threshold Detection
    MRL_output = StormSim_MRL(ecdf_y,Nyrs);
    
    
    %% MRL GPD Threshold condition
    switch GPD_TH_crit
        case 0 %user wants the 3 THs
            j=1:2;%length(TH);
            mrl_th = MRL_output.Selection.Threshold;%mrl(TH,1);
        case 1 %TH selected by lambda criterion
            j=1;
            mrl_th = MRL_output.Selection.Threshold(2);%mrl(TH,1);
        case 2 %TH selected by minimum weight criterion
            j=1;
            mrl_th = MRL_output.Selection.Threshold(1);
    end
    if isnan(mrl_th)
        j=1;
        MRL_output.Selection.Criterion(1) = {'Default'};
    else
        if iscolumn(mrl_th)
            mrl_th = mrl_th';
        end
        mrl_th = repmat(mrl_th,Nsim,1);
    end
    
    
    %% Take parameters and preallocate
    ecdf_x_adj=HC_emp(:,4); % for hazard computations, this must be AEF
    
    %pre-allocation for speed
    Resp_boot_plt=NaN(Nsim,length(HC_plt_x));
    pd_k_wOut=NaN(Nsim,1);
    pd_sigma=pd_k_wOut;
    pd_k_mod=pd_k_wOut;
    pd_TH_wOut=pd_k_wOut;
    
    HC_plt(j(end)).out=[];
    HC_plt(j(end)).MRL_Crit=[];
    HC_tbl=HC_plt;
    HC_tbl_rsp_x=HC_plt;
    
    
    %% Perform SST
    str = 'Annual Exceedance Frequency (yr^{-1})';
    if use_AEP %Convert to AEP
        HC_emp(:,4) = aef2aep(HC_emp(:,4));
        HC_plt_x2 = aef2aep(HC_plt_x);
        str = 'Annual Exceedance Probability';
    else %Convert to AEF
        HC_plt_x2 = HC_plt_x;
    end
    
    for i=j %MRL GPD Threshold loop
        
        % If the MRL didn't returned a threshold value, compute it from the
        % bootstrap samples. Then identify values above it.
        if isnan(mrl_th)
            sz = size(boot,1); %total events above threshold per simulation
            idx = ones(sz,Nsim)==1; %index of events above threshold per simulation
            idx2 = zeros(sz,Nsim)==1; %index of events below threshold per simulation
            sz = repmat(sz,1,Nsim);
            mrl_th = 0.99*min(boot,[],1,'omitnan')';
            MRL_output.Selection.Threshold(1) = mean(mrl_th,'omitnan');
            eMsg = 'No threshold found by MRL method. Default criterion applied: GPD threshold set to 0.99 times the minimum value of the bootstrap sample.';
        else
            idx = boot>mrl_th(1,i); %index of events above threshold per simulation
            sz = sum(idx,1); %total events above threshold per simulation
            idx2 = boot<=mrl_th(1,i); %index of events below threshold per simulation
            eMsg ='';
        end
        Lambda_mrl = sz/Nyrs; %annual rate of events
        
        % GPD fitting
        for k = 1:Nsim
            PEAKS_rnd = boot(:,k);
            try
                % Fit the GPD to a bootstrap data sample
                u = PEAKS_rnd(idx(:,k),1);
                pd = fitdist(u,'GeneralizedPareto','theta',mrl_th(k,i));
                pd_TH_wOut(k,i) = pd.theta;
                pd_k_wOut(k,i) = pd.k;
                pd_sigma(k,i) = pd.sigma;
            catch
            end
        end
        
        % Correction of GPD shape parameter values. Limits determined by NCNC.
        pd_k_mod(:,i) = pd_k_wOut(:,i);
        k_min=-0.5; k_max=0.3;
        pd_k_mod(pd_k_mod(:,i)<k_min,i) = k_min;
        pd_k_mod(pd_k_mod(:,i)>k_max,i) = k_max;
        
        for k = 1:Nsim
            try
                PEAKS_rnd = boot(:,k);
                
                % Compute the AEF from the GPD fit
                Resp_gpd = icdf('Generalized Pareto',1-HC_plt_x/Lambda_mrl(k),pd_k_mod(k,i),pd_sigma(k,i),pd_TH_wOut(k,i));
                AEF_gpd = HC_plt_x(~isnan(Resp_gpd));
                Resp_gpd(isnan(Resp_gpd))=[];
                
                % Take the empirical AEF from bootstrap sample
                Resp_ecdf = PEAKS_rnd(idx2(:,k));
                AEF_ecdf = ecdf_x_adj(idx2(:,k));
                
                % Merge the AEFs (empirical + fitted GPD)
                y_comb = [Resp_gpd;Resp_ecdf];
                x_comb = [AEF_gpd;AEF_ecdf];
                
                % Comvert to AEP?
                if use_AEP
                    x_comb = aef2aep(x_comb);
                end
                
                % Delete duplicates
                [~,ia,~]=unique(x_comb,'stable');
                x_comb=x_comb(ia,:);
                y_comb=y_comb(ia,:);
                
                [~,ia,~]=unique(y_comb,'stable');
                y_comb=y_comb(ia,:);
                x_comb=x_comb(ia,:);
                
                % Interpolate HC curve for table and plot
                Resp_boot_plt(k,:) = interp1(log(x_comb),y_comb,log(HC_plt_x2));
            catch
            end
        end
        
        
        %% Sample Plot for bootstrap process
        figure('Color',[1 1 1],'visible','off')
        axes('xscale','log','XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',12);
        if use_AEP
            xlim([1e-4 1]);
            XTick=[1e-4 1e-3 1e-2 1e-1 1];
        else
            xlim([1e-4 10]);
            XTick=[1e-4 1e-3 1e-2 1e-1 1 10];
        end
        if ~isempty(yaxis_Limits)
            ylim([yaxis_Limits(1) yaxis_Limits(2)]);
        end
        set(gca,'XDir','reverse','XTick',XTick)
        hold on
        
        % Historical
        scatter(HC_emp(:,4),ecdf_y,15,'g','filled','MarkerEdgeColor','k');
        
        % Resampling (last bootstrap sample)
        scatter(HC_emp(:,4),PEAKS_rnd,15,'r','filled','MarkerEdgeColor','k');
        
        % Comvert to AEP?
        if use_AEP
            AEF_ecdf = aef2aep(AEF_ecdf);
            AEF_gpd = aef2aep(AEF_gpd);
        end
        
        % Empirical
        plot(AEF_ecdf,Resp_ecdf,'y-','LineWidth',2)
        
        % MRL threshold
        th_x = interp1(y_comb,x_comb,pd_TH_wOut(k,i));
        scatter(th_x,pd_TH_wOut(k,i),15,'w','filled','MarkerEdgeColor','b');
        
        plot(AEF_gpd,Resp_gpd,'b-','LineWidth',2) % GPD with MRL
        if length(staID)==1
            title({'StormSim-SST ';['Station: ',staID{1}]},'FontSize',12);
        else
            title({'StormSim-SST ';['Station: ',staID{1},' ',staID{2}]},'FontSize',12);
        end
        xlabel(str,'FontSize',12);
        ylabel({yaxis_Label},'FontSize',12);
        legend({'Historical','Resampling','Empirical','MRL threshold','GPD'},...
            'Location','southoutside','Orientation','horizontal','NumColumns',4,'FontSize',10);
        hold off
        fname = [path_out,'SST_HC_bootCheck_',staID{1},'_TH_',MRL_output.Selection.Criterion{i},'.png'];
        saveas(gcf,fname,'png')
        
        
        %% Compute mean and percentiles
        Boot_mean_plt = mean(Resp_boot_plt,1,'omitnan');
        Boot_plt = prctile(Resp_boot_plt,prc,1);
        
        %Return an error when mean HC >= 1e3 meters
        if ind_Skew && max(Boot_mean_plt,[],'omitnan')>=1e3
            error('Values above 10^3 found in mean hazard curve')
        end
        
        Mean_p_CLs = [Boot_mean_plt;Boot_plt];
        
        % Monotonic adjustment
        for kk=1:size(Mean_p_CLs,1)
            Mean_p_CLs(kk,:) = Monotonic_adjustment(HC_plt_x2,Mean_p_CLs(kk,:));
        end
        
        % Store results for output
        HC_plt(i).out = Mean_p_CLs;
        HC_plt(i).MRL_Crit = MRL_output.Selection.Criterion{i};
        
        
        %% Interpolation to create hazard tables
        
        % preallocate
        HC_tbl_rsp_x2 = NaN(size(Mean_p_CLs,1),length(HC_tbl_rsp_y));
        HC_tbl_y2 = NaN(size(Mean_p_CLs,1),length(HC_tbl_x));
        HCmn = NaN(size(Mean_p_CLs,1),1);
        for kk=1:size(Mean_p_CLs,1)
            
            % Delete duplicates
            [~,ia,~] = unique(Mean_p_CLs(kk,:),'stable');
            dm1 = Mean_p_CLs(kk,ia); dm2 = log(HC_plt_x2(ia));
            
            % Delete NaN/Inf
            ia=isnan(dm1)|isinf(dm1); dm1(ia)=[]; dm2(ia)=[];
            
            % Interpolate
            HC_tbl_rsp_x2(kk,:) = exp(interp1(dm1,dm2,HC_tbl_rsp_y','linear','extrap'));
            HC_tbl_y2(kk,:) = interp1(dm2,dm1,log(HC_tbl_x),'linear','extrap');
            
            %interpol for 0.1 aep/aef
            HCmn(kk) = interp1(dm2,dm1,log(0.1),'linear','extrap');
            
        end
        
        % Change negatives to NaN
        HC_tbl_y2(HC_tbl_y2<0)=NaN;
        %         HC_tbl_rsp_x2(HC_tbl_rsp_x2<0)=NaN;
        HC_tbl_rsp_x2(HC_tbl_rsp_x2<1e-4)=NaN;
        
        % Compare if mean HC > 1.75* emp HC at 0.1 AEP/AEF,
        HCep = interp1(log(HC_emp(:,4)),HC_emp(:,1),log(0.1),'linear','extrap');
        str1 = {''};
        if HCmn(1)>1.75*HCep
            str1 = {'Warning: At 0.1 AEP/AEF, best estimate HC value is greater than 1.75 times the empirical HC value. Manual verification is recommended.'};
        end
        HC_plt(i).Warning = str1;
        
        % Store
        HC_tbl_rsp_x(i).out = HC_tbl_rsp_x2;
        HC_tbl_rsp_x(i).MRL_Crit = MRL_output.Selection.Criterion{i};
        HC_tbl(i).out = HC_tbl_y2;
        HC_tbl(i).MRL_Crit = MRL_output.Selection.Criterion{i};
    end
    
    
else %dont apply GPD, but compute HC + prc with empirical only
    
    
    %% Take parameters and preallocate
    Resp_boot_plt=NaN(Nsim,length(HC_plt_x));
    pd_k_wOut=NaN;
    pd_sigma=pd_k_wOut;
    pd_k_mod=pd_k_wOut;
    pd_TH_wOut=pd_k_wOut;
    HC_plt(1).out=[];
    HC_plt(1).MRL_Crit=[];
    HC_tbl=HC_plt;
    HC_tbl_rsp_x=HC_plt;
    eMsg = 'GPD not fit: POT sample size <20 and RL <20 years.';
    
    
    %% Conversion to AEF or AEP
    if use_AEP %Convert to AEP
        HC_emp(:,4) = aef2aep(HC_emp(:,4));
        HC_plt_x2 = aef2aep(HC_plt_x);
    else %Convert to AEF
        HC_plt_x2 = HC_plt_x;
    end
    
    
    %% Compute HC
    for k = 1:Nsim
        y_comb = boot(:,k);
        x_comb = HC_emp(:,4);
        
        % Delete duplicates
        [~,ia,~]=unique(x_comb,'stable');
        x_comb=x_comb(ia,:);
        y_comb=y_comb(ia,:);
        
        [~,ia,~]=unique(y_comb,'stable');
        y_comb=y_comb(ia,:);
        x_comb=x_comb(ia,:);
        
        % Interpolate HC curve for plot: lineal inside empirical, nearest neighbor outside
        HC_plt_x2a = HC_plt_x2(HC_plt_x2>=min(x_comb));
        HC_plt_x2b = HC_plt_x2(HC_plt_x2<min(x_comb));
        
        Resp_boot_plta = interp1(log(x_comb),y_comb,log(HC_plt_x2a),'linear');
        Resp_boot_pltb = interp1(log(x_comb),y_comb,log(HC_plt_x2b),'nearest','extrap');
        
        % merge results
        Resp_boot_plt(k,:) = [Resp_boot_pltb' Resp_boot_plta'];
    end
    
    
    %% Compute mean and percentiles
    Boot_mean_plt = mean(Resp_boot_plt,1,'omitnan');
    Boot_plt = prctile(Resp_boot_plt,prc,1);
    
    %For this application only: delete results if WL >= 1e3 meters
    if ind_Skew && max(Boot_mean_plt,[],'omitnan')>=1e3
        error('Values above 10^3 found in mean hazard curve')
    end
    
    Mean_p_CLs = [Boot_mean_plt;Boot_plt];
    
    % Monotonic adjustment
    for kk=1:size(Mean_p_CLs,1)
        Mean_p_CLs(kk,:) = Monotonic_adjustment(HC_plt_x2,Mean_p_CLs(kk,:));
    end
    
    % Store results for output
    HC_plt(1).out = Mean_p_CLs;
    HC_plt(1).MRL_Crit = 'None';
    
    
    %% Interpolation to create hazard tables
    
    % preallocate
    HC_tbl_rsp_x2 = NaN(size(Mean_p_CLs,1),length(HC_tbl_rsp_y));
    HC_tbl_y2 = NaN(size(Mean_p_CLs,1),length(HC_tbl_x));
    HCmn = NaN(size(Mean_p_CLs,1),1);
    for kk=1:size(Mean_p_CLs,1)
        
        % Delete duplicates
        [~,ia,~] = unique(Mean_p_CLs(kk,:),'stable');
        dm1 = Mean_p_CLs(kk,ia); dm2 = log(HC_plt_x2(ia));
        
        % Delete NaN/Inf
        ia=isnan(dm1)|isinf(dm1); dm1(ia)=[]; dm2(ia)=[];
        
        % Interpolate
        HC_tbl_rsp_x2(kk,:) = exp(interp1(dm1,dm2,HC_tbl_rsp_y','nearest','extrap'));
        HC_tbl_y2(kk,:) = interp1(dm2,dm1,log(HC_tbl_x),'nearest','extrap');
        
        %interpol for 0.1 aep/aef
        HCmn(kk) = interp1(dm2,dm1,log(0.1),'linear','extrap');
        
    end
    
    % Change negatives to NaN
    HC_tbl_y2(HC_tbl_y2<0)=NaN;
    HC_tbl_rsp_x2(HC_tbl_rsp_x2<1e-4)=NaN;
    
    % Compare if mean HC > 1.75* emp HC at 0.1 AEP/AEF,
    HCep = interp1(log(HC_emp(:,4)),HC_emp(:,1),log(0.1),'linear','extrap');
    str1 = {''};
    if HCmn(1)>1.75*HCep
        str1 = {'Warning: At 0.1 AEP/AEF, best estimate HC value is greater than 1.75 times the empirical HC value. Manual verification is recommended.'};
    end
    HC_plt(1).Warning = str1;
    
    % Store
    HC_tbl_rsp_x(1).out = HC_tbl_rsp_x2;
    HC_tbl_rsp_x(1).MRL_Crit = 'None';
    HC_tbl(1).out = HC_tbl_y2;
    HC_tbl(1).MRL_Crit = 'None';
end

%% Store parameters needed for manual GPD Shape parameter evaluation
MRL_output.pd_TH_wOut = pd_TH_wOut;
MRL_output.pd_k_wOut = pd_k_wOut;
MRL_output.pd_sigma = pd_sigma;
MRL_output.pd_k_mod = pd_k_mod;
MRL_output.Status = eMsg;
HC_emp = array2table(HC_emp,'VariableNames',{'Response','Rank','CCDF','Hazard','ARI'});

end

%% StormSim_SST_Fit_Simple.m
%{
SOFTWARE NAME:
    StormSim-SST-Fit-Simple (Statistics)

DESCRIPTION:
   This script performs the Stochastic Simulation Technique (SST) to
   generate the hazard curve (HC) of a Peaks-Over-Threshold (POT) sample.
   The HC can be expressed in terms of annual exceedance frequencies (AEF)
   or equivalent annual exceedance probabilities (AEP).

   The script is a simplified version of StormSim_SST_Fit.m, that will be
   executed when GPD_TH_crit is set to 1 or 2. It is NOT parallelized.

INPUT ARGUMENTS:
  - POT_samp: POT values (no timestamps) or second column from the first
      output argument of StormSim_POT.m; specified as a column vector.
  - Nyrs: record length in years; specified as a positive scalar.
  - HC_plt_x: predefined AEF values used to plot the HC; specified as
      a vector of positive values. The script will convert it to AEP when
      use_AEP = 1.
  - HC_tbl_x: predefined AEP values used to summarize the HC; specified as
      a vector of positive values. The script will convert it to AEF when
      use_AEP = 0.
  - HC_tbl_rsp_y: predefined response values used to summarize the HC;
      specified as a vector.
  - prc: percentage values for computing the percentiles; specified as a
      scalar or vector of positive values. Leave empty [] to apply default
      values 2.28%, 15.87%, 84.13%, 97.72%. User can enter 1 to 4 values.
      Example: prc = [2 16 84 98];
  - use_AEP: indicator for expressing the hazard as AEF or AEP. Use 1 for
      AEP, 0 for AEF. Example: use_AEP = 1;
  - GPD_TH_crit: indicator for specifying the GPD threshold option of the
      Mean Residual Life (MRL) selection process; use as follows:
       GPD_TH_crit = 0: to evaluate all thresholds identified by the MRL method (up to 3)
       GPD_TH_crit = 1: to only evaluate the MRL threshold selected by the lambda criterion
       GPD_TH_crit = 2: to only evaluate the MRL threshold selected by the minimum weight criterion
  - TH_otlr: 'ThresholdFactor' for the fillsoutlier function used to modify
      the GPD shape parameter collection resulting from the bootstrap process.
      Specified as a nonnegative scalar. Leave empty [] to apply default value: 3.
      Refer to MATLAB Documentation for 'fillsoutlier', 'ThresholdFactor' option.
      Example: TH_otlr = 3;
  - ind_Skew: indicator for computing/adding skew tides to the storm surge.
      This applies when the input is a storm surge dataset with/without SLC.
      Use as follows:
       1 for the tool to compute and add the skew tides
       0 otherwise
      Example: ind_Skew = 1;
  - POT_samp2: POT values (no timestamps) with tides, to replace POT_samp
      without tides when evaluating storm surge with skew tides; specified
      as a column vector. Cannot be empty when ind_Skew = 1.
  - SLC: magnitude of the sea level change implicit in the storm surge input
      dataset; specified as a positive scalar. Must have same units as the
      input dataset. Example: SLC = 1.8;
  - gprMdl: Gaussian process regression (GPR) model created with the MATLAB
function 'fitrgp'; specified as an object. Refer to MATLAB
      Documentation for 'fitrgp'. Train the GPR with storm surge as a
      predictor and skew tides as the response. Leave empty when ind_Skew = 0.
  - apply_GPD_to_SS: indicates if the hazard of small POT samples should evaluated
      with either the empirical or the GPD. A small sample has a sample size < 20 events
      and a record length < 20 years. Use 1 for GPD, 0 for empirical. Example: apply_GPD_to_SS = 1;
 
OUTPUT ARGUMENTS:
  - HC_emp: empirical HC; as a table array with column headings as follows:
      Response = Response vector sorted in descending order
      Rank = Rank or Weibull's plotting position
      CCDF = Complementary cumulative distribution function (CCDF)
      Hazard = hazard as AEF or AEP
      ARI = annual recurrence interval (ARI)

  - HC_plt: full HC; as a numerical array with the following format:
      row(01) mean values
      row(02) values of 2% percentile or 1st percentage of input "prc"
      row(03) values of 16% percentile or 2nd percentage of input "prc"
      row(04) values of 84% percentile or 3rd percentage of input "prc"
      row(05) values of 98% percentile or 4th percentage of input "prc"

  - HC_plt_x2: same as HC_plt_x but converted to either AEP or AEF.
  - HC_tbl: summarized HC; as a numerical array with the following format:
      row(01) mean values
      row(02) values of 2% percentile or 1st percentage of input "prc"
      row(03) values of 16% percentile or 2nd percentage of input "prc"
      row(04) values of 84% percentile or 3rd percentage of input "prc"
      row(05) values of 98% percentile or 4th percentage of input "prc"

  - HC_tbl_rsp_x: hazard values interpolated from the HC, that correspond
     to the responses in HC_tbl_rsp_y; as a numerical array with the following format:
      row(01) mean values
      row(02) values of 2% percentile or 1st percentage of input "prc"
      row(03) values of 16% percentile or 2nd percentage of input "prc"
      row(04) values of 84% percentile or 3rd percentage of input "prc"
      row(05) values of 98% percentile or 4th percentage of input "prc"

  - MRL_output: output of Mean Residual Life (MRL) function; as a structure array with fields:
     Summary = summary of the threshold selection results, as a table array with format:
        Threshold: list of threshold values evaluated
        MeanExcess: mean excess
        Weight: weights
        WMSE: weighted mean square error (MSE)
        GPD_Shape: GPD shape parameter
        GPD_Scale: GPD scale parameter
        Events: number of events above each threshold
        Rate: sample intensity (mean annual rate of events using sample of events above threshold)

     Selection = selected threshold with other parameters, as a table array with format:
        Criterion: criterion applied by MRL method for selecting the threshold
        Threshold: selected threshold value
        id_Summary: location (row ID) of selected threshold in the Summary field
        Events: number of events above the selected threshold
        Rate: sample intensity (mean annual rate of events using sample of events above threshold)

     pd_TH_wOut = GPD threshold parameter values used in the bootstrap process
     pd_k_wOut = initial GPD shape parameter values obtained in the bootstrap process
     pd_sigma = GPD scale parameter values used in the bootstrap process
     pd_k_mod = modified GPD shape parameter values used in the bootstrap process
     eMsg = Status message indicating when the MRL methodology was not able
        to objectively determine a theshold for the GPD. In this case, the
        GPD threshold is set to 0.99 times the minimum value of the bootstrap sample.
        Otherwise, eMsg will be empty.

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

CONTRIBUTORS:
    Victor M. Gonzalez-Nieves, PE

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: revised. Updated documentation.
20201026-ERS: Added capability to include skew surge. Updated documentation and inputs.
20201125-ERS: Reviewed script. Updated documentation.
20201201-ERS: Reviewed script. Updated documentation.
20201222-ERS: Applied preprocessing rules of 1st POT sample to 2nd input POT sample.
20210311-ERS: alpha version 3: created and applied a function to adjust the
    hazard mean and CLs, to make sure they are monotonic. Also corrected
    the interpolation scheme to compute the hazard tables.
20210324-ERS: alpha version 4: removed conversion of hazard table values
    since will be input in the correct definition.
20210325-ERS: alpha version 4: adjusted the script to read new format of
    MRL output. Updated description of output. Negatives are changed to NaN
    in both hazard tables.
20210406-ERS: alpha v0.4: now storing the default GPD threshold inside MRL_output.
20210407-ERS: alpha v0.4: only bootstrap samples with >1 value are
    evaluated to avoid error of fitdist().
20210412-ERS: alpha v0.4: now discarding random samples from bootstrap
    with spurious values (>100) when skew tides are used (peak data from ADCIRC).
20210413-ERS: alpha v0.4: removed the patck of 20210413.
20210503-ERS: alpha v0.4: added input argument apply_GPD_to_SS as an option to evaluate
    small POT samples with either GPD or empirical distribution. Also applied GPD
    Shape Parameter limits, resulting in representative fits of samples with less
    than 20 events.
20210505-ERS: alpha v0.4: added correction for the evaluation of skew tides
20210513-ERS: alpha v0.4: re-enabled the patch of 20210412. A try/catch
    statement was applied to ignore bad random samples. The preallocation of
    the arrays will make those NaN.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function [HC_emp,HC_plt,HC_plt_x2,HC_tbl_y,HC_tbl_rsp_x,MRL_output,str1] = StormSim_SST_Fit_Simple(POT_samp,Nyrs,HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,POT_samp2,SLC,gprMdl,apply_GPD_to_SS)

%% Bootstrap Input Parameters
Nsim = 1e3; %Number of simulations (no less than 10,000)

%% Develop empirical CDF (using output of POT function)
POT_samp = sort(POT_samp,'descend'); %Sort POT sample in descending order
HC_emp = POT_samp; %Response vector sorted in descending order; positive values only

% Weibull Plotting Position, P = m/(n+1)
% where P = Exceedance probability
%       m = rank of descending response values, with largest equal to 1
%       n = number of response values
Nstrm_hist = length(HC_emp); % Weibull's "n"
HC_emp(:,2) = (1:Nstrm_hist)'; %Rank of descending response values
HC_emp(:,3) = HC_emp(:,2)/(Nstrm_hist+1); %Webull's "P"; NEEDS Lambda Correction

% Lambda Correction - Required for Partial Duration Series (PDS)
Lambda_hist = Nstrm_hist/Nyrs; % Lambda = sample intensity = events/year
HC_emp(:,4) = HC_emp(:,3)*Lambda_hist;

% Compute Annual Recurrence Interval (ARI)
HC_emp(:,5) = 1./HC_emp(:,4);
ecdf_y = HC_emp(:,1);


%% Perform bootstrap
% rng('default');

if  ind_Skew && ~isempty(POT_samp2) && isobject(gprMdl)
    ecdf_y = ecdf_y - SLC;
    boot = ecdf_boot(ecdf_y,Nsim)';
    
    % Compute and add skew tides to bootstrap sample or surge
    for i=1:Nsim
        [skew_tide_mean,skew_tide_sd] = predict(gprMdl,boot(:,i));
        skew_tide_pred = skew_tide_mean + randn(length(skew_tide_mean),1).*skew_tide_sd;
        boot(:,i) = boot(:,i) + skew_tide_pred; %this is water level
        boot(:,i) = sort(boot(:,i),'descend');
    end
    
    % Add the SLC amount % User SLC specified? For now, this is (scalar) value is constant throughout the space
    boot = boot + SLC;
    
    %Substitute the empirical dist with user supplied surge + tides + SLC (WL)
    ecdf_y = sort(POT_samp2,'descend'); %Sort POT sample in descending order
    
    % Resizing
    szH = size(HC_emp,1); szy = length(ecdf_y);
    if szH>szy
        HC_emp = HC_emp(1:szy,:);
    elseif szH<szy
        ecdf_y = ecdf_y(1:szH);
    end
    HC_emp(:,1)=ecdf_y;
    boot = boot(1:length(ecdf_y),:);
else
    boot = ecdf_boot(ecdf_y,Nsim)';
end
boot(boot<0)=NaN;


%% Apply the GPD when empirical POT sample size is >20 and RL >20 yrs. Otherwise, compute HC using empirical.
if (length(ecdf_y)>=20 && Nyrs>=20) || apply_GPD_to_SS %apply GPD
    
    
    %% Apply "Mean Residual Life" Automated Threshold Detection
    MRL_output = StormSim_MRL(ecdf_y,Nyrs);
    
    
    %% MRL GPD Threshold condition
    if ~isnan(MRL_output.Summary.Threshold)
        switch GPD_TH_crit
            case 1 %TH selected by lambda criterion
                mrl_th = MRL_output.Selection.Threshold(2);%mrl(TH,1);
            case 2 %TH selected by minimum weight criterion
                mrl_th = MRL_output.Selection.Threshold(1);
        end
    else
        mrl_th = MRL_output.Selection.Threshold(1); %will be NaN
    end
    
    
    %% Take parameters and preallocate
    ecdf_x_adj=HC_emp(:,4);
    
    %pre-allocation for speed
    Resp_boot_plt=NaN(Nsim,length(HC_plt_x));
    pd_k_wOut=NaN(Nsim,1);
    pd_sigma=pd_k_wOut;
    pd_TH_wOut=pd_k_wOut;
    
    
    %% Perform SST
    if use_AEP %Convert to AEP
        HC_emp(:,4) = aef2aep(HC_emp(:,4));
        HC_plt_x2 = aef2aep(HC_plt_x);
    else %Convert to AEF
        HC_plt_x2 = HC_plt_x;
    end
    
    % If the MRL didn't returned a threshold value, compute it from the
    % bootstrap samples. Then identify the peaks.
    if isnan(mrl_th)
        sz = size(boot,1); %total events above threshold per simulation
        idx = ones(sz,Nsim)==1; %index of events above threshold per simulation
        idx2 = zeros(sz,Nsim)==1; %index of events below threshold per simulation
        sz = repmat(sz,1,Nsim);
        mrl_th = 0.99*min(boot,[],1,'omitnan');
        MRL_output.Selection.Threshold(1) = mean(mrl_th,'omitnan');
        MRL_output.Selection.Criterion(1) = {'Default'};
        eMsg = 'No threshold found by MRL method. Default criterion applied: GPD threshold set to 0.99 times the minimum value of the bootstrap sample.';
    else
        idx = boot>mrl_th; %index of events above threshold per simulation
        sz = sum(idx,1); %total events above threshold per simulation
        idx2 = boot<=mrl_th; %index of events below threshold per simulation
        mrl_th = repmat(mrl_th,1,Nsim);
        eMsg ='';
    end
    Lambda_mrl = sz/Nyrs; %annual rate of events
    
    % GPD fitting
    for k = 1:Nsim
        try
            PEAKS_rnd = boot(:,k);
            
            % Fit the GPD to a bootrap data sample
            u = PEAKS_rnd(idx(:,k));
            pd = fitdist(u,'GeneralizedPareto','theta',mrl_th(k));
            pd_TH_wOut(k,1) = pd.theta;
            pd_k_wOut(k,1) = pd.k;
            pd_sigma(k,1) = pd.sigma;
        catch
        end
    end
    
    % Correction of GPD shape parameter values. Limits determined by NCNC.
    pd_k_mod=pd_k_wOut;
    k_min=-0.5; k_max=0.3;
    pd_k_mod(pd_k_mod<k_min) = k_min;
    pd_k_mod(pd_k_mod>k_max) = k_max;
    
    for k = 1:Nsim
        try
            PEAKS_rnd = boot(:,k);
            
            % Compute the AEF from the GPD fit
            Resp_gpd = icdf('Generalized Pareto',1-HC_plt_x/Lambda_mrl(k),pd_k_mod(k,1),pd_sigma(k,1),pd_TH_wOut(k,1));
            AEF_gpd = HC_plt_x(~isnan(Resp_gpd));
            Resp_gpd(isnan(Resp_gpd))=[];
            
            % Compute the empirical AEF
            Resp_ecdf = PEAKS_rnd(idx2(:,k));
            AEF_ecdf = ecdf_x_adj(idx2(:,k));
            
            % Merge the AEFs (empirical + fitted GPD)
            y_comb = [Resp_gpd;Resp_ecdf];
            x_comb = [AEF_gpd;AEF_ecdf];
            
            % Comvert to AEP?
            if use_AEP
                x_comb = aef2aep(x_comb);
            end
            
            % Delete duplicates
            [~,ia,~]=unique(x_comb,'stable');
            x_comb=x_comb(ia,:);
            y_comb=y_comb(ia,:);
            
            [~,ia,~]=unique(y_comb,'stable');
            y_comb=y_comb(ia,:);
            x_comb=x_comb(ia,:);
            
            % Interpolate AEF curve for table and plot
            Resp_boot_plt(k,:) = interp1(log(x_comb),y_comb,log(HC_plt_x2));
        catch
        end
    end
    
    
    %% Compute mean and percentiles
    Boot_mean_plt = mean(Resp_boot_plt,1,'omitnan');
    Boot_plt = prctile(Resp_boot_plt,prc,1);
    
    %For this application only: delete results if WL >= 1e3 meters
    if ind_Skew && max(Boot_mean_plt,[],'omitnan')>=1e3
        error('Values above 10^3 found in mean hazard curve')
    end
    
    HC_plt = [Boot_mean_plt;Boot_plt];
    
    % Monotonic adjustment
    for kk=1:size(HC_plt,1)
        HC_plt(kk,:) = Monotonic_adjustment(HC_plt_x2,HC_plt(kk,:));
    end
    
    
    %% Interpolation to create response hazard table
    
    % preallocate
    HC_tbl_rsp_x = NaN(size(HC_plt,1),length(HC_tbl_rsp_y));
    HC_tbl_y = NaN(size(HC_plt,1),length(HC_tbl_x));
    HCmn = NaN(size(HC_plt,1),1);
    for kk=1:size(HC_plt,1)
        
        % Delete duplicates
        [~,ia,~] = unique(HC_plt(kk,:),'stable');
        dm1 = HC_plt(kk,ia); dm2 = log(HC_plt_x2(ia));
        
        % Delete NaN/Inf
        ia=isnan(dm1)|isinf(dm1); dm1(ia)=[]; dm2(ia)=[];
        
        % Interpolate
        HC_tbl_rsp_x(kk,:) = exp(interp1(dm1,dm2,HC_tbl_rsp_y','linear','extrap'));
        HC_tbl_y(kk,:) = interp1(dm2,dm1,log(HC_tbl_x),'linear','extrap');
        
        %interpol for 0.1 aep/aef
        HCmn(kk) = interp1(dm2,dm1,log(0.1),'linear','extrap');
    end
    
    % Compare if mean HC > 1.75* emp HC at 0.1 AEP/AEF,
    HCep = interp1(log(HC_emp(:,4)),HC_emp(:,1),log(0.1),'linear','extrap');
    str1 = {''};
    if HCmn(1)>1.75*HCep
        str1 = {'Warning: At 0.1 AEP/AEF, best estimate HC value is greater than 1.75 times the empirical HC value. Manual verification is recommended.'};
    end
    
    % Change negatives to NaN
    HC_tbl_y(HC_tbl_y<0)=NaN;
    HC_tbl_rsp_x(HC_tbl_rsp_x<1e-4)=NaN;
    
else %dont apply GPD, but compute HC + prc with empirical only
    
    
    %% Take parameters and preallocate
    Resp_boot_plt=NaN(Nsim,length(HC_plt_x));
    pd_k_wOut=NaN;
    pd_sigma=pd_k_wOut;
    pd_k_mod=pd_k_wOut;
    pd_TH_wOut=pd_k_wOut;
    eMsg = 'GPD not fit: POT sample size <20 and RL <20 years.';
    
    
    %% Conversion to AEF or AEP
    if use_AEP %Convert to AEP
        HC_emp(:,4) = aef2aep(HC_emp(:,4));
        HC_plt_x2 = aef2aep(HC_plt_x);
    else %Convert to AEF
        HC_plt_x2 = HC_plt_x;
    end
    
    
    %% Compute HC
    for k = 1:Nsim
        y_comb = boot(:,k);
        x_comb = HC_emp(:,4);
        
        % Delete duplicates
        [~,ia,~]=unique(x_comb,'stable');
        x_comb=x_comb(ia,:);
        y_comb=y_comb(ia,:);
        
        [~,ia,~]=unique(y_comb,'stable');
        y_comb=y_comb(ia,:);
        x_comb=x_comb(ia,:);
        
        % Interpolate HC curve for plot: lineal inside empirical, nearest neighbor outside
        HC_plt_x2a = HC_plt_x2(HC_plt_x2>=min(x_comb));
        HC_plt_x2b = HC_plt_x2(HC_plt_x2<min(x_comb));
        
        Resp_boot_plta = interp1(log(x_comb),y_comb,log(HC_plt_x2a),'linear');
        Resp_boot_pltb = interp1(log(x_comb),y_comb,log(HC_plt_x2b),'nearest','extrap');
        
        % merge results
        Resp_boot_plt(k,:) = [Resp_boot_pltb' Resp_boot_plta'];
    end
    
    
    %% Compute mean and percentiles
    Boot_mean_plt = mean(Resp_boot_plt,1,'omitnan');
    Boot_plt = prctile(Resp_boot_plt,prc,1);
    
    %For this application only: delete results if WL >= 1e3 meters
    if ind_Skew && max(Boot_mean_plt,[],'omitnan')>=1e3
        error('Values above 10^3 found in mean hazard curve')
    end
    
    HC_plt = [Boot_mean_plt;Boot_plt];
    
    % Monotonic adjustment
    for kk=1:size(HC_plt,1)
        HC_plt(kk,:) = Monotonic_adjustment(HC_plt_x2,HC_plt(kk,:));
    end
    
    
    %% Interpolation to create hazard tables
    
    % preallocate
    HC_tbl_rsp_x = NaN(size(HC_plt,1),length(HC_tbl_rsp_y));
    HC_tbl_y = NaN(size(HC_plt,1),length(HC_tbl_x));
    HCmn = NaN(size(HC_plt,1),1);
    for kk=1:size(HC_plt,1)
        
        % Delete duplicates
        [~,ia,~] = unique(HC_plt(kk,:),'stable');
        dm1 = HC_plt(kk,ia); dm2 = log(HC_plt_x2(ia));
        
        % Delete NaN/Inf
        ia=isnan(dm1)|isinf(dm1); dm1(ia)=[]; dm2(ia)=[];
        
        % Interpolate
        HC_tbl_rsp_x(kk,:) = exp(interp1(dm1,dm2,HC_tbl_rsp_y','nearest','extrap'));
        HC_tbl_y(kk,:) = interp1(dm2,dm1,log(HC_tbl_x),'nearest','extrap');
        
        %interpol for 0.1 aep/aef
        HCmn(kk) = interp1(dm2,dm1,log(0.1),'linear','extrap');
    end
    
    % Compare if mean HC > 1.75* emp HC at 0.1 AEP/AEF,
    HCep = interp1(log(HC_emp(:,4)),HC_emp(:,1),log(0.1),'linear','extrap');
    str1 = {''};
    if HCmn(1)>1.75*HCep
        str1 = {'Warning: At 0.1 AEP/AEF, best estimate HC value is greater than 1.75 times the empirical HC value. Manual verification is recommended.'};
    end
    
    % Change negatives to NaN
    HC_tbl_y(HC_tbl_y<0)=NaN;
    HC_tbl_rsp_x(HC_tbl_rsp_x<1e-4)=NaN;
end

%% Store parameters needed for manual GPD Shape parameter evaluation
MRL_output.pd_TH_wOut = pd_TH_wOut;
MRL_output.pd_k_wOut = pd_k_wOut;
MRL_output.pd_sigma = pd_sigma;
MRL_output.pd_k_mod = pd_k_mod;
MRL_output.Status = eMsg;
HC_emp = array2table(HC_emp,'VariableNames',{'Response','Rank','CCDF','Hazard','ARI'});
end

%% StormSim_SST_Fit_SimplePar.m
%{
SOFTWARE NAME:
    StormSim-SST-Fit-Simple (Statistics)

DESCRIPTION:
   This script performs the Stochastic Simulation Technique (SST) to
   generate the hazard curve (HC) of a Peaks-Over-Threshold (POT) sample.
   The HC can be expressed in terms of annual exceedance frequencies (AEF)
   or equivalent annual exceedance probabilities (AEP).

   The script is a simplified version of StormSim_SST_Fit.m, that will be
   executed when GPD_TH_crit is set to 1 or 2. It is parallelized for faster
   execution.

INPUT ARGUMENTS:
  - POT_samp: POT values (no timestamps) or second column from the first
      output argument of StormSim_POT.m; specified as a column vector.
  - Nyrs: record length in years; specified as a positive scalar.
  - HC_plt_x: predefined AEF values used to plot the HC; specified as
      a vector of positive values. The script will convert it to AEP when
      use_AEP = 1.
  - HC_tbl_x: predefined AEP values used to summarize the HC; specified as
      a vector of positive values. The script will convert it to AEF when
      use_AEP = 0.
  - HC_tbl_rsp_y: predefined response values used to summarize the HC;
      specified as a vector.
  - prc: percentage values for computing the percentiles; specified as a
      scalar or vector of positive values. Leave empty [] to apply default
      values 2.28%, 15.87%, 84.13%, 97.72%. User can enter 1 to 4 values.
      Example: prc = [2 16 84 98];
  - use_AEP: indicator for expressing the hazard as AEF or AEP. Use 1 for
      AEP, 0 for AEF. Example: use_AEP = 1;
  - GPD_TH_crit: indicator for specifying the GPD threshold option of the
      Mean Residual Life (MRL) selection process; use as follows:
       0 to evaluate all thresholds identified by the MRL method (up to 3)
       1 to only evaluate the MRL threshold selected by the lambda criterion
       2 to only evaluate the MRL threshold selected by the minimum weight criterion
      Example: GPD_TH_crit = 0;
  - TH_otlr: 'ThresholdFactor' for the fillsoutlier function used to modify
      the GPD shape parameter collection resulting from the bootstrap process.
      Specified as a nonnegative scalar. Leave empty [] to apply default value: 3.
      Refer to MATLAB Documentation for 'fillsoutlier', 'ThresholdFactor' option.
      Example: TH_otlr = 3;
  - ind_Skew: indicator for computing/adding skew tides to the storm surge.
      This applies when the input is a storm surge dataset with/without SLC.
      Use as follows:
       1 for the tool to compute and add the skew tides
       0 otherwise
      Example: ind_Skew = 1;
  - POT_samp2: POT values (no timestamps) with tides, to replace POT_samp
      without tides when evaluating storm surge with skew tides; specified
      as a column vector. Cannot be empty when ind_Skew = 1.
  - SLC: magnitude of the sea level change implicit in the storm surge input
      dataset; specified as a positive scalar. Must have same units as the
      input dataset. Example: SLC = 1.8;
  - gprMdl: Gaussian process regression (GPR) model created with the MATLAB
function 'fitrgp'; specified as an object. Refer to MATLAB
      Documentation for 'fitrgp'. Train the GPR with storm surge as a
      predictor and skew tides as the response. Leave empty when ind_Skew = 0.
  - apply_GPD_to_SS: indicates if the hazard of small POT samples should evaluated
      with either the empirical or the GPD. A small sample has a sample size < 20 events
      and a record length < 20 years. Use 1 for GPD, 0 for empirical. Example: apply_GPD_to_SS = 1;

OUTPUT ARGUMENTS:
  - HC_emp: empirical HC; as a table array with column headings as follows:
      Response = Response vector sorted in descending order
      Rank = Rank or Weibull's plotting position
      CCDF = Complementary cumulative distribution function (CCDF)
      Hazard = hazard as AEF or AEP
      ARI = annual recurrence interval (ARI)

  - HC_plt: full HC; as a numerical array with the following format:
      row(01) mean values
      row(02) values of 2% percentile or 1st percentage of input "prc"
      row(03) values of 16% percentile or 2nd percentage of input "prc"
      row(04) values of 84% percentile or 3rd percentage of input "prc"
      row(05) values of 98% percentile or 4th percentage of input "prc"

  - HC_plt_x2: same as HC_plt_x but converted to either AEP or AEF.
  - HC_tbl: summarized HC; as a numerical array with the following format:
      row(01) mean values
      row(02) values of 2% percentile or 1st percentage of input "prc"
      row(03) values of 16% percentile or 2nd percentage of input "prc"
      row(04) values of 84% percentile or 3rd percentage of input "prc"
      row(05) values of 98% percentile or 4th percentage of input "prc"

  - HC_tbl_rsp_x: hazard values interpolated from the HC, that correspond
     to the responses in HC_tbl_rsp_y; as a numerical array with the following format:
      row(01) mean values
      row(02) values of 2% percentile or 1st percentage of input "prc"
      row(03) values of 16% percentile or 2nd percentage of input "prc"
      row(04) values of 84% percentile or 3rd percentage of input "prc"
      row(05) values of 98% percentile or 4th percentage of input "prc"

  - MRL_output: output of Mean Residual Life (MRL) function; as a structure array with fields:
     Summary = summary of the threshold selection results, as a table array with format:
        Threshold: list of threshold values evaluated
        MeanExcess: mean excess
        Weight: weights
        WMSE: weighted mean square error (MSE)
        GPD_Shape: GPD shape parameter
        GPD_Scale: GPD scale parameter
        Events: number of events above each threshold
        Rate: sample intensity (mean annual rate of events using sample of events above threshold)

     Selection = selected threshold with other parameters, as a table array with format:
        Criterion: criterion applied by MRL method for selecting the threshold
        Threshold: selected threshold value
        id_Summary: location (row ID) of selected threshold in the Summary field
        Events: number of events above the selected threshold
        Rate: sample intensity (mean annual rate of events using sample of events above threshold)

     pd_TH_wOut = GPD threshold parameter values used in the bootstrap process
     pd_k_wOut = initial GPD shape parameter values obtained in the bootstrap process
     pd_sigma = GPD scale parameter values used in the bootstrap process
     pd_k_mod = modified GPD shape parameter values used in the bootstrap process
     eMsg = Status message indicating when the MRL methodology was not able
        to objectively determine a theshold for the GPD. In this case, the
        GPD threshold is set to 0.99 times the minimum value of the bootstrap sample.
        Otherwise, eMsg will be empty.

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

CONTRIBUTORS:
    Victor M. Gonzalez-Nieves, PE

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: revised. Updated documentation.
20201026-ERS: Added capability to include skew surge. Updated documentation and inputs.
20201125-ERS: Reviewed script. Updated documentation.
20201201-ERS: Reviewed script. Updated documentation.
20201222-ERS: Applied preprocessing rules of 1st POT sample to 2nd input POT sample.
20210311-ERS: alpha version 3: created and applied a function to adjust the
    hazard mean and CLs, to make sure they are monotonic. Also corrected
    the interpolation scheme to compute the hazard tables.
20210324-ERS: alpha version 4: removed conversion of hazard table values
    since will be input in the correct definition.
20210325-ERS: alpha version 4: adjusted the script to read new format of
    MRL output. Updated description of output. Negatives are changed to NaN
    in both hazard tables.
20210406-ERS: alpha v0.4: now storing the default GPD threshold inside MRL_output.
20210407-ERS: alpha v0.4: only bootstrap samples with >1 value are
    evaluated to avoid error of fitdist().
20210412-ERS: alpha v0.4: now discarding random samples from bootstrap
    with spurious values (>100) when skew tides are used (peak data from ADCIRC).
20210413-ERS: alpha v0.4: removed the patch of 20210413.
20210503-ERS: alpha v0.4: added input argument apply_GPD_to_SS as an option to evaluate
    small POT samples with either GPD or empirical distribution. Also applied GPD
    Shape Parameter limits.
20210505-ERS: alpha v0.4: added correction for the evaluation of skew tides
20210513-ERS: alpha v0.4: re-enabled the patch of 20210412. A try/catch
    statement was applied to ignore bad random samples. The preallocation of
    the arrays will make those NaN.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function [HC_emp,HC_plt,HC_plt_x2,HC_tbl_y,HC_tbl_rsp_x,MRL_output,str1] = StormSim_SST_Fit_SimplePar(POT_samp,Nyrs,HC_plt_x,HC_tbl_x,HC_tbl_rsp_y,prc,use_AEP,GPD_TH_crit,ind_Skew,POT_samp2,SLC,gprMdl,apply_GPD_to_SS)

%% Bootstrap Input Parameters
Nsim = 1e3; %Number of simulations (no less than 10,000)

%% Develop empirical CDF (using output of POT function)
POT_samp = sort(POT_samp,'descend'); %Sort POT sample in descending order
HC_emp = POT_samp; %Response vector sorted in descending order; positive values only

% Weibull Plotting Position, P = m/(n+1)
% where P = Exceedance probability
%       m = rank of descending response values, with largest equal to 1
%       n = number of response values
Nstrm_hist = length(HC_emp); % Weibull's "n"
HC_emp(:,2) = (1:Nstrm_hist)'; %Rank of descending response values
HC_emp(:,3) = HC_emp(:,2)/(Nstrm_hist+1); %Webull's "P"; NEEDS Lambda Correction

% Lambda Correction - Required for Partial Duration Series (PDS)
Lambda_hist = Nstrm_hist/Nyrs; % Lambda = sample intensity = events/year
HC_emp(:,4) = HC_emp(:,3)*Lambda_hist;

% Compute Annual Recurrence Interval (ARI)
HC_emp(:,5) = 1./HC_emp(:,4);
ecdf_y = HC_emp(:,1);


%% Perform bootstrap
% rng('default');

if  ind_Skew && ~isempty(POT_samp2) && isobject(gprMdl)
    ecdf_y = ecdf_y - SLC;
    boot = ecdf_boot(ecdf_y,Nsim)';
    
    % Compute and add skew tides to bootstrap sample or surge
    parfor i=1:Nsim
        [skew_tide_mean,skew_tide_sd] = predict(gprMdl,boot(:,i));
        skew_tide_pred = skew_tide_mean + randn(length(skew_tide_mean),1).*skew_tide_sd;
        boot(:,i) = boot(:,i) + skew_tide_pred; %this is water level
        boot(:,i) = sort(boot(:,i),'descend');
    end
    
    % Add the SLC amount % User SLC specified? For now, this is (scalar) value is constant throughout the space
    boot = boot + SLC;
    
    %Substitute the empirical dist with user supplied surge + tides + SLC (WL)
    ecdf_y= sort(POT_samp2,'descend'); %Sort POT sample in descending order
    
    szH = size(HC_emp,1); szy = length(ecdf_y);
    if szH>szy
        HC_emp = HC_emp(1:szy,:);
    elseif szH<szy
        ecdf_y = ecdf_y(1:szH);
    end
    HC_emp(:,1)=ecdf_y;
    boot = boot(1:length(ecdf_y),:);
else
    boot = ecdf_boot(ecdf_y,Nsim)';
end
boot(boot<0)=NaN;


%% Apply the GPD when empirical POT sample size is >20 and RL >20 yrs. Otherwise, compute HC using empirical.
if (length(ecdf_y)>=20 && Nyrs>=20) || apply_GPD_to_SS %apply GPD
    
    %% Apply "Mean Residual Life" Automated Threshold Detection
    MRL_output = StormSim_MRL(ecdf_y,Nyrs);
    
    
    %% MRL GPD Threshold condition
    if ~isnan(MRL_output.Summary.Threshold)
        switch GPD_TH_crit
            case 1 %TH selected by lambda criterion
                mrl_th = MRL_output.Selection.Threshold(2);%mrl(TH,1);
            case 2 %TH selected by minimum weight criterion
                mrl_th = MRL_output.Selection.Threshold(1);
        end
    else
        mrl_th = MRL_output.Selection.Threshold(1); %will be NaN
    end
    
    
    %% Take parameters and preallocate
    ecdf_x_adj = parallel.pool.Constant(HC_emp(:,4));
    
    %pre-allocation for speed
    Resp_boot_plt=NaN(Nsim,length(HC_plt_x));
    pd_k_wOut=NaN(Nsim,1);
    pd_sigma=pd_k_wOut;
    pd_TH_wOut=pd_k_wOut;
    
    
    %% Perform SST
    if use_AEP %Convert to AEP
        HC_emp(:,4) = aef2aep(HC_emp(:,4));
        HC_plt_x2 = aef2aep(HC_plt_x);
    else %Convert to AEF
        HC_plt_x2 = HC_plt_x;
    end
    
    % If the MRL didn't returned a threshold value, compute it from the
    % bootstrap samples. Then identify values above it.
    if isnan(mrl_th)
        sz = size(boot,1); %total events above threshold per simulation
        idx = ones(sz,Nsim)==1; %index of events above threshold per simulation
        idx2 = zeros(sz,Nsim)==1; %index of events below threshold per simulation
        sz = repmat(sz,1,Nsim);
        mrl_th = 0.99*min(boot,[],1,'omitnan');
        MRL_output.Selection.Threshold(1) = mean(mrl_th,'omitnan');
        MRL_output.Selection.Criterion(1) = {'Default'};
        eMsg = 'No threshold found by MRL method. Default criterion applied: GPD threshold set to 0.99 times the minimum value of the bootstrap sample.';
    else
        idx = boot>mrl_th; %index of events above threshold per simulation
        sz = sum(idx,1); %total events above threshold per simulation
        idx2 = boot<=mrl_th; %index of events below threshold per simulation
        mrl_th = repmat(mrl_th,1,Nsim);
        eMsg ='';
    end
    Lambda_mrl = sz/Nyrs; %annual rate of events
    
    % GPD fitting
    parfor k = 1:Nsim
        try
            PEAKS_rnd = boot(:,k);
            
            % Fit the GPD to a bootrap data sample
            u = PEAKS_rnd(idx(:,k));
            pd = fitdist(u,'GeneralizedPareto','theta',mrl_th(k));
            pd_TH_wOut(k,1) = pd.theta;
            pd_k_wOut(k,1) = pd.k;
            pd_sigma(k,1) = pd.sigma;
        catch
        end
    end
    
    % Correction of GPD shape parameter values. Limits determined by NCNC.
    pd_k_mod = pd_k_wOut;
    k_min=-0.5; k_max=0.3;
    pd_k_mod(pd_k_mod<k_min) = k_min;
    pd_k_mod(pd_k_mod>k_max) = k_max;
    
    parfor k = 1:Nsim
        try
            PEAKS_rnd = boot(:,k);
            
            % Compute the AEF from the GPD fit
            Resp_gpd = icdf('Generalized Pareto',1-HC_plt_x/Lambda_mrl(k),pd_k_mod(k,1),pd_sigma(k,1),pd_TH_wOut(k,1));
            AEF_gpd = HC_plt_x(~isnan(Resp_gpd));
            Resp_gpd(isnan(Resp_gpd))=[];
            
            % Compute the empirical AEF
            Resp_ecdf = PEAKS_rnd(idx2(:,k));
            AEF_ecdf = ecdf_x_adj.Value(idx2(:,k));
            
            % Merge the AEFs (empirical + fitted GPD)
            y_comb = [Resp_gpd;Resp_ecdf];
            x_comb = [AEF_gpd;AEF_ecdf];
            
            % Comvert to AEP?
            if use_AEP
                x_comb = aef2aep(x_comb);
            end
            
            % Delete duplicates
            [~,ia,~]=unique(x_comb,'stable');
            x_comb=x_comb(ia,:);
            y_comb=y_comb(ia,:);
            
            [~,ia,~]=unique(y_comb,'stable');
            y_comb=y_comb(ia,:);
            x_comb=x_comb(ia,:);
            
            % Interpolate AEF curve for table and plot
            Resp_boot_plt(k,:) = interp1(log(x_comb),y_comb,log(HC_plt_x2));
        catch
        end
    end
    
    
    %% Compute mean and percentiles
    Boot_mean_plt = mean(Resp_boot_plt,1,'omitnan');
    Boot_plt = prctile(Resp_boot_plt,prc,1);
    
    %For this application only: delete results if WL >= 1e3 meters
    if ind_Skew && max(Boot_mean_plt,[],'omitnan')>=1e3
        error('Values above 10^3 found in mean hazard curve')
    end
    
    HC_plt = [Boot_mean_plt;Boot_plt];
    
    % Monotonic adjustment
    parfor kk=1:size(HC_plt,1)
        HC_plt(kk,:) = Monotonic_adjustment(HC_plt_x2,HC_plt(kk,:));
    end
    
    
    %% Interpolation to create response hazard table
    
    % preallocate
    HC_tbl_rsp_x = NaN(size(HC_plt,1),length(HC_tbl_rsp_y));
    HC_tbl_y = NaN(size(HC_plt,1),length(HC_tbl_x));
    HCmn = NaN(size(HC_plt,1),1);
    for kk=1:size(HC_plt,1)
        
        % Delete duplicates
        [~,ia,~] = unique(HC_plt(kk,:),'stable');
        dm1 = HC_plt(kk,ia); dm2 = log(HC_plt_x2(ia));
        
        % Delete NaN/Inf
        ia=isnan(dm1)|isinf(dm1); dm1(ia)=[]; dm2(ia)=[];
        
        % Interpolate
        HC_tbl_rsp_x(kk,:) = exp(interp1(dm1,dm2,HC_tbl_rsp_y','linear','extrap'));
        HC_tbl_y(kk,:) = interp1(dm2,dm1,log(HC_tbl_x),'linear','extrap');
        
        %interpol for 0.1 aep/aef
        HCmn(kk) = interp1(dm2,dm1,log(0.1),'linear','extrap');
    end
    
    % Compare if mean HC > 1.75* emp HC at 0.1 AEP/AEF,
    HCep = interp1(log(HC_emp(:,4)),HC_emp(:,1),log(0.1),'linear','extrap');
    str1 = {''};
    if HCmn(1)>1.75*HCep
        str1 = {'Warning: At 0.1 AEP/AEF, best estimate HC value is greater than 1.75 times the empirical HC value. Manual verification is recommended.'};
    end
    
    % Change negatives to NaN
    HC_tbl_y(HC_tbl_y<0)=NaN;
    HC_tbl_rsp_x(HC_tbl_rsp_x<1e-4)=NaN;
    
else %dont apply GPD, but compute HC + prc with empirical only
    
    %% Take parameters and preallocate
    Resp_boot_plt=NaN(Nsim,length(HC_plt_x));
    pd_k_wOut=NaN;
    pd_sigma=pd_k_wOut;
    pd_k_mod=pd_k_wOut;
    pd_TH_wOut=pd_k_wOut;
    eMsg = 'GPD not fit: POT sample size <20 and RL <20 years.';
    
    
    %% Conversion to AEF or AEP
    if use_AEP %Convert to AEP
        HC_emp(:,4) = aef2aep(HC_emp(:,4));
        HC_plt_x2 = aef2aep(HC_plt_x);
    else %Convert to AEF
        HC_plt_x2 = HC_plt_x;
    end
    
    
    %% Compute HC
    p_HC_emp = parallel.pool.Constant(HC_emp);
    p_HC_plt_x2 = parallel.pool.Constant(HC_plt_x2);
    parfor k = 1:Nsim
        y_comb = boot(:,k);
        x_comb = p_HC_emp.Value(:,4);
        
        % Delete duplicates
        [~,ia,~]=unique(x_comb,'stable');
        x_comb=x_comb(ia,:);
        y_comb=y_comb(ia,:);
        
        [~,ia,~]=unique(y_comb,'stable');
        y_comb=y_comb(ia,:);
        x_comb=x_comb(ia,:);
        
        % Interpolate HC curve for plot: lineal inside empirical, nearest neighbor outside
        HC_plt_x2a = p_HC_plt_x2.Value(p_HC_plt_x2.Value>=min(x_comb));
        HC_plt_x2b = p_HC_plt_x2.Value(p_HC_plt_x2.Value<min(x_comb));
        
        Resp_boot_plta = interp1(log(x_comb),y_comb,log(HC_plt_x2a),'linear');
        Resp_boot_pltb = interp1(log(x_comb),y_comb,log(HC_plt_x2b),'nearest','extrap');
        
        % merge results
        Resp_boot_plt(k,:) = [Resp_boot_pltb' Resp_boot_plta'];
    end
    
    
    %% Compute mean and percentiles
    Boot_mean_plt = mean(Resp_boot_plt,1,'omitnan');
    Boot_plt = prctile(Resp_boot_plt,prc,1);
    
    %For this application only: delete results if WL >= 1e3 meters
    if ind_Skew && max(Boot_mean_plt,[],'omitnan')>=1e3
        error('Values above 10^3 found in mean hazard curve')
    end
    
    HC_plt = [Boot_mean_plt;Boot_plt];
    
    % Monotonic adjustment
    parfor kk=1:size(HC_plt,1)
        HC_plt(kk,:) = Monotonic_adjustment(HC_plt_x2,HC_plt(kk,:));
    end
    
    
    %% Interpolation to create hazard tables
    
    % preallocate
    HC_tbl_rsp_x = NaN(size(HC_plt,1),length(HC_tbl_rsp_y));
    HC_tbl_y = NaN(size(HC_plt,1),length(HC_tbl_x));
    HCmn = NaN(size(HC_plt,1),1);
    for kk=1:size(HC_plt,1)
        
        % Delete duplicates
        [~,ia,~] = unique(HC_plt(kk,:),'stable');
        dm1 = HC_plt(kk,ia); dm2 = log(HC_plt_x2(ia));
        
        % Delete NaN/Inf
        ia=isnan(dm1)|isinf(dm1); dm1(ia)=[]; dm2(ia)=[];
        
        % Interpolate
        HC_tbl_rsp_x(kk,:) = exp(interp1(dm1,dm2,HC_tbl_rsp_y','nearest','extrap'));
        HC_tbl_y(kk,:) = interp1(dm2,dm1,log(HC_tbl_x),'nearest','extrap');
        
        %interpol for 0.1 aep/aef
        HCmn(kk) = interp1(dm2,dm1,log(0.1),'linear','extrap');
    end
    
    % Compare if mean HC > 1.75* emp HC at 0.1 AEP/AEF,
    HCep = interp1(log(HC_emp(:,4)),HC_emp(:,1),log(0.1),'linear','extrap');
    str1 = {''};
    if HCmn(1)>1.75*HCep
        str1 = {'Warning: At 0.1 AEP/AEF, best estimate HC value is greater than 1.75 times the empirical HC value. Manual verification is recommended.'};
    end
    
    % Change negatives to NaN
    HC_tbl_y(HC_tbl_y<0)=NaN;
    HC_tbl_rsp_x(HC_tbl_rsp_x<1e-4)=NaN;
end

%% Store parameters needed for manual GPD Shape parameter evaluation
MRL_output.pd_TH_wOut = pd_TH_wOut;
MRL_output.pd_k_wOut = pd_k_wOut;
MRL_output.pd_sigma = pd_sigma;
MRL_output.pd_k_mod = pd_k_mod;
MRL_output.Status = eMsg;
HC_emp = array2table(HC_emp,'VariableNames',{'Response','Rank','CCDF','Hazard','ARI'});
end

%% StormSim_SST_Plot.m
%{
SOFTWARE NAME:
    StormSim-SST-Plot (Plot Statistics)

DESCRIPTION:
   This script plots the outputs of the StormSim_MRL.m function and the
   hazard curve generated with the StormSim_SST_Fit.m

INPUT ARGUMENTS:
  - HC_emp: output from StormSim_SST_Fit.m
  - HC_plt: output from StormSim_SST_Fit.m
  - HC_plt_x: output from StormSim_SST_Fit.m
  - MRL_out: output from StormSim_SST_Fit.m
  - prc: percentage values for computing the percentiles; specified as a
      scalar or vector of positive values. Leave empty [] to apply default
      values 2%,16%,84%,98%. User can enter 1 to 4 values.
      Example: prc = [2 16 84 98];
  - use_AEP: indicator for expressing the hazard as AEF or AEP. Use 1 for
      AEP, 0 for AEF. Example: use_AEP = 1;
  - staID: gauge station information; specified as a cell array with format:
      Col(01): gauge station ID number; as a character vector. Example: '8770570'
      Col(02): station name (OPTIONAL); as a character vector. Example: 'Sabine Pass North TX'
  - yaxis_Label: parameter name/units/datum for label of the plot y-axis;
      specified as a character vector. Example: 'Still Water Level (m, MSL)'
  - path_out: path to output folder; specified as a character vector. Leave
      empty [] to apply default: '.\SST_plots\'
  - yaxis_Limits: lower and upper limits for the plot y-axis; specified as a
      vector. Leave empty [] otherwise. Example: yaxis_Limits = [0 10];
  - GPD_TH_crit: indicator for specifying the GPD threshold option of the
      Mean Residual Life (MRL) selection process; specified as a scalar. Use as follows:
      > GPD_TH_crit = 0: to evaluate all thresholds identified by the MRL method
      > GPD_TH_crit = 1: to evaluate the MRL threshold selected by the Sample Intensity criterion
      > GPD_TH_crit = 2: to evaluate the MRL threshold selected by the minimum WMSE criterion
  - a: MATLAB release; used to select either "exportgraphics" or "saveas"
      as the method for saving the plots.

OUTPUT ARGUMENTS:
   The output are plots automatically stored in the output path. The plots
   are of hazard curves with confidence levels the specified in input "prc".

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: revised. Updated documentation.
20210324-ERS: alpha version 4: updated the x label for the hazard plots.
20210325-ERS: alpha version 4: adjusted the script to read new format of
    MRL output. Hazard plots will now apply different x-axis labels based on
    value of use_AEP. Plot filename identifies by MRL criteria.
20210406-ERS: alpha v0.4: modified to account for the default GPD threshold.
20210430-ERS: alpha v0.4: updated.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function StormSim_SST_Plot(HC_plt,HC_emp,MRL_output,prc,staID,yaxis_Label,path_out,yaxis_Limits,use_AEP,GPD_TH_crit,a,HC_plt_x)

% Colors for percentiles plots
cs={'r-.','b--','b--','r-.'};
if length(prc)<4,cs={'r-.','b-.','m-.'};end
prc=round(prc);

if ~strcmp(HC_plt(1).MRL_Crit,'None') %plot these only when the GPD fit occurred
    
    %% Plot Mean residual life results
    mrl = MRL_output.Summary;
    TH = MRL_output.Selection.Threshold;
    crit = MRL_output.Selection.Criterion;
    if strcmp(MRL_output.Status,'')
        
        figure('units','inches','Position',[1 1 6 6],'Color',[1 1 1],'visible','off');
        axes('XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',12);
        
        subplot(5,1,1)
        hold on
        ylim([round(min(mrl.Rate)-0.01,2) round(max(mrl.Rate)+0.01,2)]);
        yl=ylim;
        plot(mrl.Threshold,mrl.Rate,'g-','LineWidth',2);
        str={'k:','k--'};
        if ~isempty(TH)
            for i=1:length(TH)
                h(i).p = plot([TH(i) TH(i)],[yl(1) yl(2)],str{i},'LineWidth',2); %#ok<AGROW>
            end
            legend([h.p],crit,'Location','NorthEast')
        end
        
        if length(staID)==1
            title({'StormSim-SST - Mean Residual Life';['Station: ',staID{1}]},'FontSize',12);
        else
            title({'StormSim-SST ';['Station: ',staID{1},' ',staID{2}]},'FontSize',12);
        end
        
        ylabel({'Events';'per year'},'FontSize',12);
        hold off
        
        subplot(5,1,2)
        hold on
        ylim([round(min(mrl.WMSE)-0.01,2) round(max(mrl.WMSE)+0.01,2)]);
        yl=ylim;
        plot(mrl.Threshold,mrl.WMSE,'LineWidth',2,'Color','m');
        
        if ~isempty(TH)
            for i=1:length(TH)
                plot([TH(i) TH(i)],[yl(1) yl(2)],str{i},'LineWidth',2);
            end
        end
        ylabel({'Weighted';'MSE'},'FontSize',12);
        hold off
        
        subplot(5,1,3)
        hold on
        ylim([round(min(mrl.MeanExcess)-0.01,2) round(max(mrl.MeanExcess)+0.01,2)]);
        yl = ylim;
        plot(mrl.Threshold,mrl.MeanExcess,'LineWidth',2,'Color','k');
        
        if ~isempty(TH)
            for i=1:length(TH)
                plot([TH(i) TH(i)],[yl(1) yl(2)],str{i},'LineWidth',2);
            end
        end
        
        ylabel({'Mean';'Excess'},'FontSize',12);
        hold off
        
        subplot(5,1,4)
        hold on
        ylim([round(min(mrl.GPD_Shape)-0.01,2) round(max(mrl.GPD_Shape)+0.01,2)]);
        yl = ylim;
        plot(mrl.Threshold,mrl.GPD_Shape,'LineWidth',2,'Color','b');
        if ~isempty(TH)
            for i=1:length(TH)
                plot([TH(i) TH(i)],[yl(1) yl(2)],str{i},'LineWidth',2);
            end
        end
        ylabel({'Shape';'Parameter'},'FontSize',12);
        hold off
        
        subplot(5,1,5)
        hold on
        ylim([round(min(mrl.GPD_Scale)-0.01,2) round(max(mrl.GPD_Scale)+0.01,2)]);
        yl = ylim;
        plot(mrl.Threshold,mrl.GPD_Scale,'LineWidth',2,'Color','r');
        if ~isempty(TH)
            for i=1:length(TH)
                plot([TH(i) TH(i)],[yl(1) yl(2)],str{i},'LineWidth',2);
            end
        end
        xlabel(['Threshold: ',yaxis_Label],'FontSize',12);
        ylabel({'Scale';'Parameter'},'FontSize',12);
        hold off
        
        fname = [path_out,'SST_','MRL_',staID{1},'.png'];
        switch a
            case 0
                exportgraphics(gcf,fname,'Resolution',150)
            case 1
                saveas(gcf,fname,'png')
        end
    end
    
    %% Plot GPD parameters from bootstrap per threshold
    mn_k = mean(MRL_output.pd_k_wOut,1,'omitnan');
    mn_k2 = mean(MRL_output.pd_k_mod,1,'omitnan');
    if isempty(TH),sz=1;else,sz=length(TH);end
    pObj(sz).p=[];
    
    figure('Color','w','visible','off')
    for j=1:sz
        switch crit{j}
            case 'CritWMSE'
                str = 'WMSE Criterion';
            case 'CritSI'
                str = 'Sample Intensity Criterion';
            case 'Default'
                str = 'Default';
            otherwise
                str = 'None';
        end
        
        pObj(j).p = subplot(sz,1,j);
        subplot(sz,1,j)
        hold on
        histogram(MRL_output.pd_k_wOut(:,j),'BinWidth',.05,'FaceColor','b') %original
        histogram(MRL_output.pd_k_mod(:,j),'BinWidth',.05,'FaceColor','r') %w/outliers filled
        y1=ylim;
        plot([mn_k(j) mn_k(j)],[y1(1) y1(2)],'b-','LineWidth',2)
        plot([mn_k2(j) mn_k2(j)],[y1(1) y1(2)],'r-.','LineWidth',2)
        xlabel(['GPD Shape Parameter - MRL Threshold - ',str],'FontSize',12);
        ylabel('Count','FontSize',12)
        hold off
    end
    
    legend(pObj(1).p,{'hist w/outliers','hist w/o outliers','mean w/outliers',...
        'mean w/o outliers'},'Location','northwest','FontSize',8)
    fname = [path_out,'SST_CompareGPDShape_',staID{1},'.png'];
    switch a
        case 0
            exportgraphics(gcf,fname,'Resolution',150)
        case 1
            saveas(gcf,fname,'png')
    end
    
    %% Plot Hazard Curve
    j=1:length(TH);
    for k=j %TH loop
        
        % Take probs and mean values
        Boot_mean_plt = HC_plt(k).out(1,:);
        
        % Take percentiles
        Boot_plt=[];
        if size(HC_plt(k).out,1)>2
            Boot_plt = HC_plt(k).out(2:end,:);
        end
        
        % Do plot
        figure('Color',[1 1 1],'visible','off')
        axes('xscale','log','XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',12);
        if use_AEP
            xlim([1e-4 1]);
            XTick=[1e-4 1e-3 1e-2 1e-1 1];
        else
            xlim([1e-4 10]);
            XTick=[1e-4 1e-3 1e-2 1e-1 1 10];
        end
        
        if ~isempty(yaxis_Limits)
            ylim([min(yaxis_Limits) max(yaxis_Limits)]);
        end
        set(gca,'XDir','reverse','XTick',XTick)
        hold on
        
        % Process the percentiles
        pObj=struct('o',[],'n',[],'L','');
        for i=1:length(prc)
            pObj(i).o = plot(HC_plt_x,Boot_plt(i,:),cs{i},'LineWidth',2);
            pObj(i).n = prc(i);
            pObj(i).L = {['CL',int2str(prc(i)),'%']};
        end
        pObj_t = struct2table(pObj);
        pObj_t = sortrows(pObj_t,'n','descend'); % sort the table by 'n'
        
        % Take percentiles
        t1 = pObj_t(pObj_t.n>50,:); t1 = table2struct(t1); %above the mean
        t2 = pObj_t(pObj_t.n<50,:); t2 = table2struct(t2); %below the mean
        
        % Mean and Empirical
        h1 = plot(HC_plt_x,Boot_mean_plt,'k-','LineWidth',2); % Mean
        h2 = scatter(HC_emp.Hazard,HC_emp.Response,10,'g','filled','MarkerEdgeColor','k'); % Historical
        
        % Title
        if length(staID)==1
            title({'StormSim-SST ';['Station: ',staID{1}]},'FontSize',12);
        else
            title({'StormSim-SST ';['Station: ',staID{1},' ',staID{2}]},'FontSize',12);
        end
        
        % Axes labels
        if use_AEP
            xlabel('Annual Exceedance Probability','FontSize',12);
        else
            xlabel('Annual Exceedance Frequency (yr^{-1})','FontSize',12);
        end
        ylabel({yaxis_Label},'FontSize',12);
        
        % Legend
        legend([[t1.o],h1,[t2.o],h2],{t1.L,'Mean',t2.L,'Empirical'},...
            'Location','southoutside','Orientation','horizontal','NumColumns',5,'FontSize',10);
        hold off
        
        % Filename
        if GPD_TH_crit==0
            fname = [path_out,'SST_','HC_',staID{1},'_TH_',crit{k},'.png'];
        else
            fname = [path_out,'SST_','HC_',staID{1},'.png'];
        end
        
        % Save figure
        switch a
            case 0
                exportgraphics(gcf,fname,'Resolution',150)
            case 1
                saveas(gcf,fname,'png')
        end
    end
else
    
    %% Plot Hazard Curve
    
    % Take probs and mean values
    Boot_mean_plt = HC_plt(1).out(1,:);
    
    % Take percentiles
    Boot_plt=[];
    if size(HC_plt(1).out,1)>2
        Boot_plt = HC_plt(1).out(2:end,:);
    end
    
    % Do plot
    figure('Color',[1 1 1],'visible','off')
    axes('xscale','log','XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',12);
    if use_AEP
        xlim([1e-4 1]);
        XTick=[1e-4 1e-3 1e-2 1e-1 1];
    else
        xlim([1e-4 10]);
        XTick=[1e-4 1e-3 1e-2 1e-1 1 10];
    end
    
    if ~isempty(yaxis_Limits)
        ylim([min(yaxis_Limits) max(yaxis_Limits)]);
    end
    set(gca,'XDir','reverse','XTick',XTick)
    hold on
    
    % Process the percentiles
    pObj=struct('o',[],'n',[],'L','');
    for i=1:length(prc)
        pObj(i).o = plot(HC_plt_x,Boot_plt(i,:),cs{i},'LineWidth',2);
        pObj(i).n = prc(i);
        pObj(i).L = {['CL',int2str(prc(i)),'%']};
    end
    pObj_t = struct2table(pObj);
    pObj_t = sortrows(pObj_t,'n','descend'); % sort the table by 'n'
    
    % Take percentiles
    t1 = pObj_t(pObj_t.n>50,:); t1 = table2struct(t1); %above the mean
    t2 = pObj_t(pObj_t.n<50,:); t2 = table2struct(t2); %below the mean
    
    % Mean and Empirical
    h1 = plot(HC_plt_x,Boot_mean_plt,'k-','LineWidth',2); % Mean
    h2 = scatter(HC_emp.Hazard,HC_emp.Response,10,'g','filled','MarkerEdgeColor','k'); % Historical
    
    % Title
    if length(staID)==1
        title({'StormSim-SST ';['Station: ',staID{1}]},'FontSize',12);
    else
        title({'StormSim-SST ';['Station: ',staID{1},' ',staID{2}]},'FontSize',12);
    end
    
    % Axes labels
    if use_AEP
        xlabel('Annual Exceedance Probability','FontSize',12);
    else
        xlabel('Annual Exceedance Frequency (yr^{-1})','FontSize',12);
    end
    ylabel({yaxis_Label},'FontSize',12);
    
    % Legend
    legend([[t1.o],h1,[t2.o],h2],{t1.L,'Mean',t2.L,'Empirical'},...
        'Location','southoutside','Orientation','horizontal','NumColumns',5,'FontSize',10);
    hold off
    
    % Filename
    fname = [path_out,'SST_','HC_',staID{1},'_TH_None.png'];
    
    % Save figure
    switch a
        case 0
            exportgraphics(gcf,fname,'Resolution',150)
        case 1
            saveas(gcf,fname,'png')
    end
end
end

%% StormSim_SST_Plot_Simple.m
%{
SOFTWARE NAME:
    StormSim-SST-Plot (Plot Statistics)

DESCRIPTION:
   This script plots the outputs of the function StormSim_MRL.m and the
   hazard curve generated with the function StormSim_SST_Fit.m

INPUT ARGUMENTS:
  - HC_emp: output from StormSim_SST_Fit.m
  - HC_plt: output from StormSim_SST_Fit.m
  - HC_plt_x: output from StormSim_SST_Fit.m
  - MRL_out: output from StormSim_SST_Fit.m
  - prc: percentage values for computing the percentiles; specified as a
      scalar or vector of positive values. Leave empty [] to apply default
      values 2.28%, 15.87%, 84.13%, 97.72%. User can enter 1 to 4 values.
      Example: prc = [2 16 84 98];
  - use_AEP: indicator for expressing the hazard as AEF or AEP. Use 1 for
      AEP, 0 for AEF. Example: use_AEP = 1;
  - staID: gauge station information; specified as a cell array with format:
      Col(01): gauge station ID number; as a character vector. Example: '8770570'
      Col(02): station name (OPTIONAL); as a character vector. Example: 'Sabine Pass North TX'
  - yaxis_Label: parameter name/units/datum for label of the plot y-axis;
      specified as a character vector. Example: 'Still Water Level (m, MSL)'
  - path_out: path to output folder; specified as a character vector. Leave
      empty [] to apply default: '.\SST_output\'
  - yaxis_Limits: lower and upper limits for the plot y-axis; specified as a
      vector. Leave empty [] otherwise. Example: yaxis_Limits = [0 10];
  - GPD_TH_crit: indicator for specifying the GPD threshold option of the
      Mean Residual Life (MRL) selection process; specified as a scalar. Use as follows:
      > GPD_TH_crit = 0: to evaluate all thresholds identified by the MRL method
      > GPD_TH_crit = 1: to evaluate the MRL threshold selected by the Sample Intensity criterion
      > GPD_TH_crit = 2: to evaluate the MRL threshold selected by the minimum WMSE criterion
  - a: MATLAB release; used to select either "exportgraphics" or "saveas"
      as the method for saving the plots.

OUTPUT ARGUMENTS:
   The output are plots automatically stored in the output path. The plots
   are of hazard curves with confidence levels the specified in input "prc".

AUTHORS:
    Norberto C. Nadal-Caraballo, PhD (NCNC)
    Efrain Ramos-Santiago (ERS)

HISTORY OF REVISIONS:
20200903-ERS: revised.
20201015-ERS: revised. Updated documentation.
20201201-ERS: Reviewed script. Updated documentation.
20201218-ERS: Script can now plot data related MRL thresholds selected by min weight criterion.
20210324-ERS: alpha version 4: Updated the x label for the hazard plots.
20210325-ERS: alpha version 4: adjusted the script to read new format of
    MRL output. Hazard plots will now apply different x-axis labels based on
    value of use_AEP. Plot filename identifies by MRL criteria.
20210406-ERS: alpha v0.4: modified to account for the default GPD threshold.
20210503-ERS: alpha v0.4: updated.

***************  ALPHA  VERSION  **  FOR INTERNAL TESTING ONLY ************
%}
function StormSim_SST_Plot_Simple(HC_emp,HC_plt,HC_plt_x,MRL_output,prc,use_AEP,staID,yaxis_Label,path_out,yaxis_Limits,a,GPD_TH_crit)
%% A few check points

% Colors for percentiles plots
cs={'r-.','b--','b--','r-.'};
if length(prc)<4,cs={'r-.','b-.','m-.'};end
prc=round(prc);

%% Plot Mean residual life results
if strcmp(MRL_output.Status,'')
    
    % Take inputs
    mrl = MRL_output.Summary;
    
    if ~isnan(MRL_output.Summary.Threshold)
        switch GPD_TH_crit
            case 1 %TH selected by lambda criterion
                mrl_th = MRL_output.Selection.Threshold(2);%mrl(TH,1);
                crit = MRL_output.Selection.Criterion(2);
                
            case 2 %TH selected by minimum weight criterion
                mrl_th = MRL_output.Selection.Threshold(1);
                crit = MRL_output.Selection.Criterion(1);
        end
    else
        mrl_th = MRL_output.Selection.Threshold(1); %will be NaN
    end
    
    figure('units','inches','Position',[1 1 6 6],'Color',[1 1 1],'visible','off');
    axes('XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',12);
    
    subplot(5,1,1)
    hold on
    ylim([round(min(mrl.Rate)-0.01,2) round(max(mrl.Rate)+0.01,2)]);
    yl=ylim;
    plot(mrl.Threshold,mrl.Rate,'g-','LineWidth',2);
    if ~isempty(mrl_th)
        h=plot([mrl_th mrl_th],[yl(1) yl(2)],'k:','LineWidth',2);
        legend(h,crit,'Location','NorthEast')
    end
    if length(staID)==1
        title({'StormSim-SST - Mean Residual Life';['Station: ',staID{1}]},'FontSize',12);
    else
        title({'StormSim-SST ';['Station: ',staID{1},' ',staID{2}]},'FontSize',12);
    end
    ylabel({'Events';'per year'},'FontSize',12);
    hold off
    
    subplot(5,1,2)
    hold on
    ylim([round(min(mrl.WMSE)-0.01,2) round(max(mrl.WMSE)+0.01,2)]);
    yl=ylim;
    plot(mrl.Threshold,mrl.WMSE,'LineWidth',2,'Color','m');
    if ~isempty(mrl_th)
        plot([mrl_th mrl_th],[yl(1) yl(2)],'k:','LineWidth',2);
    end
    ylabel({'Weighted';'MSE'},'FontSize',12);
    hold off
    
    subplot(5,1,3)
    hold on
    ylim([round(min(mrl.MeanExcess)-0.01,2) round(max(mrl.MeanExcess)+0.01,2)]);
    yl = ylim;
    plot(mrl.Threshold,mrl.MeanExcess,'LineWidth',2,'Color','k');
    if ~isempty(mrl_th)
        plot([mrl_th mrl_th],[yl(1) yl(2)],'k:','LineWidth',2);
    end
    ylabel({'Mean';'Excess'},'FontSize',12);
    hold off
    
    subplot(5,1,4)
    hold on
    ylim([round(min(mrl.GPD_Shape)-0.01,2) round(max(mrl.GPD_Shape)+0.01,2)]);
    yl = ylim;
    plot(mrl.Threshold,mrl.GPD_Shape,'LineWidth',2,'Color','b');
    if ~isempty(mrl_th)
        plot([mrl_th mrl_th],[yl(1) yl(2)],'k:','LineWidth',2);
    end
    ylabel({'Shape';'Parameter'},'FontSize',12);
    hold off
    
    subplot(5,1,5)
    hold on
    ylim([round(min(mrl.GPD_Scale)-0.01,2) round(max(mrl.GPD_Scale)+0.01,2)]);
    yl = ylim;
    plot(mrl.Threshold,mrl.GPD_Scale,'LineWidth',2,'Color','r');
    if ~isempty(mrl_th)
        plot([mrl_th mrl_th],[yl(1) yl(2)],'k:','LineWidth',2);
    end
    xlabel(['Threshold: ',yaxis_Label],'FontSize',12);
    ylabel({'Scale';'Parameter'},'FontSize',12);
    hold off
    
    fname = [path_out,'SST_','MRL_',staID{1},'.png'];
    switch a
        case 0
            exportgraphics(gcf,fname,'Resolution',150)
        case 1
            saveas(gcf,fname,'png')
    end
end

%% Plot hazard curve
Boot_mean_plt = HC_plt(1,:); % Take probs and mean values
Boot_plt=[]; if size(HC_plt,1)>1,Boot_plt=HC_plt(2:end,:);end % Take percentiles

% Do plot
figure('Color',[1 1 1],'visible','off')
axes('xscale','log','XGrid','on','XMinorTick','on','YGrid','on','YMinorTick','on','FontSize',12);
if use_AEP
    xlim([1e-4 1]);
    XTick=[1e-4 1e-3 1e-2 1e-1 1];
else
    xlim([1e-4 10]);
    XTick=[1e-4 1e-3 1e-2 1e-1 1 10];
end
if ~isempty(yaxis_Limits)
    ylim([min(yaxis_Limits) max(yaxis_Limits)]);
end
set(gca,'XDir','reverse','XTick',XTick)
hold on

% Process the percentiles
pObj=struct('o',[],'n',[],'L','');
for i=1:length(prc)
    pObj(i).o = plot(HC_plt_x,Boot_plt(i,:),cs{i},'LineWidth',2);
    pObj(i).n = prc(i);
    pObj(i).L = {['CL',int2str(prc(i)),'%']};
end
pObj_t = struct2table(pObj);
pObj_t = sortrows(pObj_t,'n','descend'); % sort the table by 'n'

% Take percentiles
t1 = pObj_t(pObj_t.n>50,:); t1 = table2struct(t1); %above the mean
t2 = pObj_t(pObj_t.n<50,:); t2 = table2struct(t2); %below the mean

% Mean and Empirical
h1 = plot(HC_plt_x,Boot_mean_plt,'k-','LineWidth',2); % Mean
h2 = scatter(HC_emp.Hazard,HC_emp.Response,10,'g','filled','MarkerEdgeColor','k'); % Historical

% Title
if length(staID)==1
    title({'StormSim-SST ';['Station: ',staID{1}]},'FontSize',12);
else
    title({'StormSim-SST ';['Station: ',staID{1},' ',staID{2}]},'FontSize',12);
end

% Axes labels
if use_AEP
    xlabel('Annual Exceedance Probability','FontSize',12);
else
    xlabel('Annual Exceedance Frequency (yr^{-1})','FontSize',12);
end
ylabel({yaxis_Label},'FontSize',12);

% Legend
legend([[t1.o],h1,[t2.o],h2],{t1.L,'Mean',t2.L,'Empirical'},...
    'Location','southoutside','Orientation','horizontal','NumColumns',5,'FontSize',10);
hold off

% Filename
fname = [path_out,'SST_','HC_',staID{1},'.png'];

% Save figure
switch a
    case 0
        exportgraphics(gcf,fname,'Resolution',150)
    case 1
        saveas(gcf,fname,'png')
end
end