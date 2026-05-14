![SSLogo](https://github.com/Coastal-Hazards-System/StormSim-Library/assets/51959561/2b532547-1716-4bdb-815e-9e646b93615a)

# *CHESS – Coastal Hazards, Engineering, and Structures System*
### *A Modular Coastal Engineering Framework for Probabilistic Coastal Hazard Analysis and Structure Design*
### *(formerly StormSim)*
 
Thank you for being a part of this journey as we unveil **CHESS - the Coastal Hazards, Engineering, and Structures System**, developed by the the U.S. Army Corps of Engineers (USACE), Coastal Hazards Group (CHG). Your interest and participation are truly appreciated.
 
## A New Name for an Expanding System
 
After more than two decades of continuous development, the suite previously known as **StormSim** has grown far beyond storm simulation. Today it spans probabilistic coastal hazard analysis, life-cycle stochastic simulation, response-based structural design, and an expanding ecosystem of engineering analysis modules. To better reflect what the system actually is, and where it is going,it is being rebranded as **CHESS: Coastal Hazards, Engineering, and Structures System**.
 
The intent is for CHESS to serve the modern USACE coastal engineering community can rely on as a common reference for coastal hazard and structure analysis. CHESS is intended to be built on modern probabilistic frameworks, high-fidelity numerical modeling, and an open, modular software architecture.
 
**What is changing:**
- The umbrella name: **StormSim → CHESS**.
- Module prefixes: **StormSim-PST → CHESS-PST**; **StormSim-JPM → CHESS-JPM**; **StormSim-LCS → CHESS-LCS**; **StormSim-PROS → CHESS-PROS**; etc.
- Repository, documentation, and presentation branding will be updated on a phased schedule.
**What is *not* changing:**
- The development team, leadership, and Coastal Hazards Group (CHG) organizational home.
- The underlying methods, e.g., JPM, EST/PST, GPD-based extremal analysis, response-based design, life-cycle simulation, and their published technical foundations.
- The codebase, input file formats, and existing user workflows. Projects, scripts, and CHS-based analyses developed under the StormSim name continue to work without modification.
- The integration with the **Coastal Hazard System (CHS)** as the primary forcing-data backbone.
- The validity of prior peer-reviewed literature, technical reports, and study citations referencing StormSim. These remain authoritative references to the same underlying tool suite.
For a transition period, documentation, releases, and presentations will use the form **“CHESS (formerly StormSim)”** to preserve continuity for existing users and the published record. A separate FAQ document addresses citation guidance, version compatibility, and repository migration details.
 
## Development Stage
CHESS tools (under both the StormSim and CHESS names) have been in development for more than 20 years and have contributed to numerous USACE projects and federal studies. While the existing tools have undergone extensive testing, pauses in funding have inevitably delayed the large-scale release of the full suite.
 
We are genuinely excited to have you on board as part of this new phase. Together, we are paving the way for a future of innovative probabilistic analysis and engineering design workflows under a unified, agency-wide system.
 
## External Dependencies
### Coastal Hazard System (CHS) regional coastal storm hazard studies and PCHA results (https://chs.erdc.dren.mil)
The CHS is a national-scale, multiagency initiative for quantifying coastal storm hazards along U.S. coastlines and other strategic locations critical to national security. Key CHS components encompass a comprehensive and robust probabilistic framework (CHS-PF) and an online database (CHS-BD) for production and deployment of PCHA results. CHS directly supports USACE’s Coastal Storm Risk Management (CSRM), Flood Risk Management (FRM), and Navigation missions, in addition to other federal, state, and local government agencies. The CHS ensures that the atmospheric and hydrodynamic forcing data, including mean sea level fluctuations and astronomical tides, used to evaluate multiple aspects related to USACE missions are based on consistent methodology, quality, and are uniquely tailored to address the primary drivers of coastal hazards specific to the regions in which assets are located, including hurricanes or tropical cyclones (TCs), extratropical cyclones (XCs), and other extreme hydrodynamic events.
 
CHS is a unique data resource that spans probability space, was developed with high-fidelity modeling, is spatially uniform and dense on a national scale, and includes aleatory and epistemic uncertainty information necessary for determining probabilistic responses. Modeling results are reported through virtual gauges (e.g., save point or node) across regional studies. These virtual gauges include modeling results for both ocean circulation (ADCIRC) and wave models (STWAVE/WAM/SWAN). CHS regional study probabilistic and numerical modeling results can be downloaded [here](https://chs.erdc.dren.mil).
 
![image](https://github.com/Coastal-Hazards-System/StormSim-Library/assets/51959561/df664577-5593-4410-99c1-b967f68ba59f)
 
The following figure summarizes the CHS data download process: (1) choose a regional study; (2) apply the modeling condition filter; (3) choose the most suitable virtual gauge for your project.
When selecting project-forcing data, consider the hazards you want to include (XC vs. TC vs. both) and the level of fidelity (Peaks vs. Timeseries vs. both). Currently, CHESS only supports base modeling conditions as a forcing input. Utilizing CHS data means any analysis will involve the full storm suite for a given regional study. The downloaded zip folder must contain Peaks and/or Timeseries files for both models associated with a savepoint (ADCIRC + Wave Model).
 
<p align="center" width="100%">
    <img width="100%" src="https://github.com/Coastal-Hazards-System/StormSim-Library/assets/51959561/f8f5a12f-3e75-44bd-a6b5-e7af0845a877">
</p>
For a more detailed guide on accessing and navigating CHS files/website, [visit this link](https://chs.erdc.dren.mil/Content/documents/CHSQuickGuide.pdf).
 
# Select List of CHESS Modules in Development
Some of the CHESS modules (formerly StormSim modules) currently in development and that will be hosted in this repository include:
 
### CHESS-PST: Probabilistic Simulation Technique *(formerly StormSim-PST)*
CHESS-PST estimates extra-tropical cyclone (XC) storm response hazards. CHESS is a component of the Coastal Hazards System (CHS https://chs.erdc.dren.mil, Nadal-Caraballo et al. 2020), a national-scale initiative to quantify coastal storm hazards. Results from high-fidelity, physics-informed numerical modeling of coastal storm events spanning the practical probability space for the U.S. coastline are stored on an online database (CHS-DB). Historical XCs were modeled in a high-fidelity coupled hydrodynamic framework (Massey et al. 2012) for regional studies to produce coastal wave and water level responses. This was completed for multiple sea level rise conditions. Both peak and timeseries values are available. CHESS can evaluate responses over entire timeseries within a probabilistic framework with relatively small computational costs. Tides can be either randomly sampled, applied as uncertainty, or applied as a skew tide to account for nonlinearities. CHESS-PST ingests modeled historical storm responses. Responses are bootstrapped sampled, including aleatory uncertainty, to encompass multiple sequences or life cycles of storm responses. Peaks-over-threshold is used to identify extreme events from sampling. These extreme events are fit to a generalized Pareto distribution (GPD), which captures the low-frequency tail of response hazards.
 
### CHESS-JPM: Joint Probability Method *(formerly StormSim-JPM)*
CHESS-JPM estimates tropical cyclone (TC) storm response hazards (Nadal-Caraballo and Melby 2014). CHESS is a component of the Coastal Hazards System (CHS https://chs.erdc.dren.mil/, Nadal-Caraballo et al. 2020), a national-scale initiative to quantify coastal storm hazards. Results from high-fidelity, physics-informed numerical modeling of coastal storm events spanning the practical probability space for the U.S. coastline are stored on an online database (CHS-DB). Synthetic TCs were modeled in a high-fidelity coupled hydrodynamic framework (Massey et al. 2012) for regional studies to produce coastal wave and water level responses. This was completed for multiple sea level rise conditions. Discrete storm weights (DSW) define TC storm probability and are used to estimate hazards. Both peak and timeseries values are available. CHESS can evaluate responses over entire timeseries within a probabilistic framework with relatively small computational costs. Tides can be either randomly sampled, applied as uncertainty, or applied as a skew tide to account for nonlinearities. TC responses, such as waves and water levels, are ingested into CHESS-JPM, aleatory uncertainty is applied to these storm responses using hundreds of thousands of storms, and associated DSWs are used to estimate best-estimate (BE) exceedance probabilities. Confidence limits are computed by applying epistemic uncertainties to BE hazard curves. Additional details about coastal storm uncertainties are in Gonzalez et al. (2019).
 
### CHESS-LCS: Life Cycle Simulation *(formerly StormSim-LCS)*
CHESS-LCS computes time-dependent stochastic analysis of coastal structure responses. Current coastal structure responses include dune and beach morphology and rubble mound armor stone damage progression. CHESS is a component of the Coastal Hazards System (CHS https://chs.erdc.dren.mil, Nadal-Caraballo et al. 2020), a national-scale initiative to quantify coastal storm hazards. Results from high-fidelity, physics-informed numerical modeling of coastal storm events spanning the practical probability space for the U.S. coastline are stored on an online database (CHS-DB). Synthetic TCs and XCs were modeled in a high-fidelity coupled hydrodynamic framework (Massey et al. 2012) for regional studies to produce coastal wave and water level responses. This was completed for multiple sea level rise conditions. Discrete storm weights (DSW) define TC storm probability and are used to estimate hazards. Both peak and timeseries values are available. CHESS can evaluate responses over entire timeseries within a probabilistic framework with relatively small computational costs. Tides can be either randomly sampled, applied as uncertainty, or applied as a skew tide to account for nonlinearities. Life cycles are generated from sampled storms, sampled using a Poisson distribution to create a life cycle reflective of the larger population of storm intensities and associated probabilities at the study location. TCs are sampled using DSWs and XCs are sampled either historically or from stochastic simulation using bivariate or multivariate Gaussian copulas. Storms may also be sampled by binning storms by intensity prior to sampling, then sampling from intensity bins to represent storm intensity probabilities more accurately throughout the entire life cycle. Aleatory uncertainties are accounted by stochastically creating hundreds of life cycles, which also ensure statistical convergence. Structure life-cycle responses can be used to estimate statistical structure reliability and performance.
 
### CHESS-PROS: Probabilistic Responses Of Structures *(formerly StormSim-PROS)*
CHESS-PROS probabilistically estimates structure response hazards using response-based methods. Coastal structure responses include floodwall, levee, and rubble mound overtopping and volume discharge, levee and rubble mound runup, floodwall hydrodynamic and hydrostatic pressures, and rubble mound stone stability. CHESS is a component of the Coastal Hazards System (CHS https://chs.erdc.dren.mil, Nadal-Caraballo et al. 2020), a national-scale initiative to quantify coastal storm hazards. Results from high-fidelity, physics-informed numerical modeling of coastal storm events spanning the practical probability space for the U.S. coastline are stored on an online database (CHS-DB). Synthetic TCs and XCs were modeled in a high-fidelity coupled hydrodynamic framework (Massey et al. 2012) for regional studies to produce coastal wave and water level responses. This was completed for multiple sea level rise conditions. Discrete storm weights (DSW) define TC storm probability and are used to estimate hazards. Both peak and timeseries values are available. CHESS can evaluate responses over entire timeseries within a probabilistic framework with relatively small computational costs. Tides can be either randomly sampled, applied as uncertainty, or applied as a skew tide to account for nonlinearities. Structure responses are computed for hundreds of thousands of storm realizations following the CHESS-JPM and CHESS-PST workflows. Exceedance probabilities of the structure responses themselves are estimated. CHESS-PROS maintains the high fidelity multivariate statistical and physical interdependencies between storm forcing parameters without assuming identical probabilities of forcing parameters and structure responses (Stehno and Melby, in review). Epistemic uncertainties associated with storm response and empirical equations are applied as confidence limits to the BE response hazards. Results from CHESS-PROS are used for coastal structure designs when designs require non-exceedance of a structure response hazard at a given probability.
 
# Citing CHESS and Legacy StormSim Work
 
Prior peer-reviewed literature, ERDC technical reports, and study deliverables that cite **StormSim** or any **StormSim-\*** module remain valid references to the same underlying tool suite and methods. No revisions to existing citations are required.
 
For new work, the recommended citation form is:
 
> *Coastal Hazards, Engineering, and Structures System (CHESS), formerly StormSim. U.S. Army Engineer Research and Development Center, Coastal and Hydraulics Laboratory, Coastal Hazards Group, Vicksburg, MS.*
 
Module-level citations should use the new prefix (e.g., **CHESS-JPM**) and may optionally include the legacy designation (e.g., *“CHESS-JPM (formerly StormSim-JPM)”*) during the transition period.
 
# Contact Information

## CHESS Feedback & Support
*Development Team:* CHESS@usace.army.mil (**email address is currently inactive**)

## CHESS Team

Fabian A. Garcia-Moreno - Fabian.A.Garcia-Moreno@usace.army.mil

Kevin C. Hodgens - Kevin.C.Hodgens@usace.army.mil

Abigail L. Stehno, PhD - Abigail.L.Stehno@usace.army.mil

Norberto C. Nadal-Caraballo, PhD - Norberto.C.Nadal-Caraballo@usace.army.mil

Madison C. Yawn - Madison.C.Yawn@usace.army.mil


# References
Gonzalez, V.M., Nadal-Caraballo, N.C., Melby, J.A., Cialone, M.A. (2019). Quantification of Uncertainty in Probabilistic Storm Surge Models: Literature Review. Technical Report ERDC/CHL TR-21-15, Vicksburg, MS: U.S. Army Engineer Research and Development Center.

Melby, J.A., Nadal, N.C., Males, R.M. (2011). CSsim: Breakwater-Harbor Time-Dependent Life-Cycle Analysis Software. Coastal Structure 2011 Conference Proceedings, World Scientific, Singapore, pp 649-658.

Melby, J.A., Massey, T.C., Das, H.S., Nadal-Caraballo, N.C., Gonzalez, V.M., Bryant, M.A., Tritinger, A.S., Provost, L.A., Owensby, M.B, and Stehno, A.L. (2021). Coastal Texas Protection and Restoration Feasibility Study: Coastal Texas flood risk assessment: hydrodynamic response and beach morphology. ERDC TR-21-11, Vicksburg, MS: U.S. Army Engineer Research and Development Center.

Melby, J.A., Massey, T.C., Stehno, A.L., Nadal-Caraballo, N.C., Misra, S., Gonzalez, V.M. (2021). Sabine Pass to Galveston Bay, TX Pre-Construction, Engineering and Design PED Hurricane Coastal Storm Surge and Wave Hazard Assessment: Report 1- Background and Approach. Technical Report ERDC/CHL TR-21-15, Vicksburg, MS: U.S. Army Engineer Research and Development Center. http://dx.doi.org/10.21079/11681/41820

Nadal-Caraballo, N. C., J. A. Melby, V. M. Gonzalez, A. T. Cox (2015). North Atlantic Coast Comprehensive Study (NACCS) Coastal Storm Hazards from Virginia to Maine. ERDC/CHL TR-15-5. U.S. Army Engineer R&D Center, Vicksburg, M.S

Nadal-Caraballo, Campbell, Gonzalez, Torres, Melby, Taflanidis, (2020). Coastal Hazards System: A Probabilistic Coastal Hazard Analysis Framework. In: Malvárez, G. and Navas, F. eds., Global Coastal Issues of 2020. Journal of Coastal Research, Special Issue No. 95, pp. 1211-1216.

Stehno, A.L., (2021). Coastal Structure Overtopping and Overflow Stochastic Simulation Method Comparison. M.S. Thesis. Mississippi College, Clinton, Mississippi. 10.13140/RG.2.2.10123.41761 
