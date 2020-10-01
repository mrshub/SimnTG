# SimnTG

This is the repository for the MRSHub submission "SimnTG package" (submitted by [Ronald Ouwerkerk, NIDDK/NIH, Bethesda MD, USA](ouwerkerkr@niddk.nih.gov)).

## Guidance

1) Sim_UFATG_PRESS.m
Matlab script for use with FID-A-master.
Generates simulated spectra of triglycerides with unsaturated fatty acids as would be acquired with a TE series with PRESS localization. The fatty acid chain length can be specified and
 Chemical shifts and J coupling data were extracted from the following publications:
K. Schaumburg and H._ J. Bernstein Lipids Vol3 1968 p 193
Eleni Alexandri et al Molecules 2017, 22, 1663
Carmen Salinero et al. Molecules 2012, 17, 6716-6727
And
http://neu.ilps.org/wp-content/uploads/2019/07/Omega_3.pdf

2) sim_outAdd.m
Matlab script for use with FID-A-master.
This sums two simulation outputs as produced by Sim_UFATG_PRESS.m and any of the other simulation scripts in FID-A-master/simulationTools. This script can be used to mix various simulated fatty acid types in order to mimic a biological oil or fat.

3) stackedoutplots.m
Matlab script for use with FID-A-master.
This function plots a series of spectra as produced by Sim_UFATG_PRESS.m with a damping factor and line broadening.

Additional: TGspinSystems.mat
A spin system as is created within Sim_UFATG_PRESS.m

## Acknowledgement

If you would be so kind to acknowledge me as the source of this software when you use it for a publication it will be much appreciated.
