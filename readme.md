
A survey and an extensive evaluation of popular audio declipping methods
========================================================================
[Pavel Záviška](https://orcid.org/0000-0003-2221-2058), [Pavel Rajmic](https://orcid.org/0000-0002-8381-4442), [Alexey Ozerov](https://orcid.org/0000-0002-7602-4610) and [Lucas Rencker](https://orcid.org/0000-0002-9332-6602)
------------------------------------------------------------------------


This readme file describes the MATLAB toolbox accompanying the article from the title.


### Requirements
The code has been developed in MATLAB version R2019b and has been tested in R2019b and R2020a. 
Several methods rely on the LTFAT toolbox (version 2.4.0 used) and Constrained OMP (C-OMP) relies on the CVX toolbox (version 2.2).

### Quick Tutorial
To use this declipping toolbox, just download all the files, add them to the MATLAB path and make sure that the LTFAT toolbox and CVX toolbox are properly installed.

The toolbox is organized into several subfolders:
   - "Evaluation_algorithms" contains the PEAQ method used in the paper to evaluate the objective quality of restoration. Unfortunately, PEMO-Q (also used in the paper) cannot be made part of the repository due to copyright issues.
   - "Methods" contains all the algorithms used in the extensive evaluation except for the NMF-based method (not publicly available).
   - "Numerical_results" contains ".mat" files with numerical results based on five different metrics (dSDR computed on the whole signal&#8212;dSDR_all, dSDR computed on the clipped samples only&#8212;dSDR_clipped, PEAQ, PEMO-Q, and Rnonlin). Suffix "RR" denotes replacing the reliable samples of the restored signals before evaluation. The specific values are arranged in a matrix, where rows represent individual audio excerpts, and columns represent input SDR values.    
   - "Sounds" contains wav-files used for testing and a mat file with the testing signals "Sound_database.mat"
   - "Tools" contains auxiliary methods for clipping the signal and computing the SDR.

The toolbox is designed to be run from the root folder of the toolbox. 
Each method in the "Methods" folder contains a "main" file that runs the declipping experiment for the selected audio signal and a clipping threshold (defined by the input SDR). It is possible to select an audio file from the preset or you can use your audio file. Furthermore, it is possible to select the input SDR and other parameters corresponding to the used algorithm and signal representation.

One main file could run one or more algorithms. The algorithm is chosen with the `algorithm` variable.
This setting applies to:
   - main_l1_relaxation.m
      * in `param.algorithm` select `'DR'` for the Douglas-Rachford algorithm (synthesis model) or `'CP'` for the Chambolle-Pock algorithm (analysis model);
      * in `param.reweighting` select "1" for reweighting of the coefficients (Rl1CC algorithms);
      * in `param.weighting` select "1" for parabola weighting of the coefficients
   - main_csl1.m;
      * in `param.algorithm` select `'csl1'` for the non-weighted variant, `'pcsl1'` for the psychoacoustically weighted variant and `'pwcsl1'` for the parabola-weighted variant;
   - main_Social_Sparsity.m
      * in `shrinkage` select `'EW'` for Empirical Wiener and `'PEW'` for Persistent Empirical Wiener 
   - main_SPADE.m
      * in `param.algorithm` select `'aspade'` for the analysis variant of SPADE and `'sspade_new'` for the (proper) synthesis variant of SPADE

The other main files (with SQAM suffix) can be run to reproduce the results from the paper&#8212;it declipps all ten audio files with the seven input distortion levels used (70 combinations in total). Note that such computation can be time consuming.

Default parameter values are identical to the values used for the experiments in the paper. For different sounds the optimum parameters could differ.

### Special Thanks
The authors would like to thank:
   - Lucas Rencker for providing the implementation of the Dictionary Learning method and for running the computations of the Dictionary Learning and the Constrained OMP.
   - Valentin Emiya for providing the implementation of the Constrained OMP and the Janssen method.
   - Matthieu Kowalski for providing the implementation of the Social Sparsity Declipper.
   - Ondřej Mokrý for running the computations of the Janssen method.

### Audio Excerpts
If you only want to listen to the audio excerpts, you can visit the web page:
https://rajmic.github.io/declipping2020/

### How to cite this toolbox
Please cite the following paper:

P. Záviška, P. Rajmic, A. Ozerov, and L. Rencker:
A survey and an extensive evaluation of popular audio declipping methods.
Submitted to IEEE Journal of Selected Topics in Signal Processing (2020).
Available at https://arxiv.org/abs/2007.07663.

### License
The code of this toolbox is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).

### References
- A. Adler, V. Emiya, M. Jafari, M. Elad, R. Gribonval, and M. Plumbley,
“A constrained matching pursuit approach to audio declipping,” 
in *2011 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)*, 
May 2011, pp. 329&#8211;332.

- S. Kitić, N. Bertin, and R. Gribonval, 
“Sparsity and cosparsity for audio declipping: a flexible non-convex approach,” 
in *12th International Conference on Latent Variable Analysis and Signal Separation (LVA/ICA)*, 
Aug. 2015, pp. 243&#8211;250.

- P. Záviška, P. Rajmic, O. Mokrý, and Z. Průša, 
“A proper version of synthesis-based sparse audio declipper,”
in *2019 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)*,
May 2019, pp. 591&#8211;595.

- P. Rajmic, P. Záviška, V. Veselý, and O. Mokrý, 
“A new generalized projection and its application to acceleration of audio declipping,”
*Axioms*, vol. 8, no. 3, 2019.

- A. J. Weinstein and M. B. Wakin, 
“Recovering a clipped signal in sparseland,” 
*Sampling Theory in Signal and Image Processing*, 
vol. 12, no. 1, pp. 55&#8211;69, 2013.

- K. Siedenburg, M. Kowalski, and M. Dorfler, 
“Audio declipping with social sparsity,” 
in *2014 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)*, 
May 2014, pp. 1577&#8211;1581.

- B. Defraene, N. Mansour, S. D. Hertogh, T. van Waterschoot, M. Diehl, and M. Moonen, 
“Declipping of audio signals using perceptual compressed sensing,” 
*IEEE Transactions on Audio, Speech, and Language Processing*, 
vol. 21, no. 12, pp. 2627&#8211;2637, 2013.

- P. Záviška, P. Rajmic, and J. Schimmel, 
“Psychoacoustically motivated audio declipping based on weighted l1 minimization,” 
in *2019 42nd International Conference on Telecommunications and Signal Processing (TSP)*, 
July 2019, pp. 338&#8211;342.

- L. Rencker, F. Bach, W. Wang, and M. D. Plumbley, 
“Consistent dictionary learning for signal declipping,” 
in *14th International Conference on Latent Variable Analysis and Signal Separation (LVA/ICA)*. 
July 2018, pp. 446&#8211;455.

- A. J. E. M. Janssen, R. N. J. Veldhuis, and L. B. Vries, 
“Adaptive interpolation of discrete-time signals that can be modeled as autoregressive processes,” 
*IEEE Transactions on Acoustics, Speech, and Signal Processing*, 
vol. 34, no. 2, pp. 317&#8211;330, 1986.


--------------------------------------------------
Pavel Záviška, Brno University of Technology, 2020
