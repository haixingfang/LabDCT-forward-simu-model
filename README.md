# [LabDCT forward simulation model](https://github.com/haixingfang/LabDCT-forward-simu-model)
Code for simulating diffraction patterns of laboratory X-ray diffraction contrast tomography (LabDCT).
Comparison between simulated and experimental projections are also available.The scripts may be continuously updated as work progresses. The code was developed by [Haixing Fang](https://orcid.org/0000-0001-8114-5276) in his postdoc project woking together with [Prof.dr. Dorte Juul Jensen](https://www.dtu.dk/english/service/phonebook/person?id=38577&tab=2&qt=dtupublicationquery) and [Dr. Yubin Zhang](https://www.dtu.dk/english/service/phonebook/person?id=50135&tab=2&qt=dtupublicationquery). The project is funded by the European Research Council (ERC) under the Europea Union's Horizon 2020 research and innovation programme (M4D - grant agreement 788567).

# Preparations for running the code
Installing external Matlab toolboxes that are freely accessible:
1. DIPimage: http://www.diplib.org/, for image processing.
2. Multi-Parametric Toolbox 3 (mpt3): https://www.mpt3.org/, for generating Voronoi cells and polyheron mesh in 3D.
3. MTEX toolbox: https://mtex-toolbox.github.io/, for analyzing and plotting crystallographic data.

For convenience, I have attached the DIPimage 2.9, mtex 5.1.1 and mpt3 scripts in the folder.
One can just clone the whole master file for running the code without the need to download these toolboxes.

It is also recommended to have 'computer vision system toolbox' installed with your own Matlab package.

All the codes have been tested executable with Matlab 2014b or above.
However, it is preferentially to run the code with a Matlab version 2018b or later.

# How to run the code
## step 1 - run [input_main.m](https://github.com/haixingfang/LabDCT-forward-simu-model/blob/master/input_main.m)
This step is to setup grain structure input, which can be either virtual rendered grain structure or from experimentally characterized structure written in h5 file. One can change values of variables 'grain_file' and 'grain_flag' to modify input option.

## step 2 - run [diffLabDCTsim_poly_3Dmesh_abs_DQE.m](https://github.com/haixingfang/LabDCT-forward-simu-model/blob/master/diffLabDCTsim_poly_3Dmesh_abs_DQE.m)
This requires input from step 1 to run the simulations of diffraction images.
Note to remember check the experimental parameter defined in [exp_parameters.m](https://github.com/haixingfang/LabDCT-forward-simu-model/blob/master/exp_parameters.m).
Images will be saved in the '\TFT\' folder and data saved in the '\DA\' folder.

If the meshed input structure already exists, one can skip the first step and directly run 'diffLabDCTsim_poly_3Dmesh_abs_DQE.m'.

# Others
If you wish to compare the simulated diffraction image with the experimental one as well as analyzing experimental LabDCT data.
You can run [diffLabDCTsim_poly_3Dmesh_comp_exp_v3_abs_DQE.m](https://github.com/haixingfang/LabDCT-forward-simu-model/blob/master/diffLabDCTsim_poly_3Dmesh_comp_exp_v3_abs_DQE.m) provided you have the meshed input of the experimental grain structure
and also the experimental diffraction images, which by default named as 'proj0000.tiff' alike.

The output images will be saved in the '\TFT_cmp\' folder and data saved in the '\DA_cmp\' folder.

Examples of input and one experimental LabDCT projection image can be found in the '\Examples\' folder.

# Remind
Always start with simulations for one projection at a certain rotation angle before running simulations for a whole dataset, e.g. 181 projections for a full rotation of 360 degrees.

# License
This package is free to use, ditribute and adapt for non-commercial use only.
See [LICENSE](https://github.com/haixingfang/LabDCT-forward-simu-model/blob/master/LICENSE) for license rights and limitations (CC BY-NC 4.0).

# Reference
[H. Fang, D. Juul Jensen, Y. Zhang, A flexible and standalone forward simulation model for laboratory X-ray diffraction contrast tomography, Acta Crystallographica Section A, 2020, vol.76.](https://doi.org/10.1107/S2053273320010852) <br>
Please cite this article if you use or get inspired by the code presented here.

## Contact via hfang@mek.dtu.dk or haixingfang868@gmail.com

