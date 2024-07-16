# Microlensing Field Diffraction Integration Code

This repository contains code for calculating the diffraction integration of microlensing fields. The implementation utilizes the component decomposition algorithm and the adaptive hierarchical algorithm.

## Files and Descriptions

1. **Micro_field_adaptive.cpp**:
   The main program for the adaptive algorithm.
2. **GetPsi_micro_field.cpp**:
   Code for generating the microlensing field time delay using the adaptive algorithm.
3. **SampleMethid/RejectAndAcceptSample.cpp**:
   Code for randomly generating microlenses based on the microlensing mass function.
4. **GetMicroDiffrac.py**:
   Program for calculating the time delay using the component decomposition algorithm.

## Example

To run a simple example, execute the following commands:

bash Example.sh
./Example




## Citation

If you find our work useful in your research, please cite our paper:


@article{Shan:2022xfx,
    author = "Shan, Xikai and Li, Guoliang and Chen, Xuechun and Zheng, Wenwen and Zhao, Wen",
    title = "{Wave effect of gravitational waves intersected with a microlens field: A new algorithm and supplementary study}",
    eprint = "2208.13566",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1007/s11433-022-1985-3",
    journal = "Sci. China Phys. Mech. Astron.",
    volume = "66",
    number = "3",
    pages = "239511",
    year = "2023"
}

@article{Shan:2024abc,
    author = "Shan, Xikai and Li, Guoliang and Chen, Xuechun, Zhao, Wen, Hu, Bin and Mao Shude",
    title = "{Wave effect of gravitational waves intersected with a microlens field II: An adaptive hierarchical tree algorithm and population study}",
    eprint = "",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "",
    journal = "Sci. China Phys. Mech. Astron.",
    volume = "",
    number = "",
    pages = "",
    year = ""
}
