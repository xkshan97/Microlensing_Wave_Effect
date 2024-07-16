# Microlensing Field Diffraction Integration Code

This repository contains code for calculating the diffraction integration of microlensing fields. The implementation utilizes the component decomposition algorithm \cite{shan1} and the adaptive hierarchical algorithm.

## Files and Descriptions

1. **Micro_field_adaptive.cpp**: The main program for the adaptive algorithm.
2. **GetPsi_micro_field.cpp**: Code for generating the microlensing field time delay using the adaptive algorithm.
3. **SampleMethid/RejectAndAcceptSample.cpp**: Code for randomly generating microlenses based on the microlensing mass function.
4. **GetMicroDiffrac.py**: Program for calculating the time delay using the component decomposition algorithm.

## Example

To run a simple example, execute the following commands:

```sh
bash Example.sh
./Example

## References

\begin{thebibliography}{9}

\bibitem{shan1} Shan, X., et al. (2023). Wave effect of gravitational waves intersected with a microlens field: A new algorithm and supplementary study. Sci.China Phys.Mech.Astron., 66 (2023). https://doi.org/10.1007/s11433-022-1985-3
\bibitem{shan2} Shan, X., et al. (2024). Wave effect of gravitational waves intersected with a microlens field
II: An adaptive hierarchical tree algorithm and population study. Sci.China Phys.Mech.Astron., 66 (2023). https://doi.org/10.1007/s11433-022-1985-3

\end{thebibliography}
