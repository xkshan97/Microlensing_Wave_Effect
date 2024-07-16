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
