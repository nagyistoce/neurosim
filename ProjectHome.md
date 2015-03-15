**Features and Stats**
  * Neuron model: Izhikevich
  * Model equation solver: adaptive Parker-Sochacki (PS) method with Newton-Raphson (NR) root search.
  * Solver precision:
    * PS: up to full floating point standard.
    * NR: up to 1E-16.
  * Floating point standard: single precision
  * Simulation step size: variable (each spike is defined according to its unique time).
  * Network size range: up to 2Mil
  * Synaptic connection size range: up to 16K
  * Performance tuning: characterization script allows to search for application configuration that provides the best performance
  * Programming language: OpenCL, C++
  * Implementation type: parallel and sequential
  * Supported device architectures:
    * GPU: AMD Cayman, AMD Tahiti
    * CPU: any x86\_64
  * OS: Windows 7, 64 bit
  * Required software:
    * Config 1: AMD APP SDK, Visual Studio Pro
    * Config 2: AMD APP SDK, VC C++ compiler, Strawberry Perl, msys
  * Verification: full results verification between simulation on CPU and GPU
  * GPU/CPU execution time speedup factor:
    * Tahiti GPU vs single thread on AMD Phenom II, 3.2 GHz CPU:
| **Network Size (neurons)** |	| **Average Synapses per Neuron** |	| **Average Events per Iteration** |	| **Average Spikes per Iteration** |	| **Total synapse count** |	| **GPU Time per Step, (ms)** |	| **CPU Time per Step, (ms)** |	| **Speedup Factor** |
|:---------------------------|:|:--------------------------------|:|:---------------------------------|:|:---------------------------------|:|:------------------------|:|:----------------------------|:|:----------------------------|:|:-------------------|
| 2^21 |	 | 90 |	 | 227152 |	 | 2522 |	 | 190190561 |	 | 13.5 |	 | 659 |	 | 48 |
| 2^17|	 | 1458 |	 | 373174 |	 | 257 |	 | 191209680 |	 | 5.7 |	 | 279 |	 | 48 |
| 2^14|	 | 11677 |	 | 296469 |	 | 25 |	 | 191323616 |	 | 3.2 |	 | 283 |	 | 88 |

<br>
<br>
<br>

<b>Literature</b>
<br>
<a href='http://code.google.com/p/neurosim/source/browse/wiki/Scalable_Multi-Precision_Simulation_of_Spiking_Neural_Networks_on_GPU_with_OpenCL.pdf'>Scalable Multi-Precision Simulation of Spiking Neural Networks on GPU with OpenCL</a>

<br>
<br>
<br>

<b>Current work in progress</b>
<ul><li>Code refactoring: object-oriented and easy to understand<br>
</li><li>Parallel kernel execution.<br>
</li><li>STDP feature.<br>
</li><li>Non-divergent root-search method as an alternative to NR.