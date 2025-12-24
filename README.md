# SVO

DESCRIPTION
-----------
These are Julia programs (version 1.5.3) and Matlab programs (version 2023b) designed for the production of a paper entitled "Evolutionary dynamics of social value orientation."

FILE
-----
<pre>
Figure2.mlx     obtain the frequency of cooperation for a range of benefit value b and SVO \psi (see Figure 2 in the main text)
Figure3.mlx     obtain the variance in node degrees for six types of networks, frequency of cooperation for a range of structural coefficients (see Figure 3 in the main text)
Figure4.jl      obtain the fixation probability of SVO psi_critical in a population of SVO psi_instance, and the fixation probability of SVO psi_instance in a population of SVO psi_critical (see Figure 4cf                  in the main text)
Figure5.jl      obtain the cooperation frequencies in the direct evolutionary processes and in the indirect processes (see Figure 5ef in the main text)
</pre>
                        
INSTALLATION
------------
<pre>
*) Download the open-source software Julia 1.5.3 or a newer version from https://julialang.org/downloads/.
*) Install Julia on your local PC. The expected installation time is approximately five minutes.
*) Download the open-source software VS Code 1.77.3 or a newer version from https://code.visualstudio.com/download.
*) Install VS Code on your local PC. The expected installation time is approximately ten minutes.
*) Install the Julia VS Code extension from https://code.visualstudio.com/docs/languages/julia. The expected installation time is approximately five minutes.
*) Run VS Code, add Figure4.jl, and install a solver package: 
   1) Type the following command in TEMINAL, and then press ENTER:
      julia
   2) Type the following command, and then press ENTER:  
      using Pkg; Pkg.add("IterativeSolvers");
*) Click "Run and Debug" to execute the program.
*) The expected outputs are two fixation probabilities for two examples, i.e., taking b=10, rho_critical for psi_critical taking up a population of -\pi/2, and rho_instance for -pi/2 taking up a population of psi_critical. The expected run time is approximately 10 hours.
*) To reproduce the entire panel of Figure 4cf, sample 20 psi_instance values from (-pi, pi) and then run the codes above. The time can be reduced by distributing the computations across multiple computers.
  
*) Add Figure5.jl, and click "Run and Debug" to execute the program. The expected outputs are two cooperation frequency, one for the direct evolutionary process (xC_direct_average) and the other for the indirect evolutionary process (xC_indirect_average). The expected run time is approximately 1 minute.
*) To reproduce the entire panel of Figure 5ef, sample 20 benefit values from (1, 20) and then run the codes above. The time can be reduced by distributing the computations across multiple computers.
  
*) Run Matlab, add Figure2.mlx and click "Run" to execute the program. The expected outputs are
   a figure showing the frequency of cooperation for a range of benefit value b and SVO \psi (see Figure 2 in the main text). The expected run time is approximately 3 minutes.

*) Run Matlab, add Figure3.mlx and click "Run" to execute the program. The expected outputs are
   a figure showing the variance in node degrees for six types of networks, frequency of cooperation for a range of structural coefficients (see Figure 3 in the main text). The expected run time is approximately 1 minutes.
</pre>

QUESTIONS
---------
For any questions about this program, please contact
Dr. Qi Su, Email: qisu@sjtu.edu.cn.
