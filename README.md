
# Usage
## Normal Use
### 1. Configure Run
For the first run change `maxdLdpdt.f90` to
<pre>
      standardParams = .false.

      if(standardParams) then
        ...
      else
         iguess = 50
         lebesgueQ = 4.0_pr
         resol = 1024
         bIterOffset = 0
         optimizationIterOffset = 0
      end if
      ...
</pre>
to continue a previous branching procedure
<pre>
      standardParams = .true.
      if(standardParams) then
         iguess = 9
         loadTempFunctionName = "u_result_q9_n1024_B022_iter00150.nc"
         call set_q_resol_bIterOffset_optimIterOffsets(loadTempFunctionName)
      else
        ...
      end if
</pre>

### <h3 id="compilingHeader">2. Compiling</h3>
`make`

### 3. Copy to Scratch
`cp prog ~/scratch/.../path/to/runFolder`

### 4. Configure Slurm File
create a `run.sh` file in `~/scratch/.../path/to/runFolder` with content like
```
#!/bin/bash
#SBATCH --account=rrg-aeinstein             # or def-aeinstein
#SBATCH --ntasks=128                        # best on nibi for 1024^3
#SBATCH --mem-per-cpu=12000M                # fine for 1024^3
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=my@email.address        # if you want to get notified about the run
#SBATCH --mail-type=begin,end,fail
srun prog
```
### 5. Submit
`sbatch run.sh`

## Special Usage
### Q-Continuation
For q continuation $q=4\to q=3.7$ configure `maxdLdpdt.f90`
<pre>
standardParams = .false.`
...
qContinuation = .true.`
...
scratchPath = "/home/MYUSERNAME/scratch/"
qContNcFileFolder = scratchPath//"8_1024_production/q4/5/output/ncFiles"
numberOfqValues = 2
qStartOffset = 0                 ! offset if wanna start for first file at different q. could be since to continue a sim that ended within a file and different q values
qContStart = 4.0_pr
qContEnd = 3.7_pr
bIterRangeStart = -1             ! specify range of used original b values, negative = all
bIterRangeEnd = -1
</pre>
then [continue like normal](compilingHeader)

### Special $B$ values
For q continuation $q=4\to q=3.7$ configure `maxdLdpdt.f90`
<pre>
standardParams = .false.`
...
!!! overwrite B values !!!
allocate( B_list_overwrite(1:17+bIterOffset) )
B_list_overwrite(1) = 1.0_pr
do B_list_iterator=2,size(B_list_overwrite)
    constraintB = B_list_overwrite(B_list_iterator-1)*10**(1.0_pr/8.0_pr)
    B_list_overwrite(B_list_iterator) = constraintB
end do
</pre>
then [continue like normal](compilingHeader)




            

# Fortran File Structure
```
maxLq/
├── maxdLpdt.f90          // main program: define vars, startup, ...
├── initialize.f90        // initializes parameters
├── optimization.f90      // main solver function: riemannian conjugate gradient
├── function_ops.f90      // math functions: R(u), nabla u, Delta^{-1} u, ...
├── global_variables.f90  // shared variables
├── fftwfunction.f90      // fftw3-mpi wrapper
├── data_ops.f90          // data import/export
```

# Algorithm Explanation

See also the [full algorithm on p. 23 of arXiv:2607.02739v1](https://arxiv.org/pdf/2607.02739v1#page=23) for more definitions.

Continuation approach used to compute branches of local maximizers $u_B^e$ for increasing values of the constraint parameter $B$.

<span style="color:teal">The corresponding Fortran variables, functions, explanations are colored</span>

## Input

- $q$ — <span style="color:teal">lebesgueQ (def in maxdLpdt.f90/main)</span> — Exponent of the Lebesgue space of interest
- $u_{\text{init}}$ — <span style="color:teal">uvec (set by function_ops.f90/initial_guess)</span> — Initial guess corresponding to the small-data limit 
- $B_{\text{init}}$ — <span style="color:teal"> B_list (maxdLpdt.f90/main)</span> — Initial value for the constraint parameter
- $B_{\max}$ — <span style="color:teal"> B_list (maxdLpdt.f90/main)</span> — Maximum value for the constraint parameter
- $\delta B$ — <span style="color:teal"> B_list (maxdLpdt.f90/main)</span> — Increment of the constraint parameter
- $\epsilon$ — <span style="color:teal">OPTIM_TOL (automatically initialize.f90/setStandardParams)</span> — Convergence tolerance
## Output

- $u_B^e$ — extreme states
- $R(u_B^e)$ — objective functional values
for $B_{\text{init}} \le B \le B_{\max}$

## Function `MainBranch()`

<pre>
B ← B_init
u0 ← u_init
repeat                                            <span style="color:teal">// maxdLpdt.f90/main (MAIN LOOP)</span>
    u_B^e ← SolveProblemHilbert(q, B, u0, ε)      <span style="color:teal">// maxdLpdt.f90/call maxdLqdt → optimization.f90/maxdLqdt </span>
    Evaluate R(u_B^e)                             <span style="color:teal">// calculated in optimization.f90/J1 = eval_Jcall </span>
                                                  <span style="color:teal">// saved in save_to_optimizationResultList → data_ops.f90/save_... </span>
    u0 ← u_B^e                                    <span style="color:grey">// Branch continuation with prev max</span>
    B ← B(1 + δB)                                 <span style="color:grey">// Increase constraint value</span>
until B ≥ B_max
</pre>

## Function `SolveProblemHilbert(q, B, u0, ε)` — <span style="color:teal">optimization.f90/maxdLqdt</span>

<pre>
k = 0
Γ_{τ_{-1} d_{-1}}(d_{-1}) = 0
repeat                                            <span style="color:teal">// optimization.f90/OPTIMIZATION LOOP</span>
    Compute the L^2 gradient ∇_L R(u_k)           <span style="color:teal">// optimization.f90/GradL2ForLq → function_ops.f90</span>
    Compute the H^{3/2 - 1/q} gradient ∇_H R(u_k) <span style="color:teal">// optimization.f90/HilbertGradient → function_ops.f90</span>
    Compute the step direction d_k using CG       <span style="color:teal">// optimization.f90/CONJUGATE GRADIENT</span>
    Find step size τ_k                            <span style="color:teal">// optimization.f90/mnbrak</span>
    Compute the vector transport Γ_{τ_k d_k}(d_k) <span style="color:teal">// optimization.f90/vectorTransport → function_ops.f90</span>
    Compute u_{k+1}                               <span style="color:teal">// optimization.f90/UPDATE U</span>
    k ← k + 1
until |R(u_{k+1}) − R(u_k)| / |R(u_k)| < ε

return u_{k+1}
</pre>

