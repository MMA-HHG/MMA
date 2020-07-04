## The plot of the code:

1. See the comment at the end of the file for deatils about the MPI-scheduler.

2. The main program (MP) reads all the parameters from an input HDF5-archive (each process independently, read-only).

3. MP decides based on parameters the type of input field. Fist implementation: stick to numerical fields from CUPRAD stored in the archive.

4. MPI-scheduler executes a simulation in every point. THe data are directly written into the output using the mutex.

5. The code finishes.

- DEVELOPMENT: wrap the call of singleTDSE, where propagation is called to test, erase this extra step after


## Extensions/features already presented in 1DTDSE
THere is a list of features we added in the code throughout the time, we don't have them in the sheduler. But it would be worthy to re-introduce them simply by modifying the output structure.

- The features already presented in TDSE:
	1. Print wavefunction (WFT)
	2. Print Gabor transformation
	3. Photoelectron spectrum (PeS, Fabrice)
  4. Various outputs modifications (e.g. ionisation filtering done for Sophie)

Gabor is computationally demanding; WF and PeS are both data-storage demandig. It should be used then only in few prescribed points.


--There is possibility of various inputs:
	1. Numerical/analytic field
	2. Computation on velocity/length gauge

We have already an analytic model of a beam (Python/MATLAB), we will then rewrite it and construct the parameters on-the-fly.
The versatility in numeric-or-analytic field length-or-velocity gauge is ideal for testing of numerical vector potential that we can use after in SFA.


## Development notes:
<pre>
1) We can use checks whether parameters exist in the input HDF5-archive. If not, we can create them with default values.
Implementation: since this is I/O operation with one file, we need r/w. Maybe read paramc only by proc 0 and the broadcast structure (see the MPI book for transfering structs).

For reading, it should be easy. ** R/W may occur simultaneously in in the MPI loop. Separate I/O at the instant or ensure it will work (R/W from independent datasets may be fine???).
https://support.hdfgroup.org/HDF5/Tutor/selectsimple.html


2) we get rid of mutexes and use rather parallel acces to files all the time.
2.develop) it seems that many-readers many-writers would be possible by HDF5 parallel since we will not modify the file much. However, we may also try stick with independent files and eventually 
https://stackoverflow.com/questions/49851046/merge-all-h5-files-using-h5py
https://portal.hdfgroup.org/display/HDF5/Collective+Calling+Requirements+in+Parallel+HDF5+Applications
</pre>
