 ### Gridsynth/NewSynth Adapted for the Lattice Surgery Compiler

Additions:

  * `LsqeccGlue/ExportCCals.hs` exposes an ABI to call gridsynth from C code
  * `include/grisynth_ccals.h` contains a function signature for this ABI for C and C++
  * `test.c`: contains an example of how to call gridsynth from C
  * `Makefile`: shows how to build and link gridsynth, using test.c as an example

Modifications:

  * `newsynth.cabal` adds the `foreign-library gridsynth_ccals`

For a more precise list of changes check the commit history.