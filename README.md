### Rotation decomposer for the LatticeSurgery compiler

Generate strings of quantum gates to approximate arbitrary phase rotations on a single qubit.

 - The directory `newsynth` contains [Selinger's Newsinth](https://www.mathstat.dal.ca/~selinger/newsynth), adapted to build on our systems. 

 - `gridsynth` is statically linked binary of newsynth compiled to run on linux machines. We were having toruble with the one provided on the authors' website

 - `gen_cached_rotations.py`produces a table of angles of the form `pi/2^n` and sequence of gates to achieve them.
