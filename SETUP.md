# SETUP NOTES

This project was cloned to a local repo on an NRCan laptop which does NOT possess admin rights. Lacking admin rights presents problems for qutip installation via pip, which needs to compile C++ code and in turn requires MS Build tools. MS Build tools require admin privileges.

Instead, you MUST use an Anaconda interpreter, and the conda package manager. This package manager will most likely need to be added to your PATH.