# NumericalMathematics

Notebooks for the course _Numerical Analysis_ or _Numerical Mathematics_ as a first (introductory) one-semester course which is standardly thought to STEM students.

Notebooks are used in the course
  _[Numerical Analysis](https://nastava.fesb.unist.hr/nastava/predmeti/8183)_ to Master's students of Computer Science at [FESBu](https://www.fesb.unist.hr/).
  Notebooks are jointly written by Jesse Barlow from the Penn State University and Ivan Slapničar from University of Split.  

Notebooks are particularly useful in on-line teaching. Notebooks are written in [Julia](https://julialang.org).

## Usage

Materials are written as [Jupyter](http://jupyter.org/) notebooks and/or [Pluto](https://github.com/fonsp/Pluto.jl) notebooks.
Notebooks can be used in the following ways:

### In Web Browser
You can view notebooks using following links:
* [Jupyter notebook viewer](http://nbviewer.ipython.org/url/github.com/ivanslapnicar/NumericalMathematics/tree/master/src/)
* [Pluto](https://ivanslapnicar.github.io/NumerickaMathematics/)

###  Local Installation and Running
* Download the notebooks (repository) using `git` command:
```
git clone https://github.com/ivanslapnicar/NumericalMathematics.git
```
If you are unfamiliar with `git` tool, you can taka a look at GitHub's [help pages](https://help.github.com/articles/set-up-git/) or you can download repository directly as zip file.
* Intall [Julia](https://julialang.org/downloads/). In Julia terminal run the commands
```
> using Pkg
> Pkg.add("IJulia")
> Pkg.add("Pluto")
```
The above commands need to be executed only once.
* Jupyter notebook server (for notebooks with extension `.ipynb`) is started with commands
```
> using IJulia
> notebook(detached=true)
```
and Pluto notebook server (for notebooks with extension `.jl`) is started with commands
```
> using Pluto
> Pluto.run()
```

You can now run notebooks which are located in the directory `NumericalMathematics/Lectures`
