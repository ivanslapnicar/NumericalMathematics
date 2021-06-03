# NumericalMathematics

Lecture notebooks for the course _Numerical Analysis_ or _Numerical Mathematics_ as a first (introductory) one-semester course which is standardly thought to STEM students.

Notebooks are used in the course
  _[Numerical Analysis](https://nastava.fesb.unist.hr/nastava/predmeti/8183)_ to Master's students of Computer Science at [FESB](https://www.fesb.unist.hr/).
  Notebooks are jointly written by Jesse Barlow from the Penn State University and Ivan Slapničar from University of Split.  

Notebooks are particularly useful in on-line teaching. Notebooks are written in [Julia](https://julialang.org) using [Pluto.jl](https://github.com/fonsp/Pluto.jl).

## Viewing the notebooks

You can view the notebooks at [https://ivanslapnicar.github.io/NumericalMathematics/](https://ivanslapnicar.github.io/NumericalMathematics/)

## Running the notebooks

You can run the notebooks in two ways:

### Running on `binder`

1. Go to [https://ivanslapnicar.github.io/NumericalMathematics/](https://ivanslapnicar.github.io/NumericalMathematics/) and choose the desired notebook.
2. Press `Edit or run this notebook` button and choose `binder`. This will read all the necessary packages and start the notebook (within several minutes).

### Running on your comoputer

1. Clone the entire repository using `git` command:
```
git clone https://github.com/ivanslapnicar/Matematika.git
```
If you are unfamiliar with the `git` tool, check GitHub [help pages](https://help.github.com/articles/set-up-git/). You can also download the repository as a zip file.

2. Install [Julia](https://julialang.org/downloads/). In Julia terminal run the commands
```
> using Pkg
> Pkg.add("Pluto")
> using Pluto
> Pluto.run()
```
This opens local Pluto server in your browser. Now you can choose the notebook and run it
(the notebboks are located in the directory `NumericalMathematics/Lectures/`).
