# Contributing to Fastscapelib

The Fastscapelib project is open source and welcomes your feedback and/or
contributions! There are many ways to contribute to the development of the
library, including bug reports, bug fixes, documentation improvements,
enhancement suggestions, and other ideas.

This short guide will help you getting started. Please check online resources
and tutorials for more details about the developer tools that are used here.

## Bug Reports, Suggestions and Discussions

Fastscapelib's source repository is hosted on
[GitHub](https://github.com/fastscape-lem/fastscapelib). This is also where you
can submit bug reports, enhancement ideas or general usage questions either in
the [issue tracker](https://github.com/fastscape-lem/fastscapelib/issues) or in
the [discussions](https://github.com/fastscape-lem/fastscapelib/discussions).

:::{tip}

Check this [stack-overflow
article](https://stackoverflow.com/help/minimal-reproducible-example) and this
[blog post](https://matthewrocklin.com/minimal-bug-reports) by Matthew Rocklin
for tips on how to write good bug reports.

GitHub supports the [markdown
syntax](https://docs.github.com/en/get-started/writing-on-github) for formatting
code, add links, etc. This is helpful for writing clear reports.

:::

## Fastscapelib's Git Repository

:::{note}

Contributing to the code or documentation requires some basic knowledge of
[Git](https://git-scm.com/). You can find many resources and tutorials online
for learning Git.

:::

Follow the steps below to get or update Fastscapelib's source code and prepare
your contribution (you can skip steps 1-2-3 if this is not your first
contribution):

1. Go on the Fastscapelib
   [GitHub](https://github.com/fastscape-lem/fastscapelib) page and click on the
   "Fork" button (you need to have an account on GitHub and be logged in).

2. Clone your Fastscapelib forked Git repository (by default, this will setup a
   remote ``origin`` repository pointing to your fork on GitHub):

```bash
$ git clone https://github.com/your-username/fastscapelib
$ cd fastscapelib
```

3. Add a new remote ``upstream`` repository pointing to the Fastscapelib main
   repository on GitHub:

```bash
$ git remote add upstream https://github.com/fastscape-lem/fastscapelib
```

4. Pull the latest updates from the Fastscapelib original repository into your
   local main branch (you can skip this step if this is your first
   contribution):

```bash
$ git checkout main
$ git pull upstream main
```

5. Create a new Git branch in which you will implement your bug fixes,
   documentation enhancement or new feature:

```bash
$ git checkout -b your-branch-name
```

(contributing-development-environment)=
## Creating a Development Environment

We strongly recommend installing all the dependencies and development tools in a
new [conda](https://conda.io/docs/) environment created specifically for this
project.

You can create a new conda environment with all the dependencies needed for both
C++ and Python development by running the following command from the
repository's root directory:

```bash
$ conda env create -n fastscapelib-dev -f environment-dev.yml -f environment-python-dev.yml
```

Note: this doesn't include a recent C++ compiler supporting the C++17 standard
(if you need to install one, you could also use conda for that).

Then activate the conda environment (needed each time you open a new terminal or
session):

```bash
$ conda activate fastscapelib-dev
```

(contributing-pre-commit)=
### Pre-Commit

Fastscapelib also provides a configuration file for setting up Git hooks via
[pre-commit](https://pre-commit.com/). This will automatically check and format
the C++ and Python code before committing, using
[Clang-Format](https://clang.llvm.org/docs/ClangFormat.html) and
[Black](https://black.readthedocs.io) respectively.

``pre-commit`` has to be installed separately in your activated conda
environment:

```bash
$ conda install pre-commit -c conda-forge
```

Then run the following command once to configure the Git hooks:

```bash
$ pre-commit install
```

## Contributing to the Code Base

Once the {ref}`development environment <contributing-development-environment>`
is created and activated, each time after editing the code you might want to
perform one or several of the following tasks depending on what has been added,
removed or modified:

- build and run the {ref}`C++ tests <run-cpp-tests>`
- {ref}`(re)-install the Python library <install-python-source>` and run the
  {ref}`Python tests <run-python-tests>`
- build and run the {ref}`benchmarks <run-benchmarks>`
- check Python type hints added in ``.pyi`` files using
  [mypy](https://mypy-lang.org/)

To run ``mypy``:

```bash
$ python -m mypy .
```

## Contributing to the Documentation

The documentation is written in various formats:

- [reStructuredText
  (reST)](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html):
  API reference and docstrings in the code,
- [Markedly Structured Text (MyST)](https://mystmd.org): user and developer guides
- [Jupyter Notebook (ipynb)](https://jupyter.org/): examples

It is built using [Sphinx](https://www.sphinx-doc.org),
[MyST-Parser](https://myst-parser.readthedocs.io),
[MyST-NB](https://myst-nb.readthedocs.io) as well as a few other Sphinx
extensions.

The API reference documentation is generated automatically from the docstrings
found in the code using [Doxygen](https://www.doxygen.nl/index.html),
[Breathe](https://www.breathe-doc.org/) and the Sphinx [autodoc
extension](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html).

The docstrings follow the [Doxygen style
conventions](https://www.doxygen.nl/manual/docblocks.html) for the C++ API and
the [Numpydoc style
guide](https://numpydoc.readthedocs.io/en/latest/format.html) for the Python
API.

You can install all the documentation development tools in a new conda
environment with the following command:

```bash
$ conda env create -f doc/environment.yml
$ conda activate fastscapelib-docs
```

Fastscapelib provides a [Makefile](https://makefiletutorial.com/) to help
building the documentation locally with the following commands:

```bash
$ cd doc
$ make html
```

This will run ``doxygen`` to generate xml files from the C++ docstrings and then
``sphinx-build`` to generate API source files, run the notebook examples and
render all the source files as html.

If the build completed successfully, you can then open `_build/html/index.html`
with your web-browser to visualize the output files.

Building the documentation for the first time may take a little while. During
the next builds, Sphinx will only rebuild the modified source files. To trigger
a full build from scratch, run:

```bash
$ make clean
$ make html
```

## Contributing Your Changes to Fastscapelib

Once you are ready to submit your fixes or changes to the Fastscapelib
repository, you can add and commit it using Git:

```bash
$ git add <files_or_paths_to_add>
$ git commit -m "replace this by a meaningful commit message"
```

:::{tip}

Run the ``git status`` command to check which files have been changed.

:::

If you configured {ref}`pre-commit <contributing-pre-commit>`, the last command
will automatically check and reformat the code (after installing the
auto-formatting tools if needed). In case some code is reformatted, you'll need
to repeat the two commands above.

At any time you can push the last commits of your local branch to your fork on
GitHub with:

```bash
$ git push origin your-branch-name
```

Then, if you visit again the Fastscapelib repository (or your fork) on GitHub
you should see an invite to create a new pull-request from your uploaded branch.

It is possible that new commits are automatically added in the pull-request by
the pre-commit {ref}`continuous integration <contributing-ci>` service to
reformat some files. It that's the case, don't forget to pull these commits to
keep your local branch synchronized:

```bash
$ git pull origin your-branch-name
```

If you need to make further changes, you can repeat the add / commit / push Git
operations. The pull-request will be automatically updated with the new changes.

(contributing-ci)=
## Continuous Integration

After each commit in a pull-request or after merging a pull request, all the
following tasks are run automatically using GitHub Actions:

- build and run the C++ tests on Linux, Mac and Windows for different
  compilers and/or versions
- build the Python package and run the Python tests on Linux, Mac and Windows
  for different Python versions
- run pre-commit to check for code formatting inconsistencies
- run mypy for checking Python type hints
- build the documentation on readthedocs
