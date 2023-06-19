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

(contributing-development-environment)=
## Creating a Development Environment

Follow the steps below to setup a new environment for Fastscapelib development.

:::{note}

Contributing to the code or documentation requires some basic knowledge of
[Git](https://git-scm.com/).

We also recommend installing all the dependencies and development tools in a new
[conda](https://conda.io/docs/) environment created specifically for this
project.

:::

1. Go on the Fastscapelib
   [GitHub](https://github.com/fastscape-lem/fastscapelib) page and click on the
   "Fork" button (you need to have an account on GitHub and be logged in).

2. Clone your Fastscapelib forked Git repository:

```
$ git clone https://github.com/your-username/fastscapelib
$ cd fastscapelib
```

3. Install Fastscapelib's required dependencies, for example in a new ``conda``
   environment. See the relevant installation sections ({ref}`C++
   <install-cpp-source>` and {ref}`Python <install-python-source>`) for more
   details.

4. Fastscapelib provides a configuration file for setting up Git hooks via
   [pre-commit](https://pre-commit.com/). This will automatically format the C++
   and Python code (using
   [Clang-Format](https://clang.llvm.org/docs/ClangFormat.html) and
   [Black](https://black.readthedocs.io), respectively) before creating a new
   commit. To install and configure ``pre-commit``, run the commands below (be
   sure your ``conda`` environment is activated):

```
$ conda install pre-commit -c conda-forge
$ pre-commit install
```

5. Optionally install [mypy](https://mypy-lang.org/) for checking the Python
   type hints:

```
$ conda install mypy -c conda-forge
```

6. Create a new Git branch in which you will implement your bug fixes,
   documentation enhancement or new feature:

```
$ git checkout -b your-branch-name
```

7. You are all set!

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
  [mypy](https://mypy-lang.org/).

To run ``mypy``:

```
$ cd python
$ python -m mypy . --exclude=build
```

## Contributing to the Documentation

The documentation is written in the [reStructuredText
(reST)](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html),
[Markedly Structured Text (MyST)](https://mystmd.org) and Jupyter Notebook
(ipynb) formats. It is built using [Sphinx](https://www.sphinx-doc.org),
[MyST-Parser](https://myst-parser.readthedocs.io),
[MyST-NB](https://myst-nb.readthedocs.io) as well as a few other Sphinx
extensions.

The API reference documentation is generated automatically from the docstrings
found in the code using [Doxygen](https://www.doxygen.nl/index.html),
[Breathe](https://www.breathe-doc.org/) and the Sphinx [autodoc
extension](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html).

You can install all the documentation development tools in a new conda
environment with the following command:

```
$ conda env create -f doc/environment.yml
$ conda activate fastscapelib-docs
```

Fastscapelib provides a [Makefile](https://makefiletutorial.com/) to help
building the documentation locally with the following commands:

```
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

```
$ make clean
$ make html
```

## Contributing Your Changes to Fastscapelib

Once you are ready to submit your fixes or changes to the Fastscapelib
repository, you can add / commit / push it using Git:

```
$ git add <files_or_paths_to_add>
$ git commit -m "replace this by a meaningful commit message"
$ git push origin your-branch-name
```

Then, if you visit again the Fastscapelib repository (or your fork) on GitHub
you should see an invite to create a new pull-request from your uploaded branch.

If you need to make further changes, you can repeat the add / commit / push Git
operations. The pull-request will be automatically updated with the new changes.

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
