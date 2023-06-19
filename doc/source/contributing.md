# Contributing to Fastscapelib

Fastscapelib is open source and welcomes your feedback and/or contributions!
There are many ways to contribute to the development of Fastscapelib, including
bug reports, bug fixes, documentation improvements, enhancement suggestions, and
other ideas.

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

:::

:::{tip}

GitHub supports the [markdown
syntax](https://docs.github.com/en/get-started/writing-on-github) for formatting
code, add links, etc. This is helpful for writing clear reports.

:::

## Creating a Development Environment

Contributing to the code or documentation requires some basic knowledge of
version control systems and more specifically [Git](https://git-scm.com/).

Follow the steps below to setup a new environment for Fastscapelib development.

1. Clone the Fastscapelib git repository:

```
$ git clone https://github.com/fastscape-lem/fastscapelib
```

2. Create a new branch in which you will implement your bug fixes, documentation
   enhancement or new feature:

```
$ git checkout -b my-bug-fix-or-feature
```

3. Install Fastscapelib's required dependencies. We recommend creating a new
   [conda](https://conda.io/docs/) environment. See the relevant Sections
   {ref}`install-cpp-source` and {ref}`install-python-source` for more details.

## Contributing to the Code Base

After editing the code, you might want to build and run the {ref}`C++ tests <run-cpp-tests>`, {ref}`(re)-install the Python library <install-python-source>` and run the {ref}`Python tests <run-python-tests>` and/or build and run the {ref}`benchmarks <run-benchmarks>`.

## Contributing to the Documentation

The documentation is written in the [reStructuredText
(reST)](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html),
[Markedly Structured Text (MyST)](https://mystmd.org) and Jupyter Notebook
(ipynb) formats. It is built using [Sphinx](https://www.sphinx-doc.org),
[MyST-Parser](https://myst-parser.readthedocs.io),
[MyST-NB](https://myst-nb.readthedocs.io) as well as a few other Sphinx
extensions.

The API reference documentation is generated automatically from the docstrings found in the code using [Doxygen](https://www.doxygen.nl/index.html), [Breathe](https://www.breathe-doc.org/) and the Sphinx [autodoc extension](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html),

You can install all the documentation development tools in a new conda environment using:

```
$ conda env create -f doc/environment.yml
$ conda activate fastscapelib-docs
```
