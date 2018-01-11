egpcrunch
=========

Number crunching module, used in cosmology research.


Installation
------------

**On Unix (Linux, macOS)**

 - clone this repository
 - `pip install ./egpcrunch`

On macOS, if you installed gcc via MacPorts, you can get compiler issues. In that case, specify the Apple clang:

 - `CXX=cc CC=cc pip install ./egpcrunch`

**On Windows (Requires Visual Studio 2015)**

 - For Python 3.5:
     - clone this repository
     - `pip install ./egpcrunch`
 - For earlier versions of Python, including Python 2.7:

   xtensor requires a C++14 compliant compiler (i.e. Visual Studio 2015 on
   Windows). Running a regular `pip install` command will detect the version
   of the compiler used to build Python and attempt to build the extension
   with it. We must force the use of Visual Studio 2015.

     - clone this repository
     - `"%VS140COMNTOOLS%\..\..\VC\vcvarsall.bat" x64`
     - `set DISTUTILS_USE_SDK=1`
     - `set MSSdk=1`
     - `pip install ./egpcrunch`

   Note that this requires the user building `egpcrunch` to have registry edition
   rights on the machine, to be able to run the `vcvarsall.bat` script.


Windows runtime requirements
----------------------------

On Windows, the Visual C++ 2015 redistributable packages are a runtime
requirement for this project. It can be found [here](https://www.microsoft.com/en-us/download/details.aspx?id=48145).

If you use the Anaconda python distribution, you may require the Visual Studio
runtime as a platform-dependent runtime requirement for you package:

```yaml
requirements:
  build:
    - python
    - setuptools
    - pybind11

  run:
   - python
   - vs2015_runtime  # [win]
```


Building the documentation
--------------------------

Documentation for the example project is generated using Sphinx. Sphinx has the
ability to automatically inspect the signatures and documentation strings in
the extension module to generate beautiful documentation in a variety formats.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

 - `egpcrunch/docs`
 - `make html`


Running the tests
-----------------

Running the tests requires `pytest`.

```bash
py.test .
```


Credits
-------

All the "boilerplate" in this package was created in a minute or two using the [xtensor-python cookiecutter template](https://github.com/QuantStack/xtensor-python-cookiecutter).