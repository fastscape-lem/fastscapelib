import os
import platform
import subprocess
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

CMAKE_OPTIONS = {"BUILD_PYTHON_MODULE": "ON"}


class CMakeExtension(Extension):
    def __init__(self, name, source_dir="", cmake_options=None):
        Extension.__init__(self, name, sources=[])

        self.source_dir = os.path.abspath(source_dir)
        self.cmake_options = cmake_options or {}


class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]
        cmake_args += ["-D{}={}".format(k, v) for k, v in ext.cmake_options.items()]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]

        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # note: need to run "pip install . -v" to see cmake output
        subprocess.check_call(
            ["cmake", ext.source_dir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


ext_modules = [
    CMakeExtension("_fastscapelib_py", source_dir="..", cmake_options=CMAKE_OPTIONS)
]


setup(
    name="fastscapelib",
    version="0.1.3",
    maintainer="Benoit Bovy",
    maintainer_email="benbovy@gmail.com",
    url="https://github.com/fastscape-lem/fastscapelib",
    description=(
        "A library of efficient algorithms for processing "
        "topographic data and landscape evolution modeling"
    ),
    long_description="",  # TODO: import README here
    packages=find_packages(),
    package_data={"fastscapelib": ["**/*.pyi"]},
    ext_modules=ext_modules,
    install_requires=["numpy"],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
