from setuptools import setup, Extension

#  https://stackoverflow.com/questions/4529555/building-a-ctypes-based-c-library-with-distutils
from distutils.command.build_ext import build_ext as build_ext_orig
class CTypesExtension(Extension): pass
class build_ext(build_ext_orig):

    def build_extension(self, ext):
        self._ctypes = isinstance(ext, CTypesExtension)
        return super().build_extension(ext)

    def get_export_symbols(self, ext):
        if self._ctypes:
            return ext.export_symbols
        return super().get_export_symbols(ext)

    def get_ext_filename(self, ext_name):
        if self._ctypes:
            return ext_name + '.so'
        return super().get_ext_filename(ext_name)


setup(
    name="p3p",
    version="1.0.2",
    py_modules = ["p3p.p3p"],
    ext_modules=[ CTypesExtension( "p3p.p3plib", ["c/p3p.c"],),],
    cmdclass={'build_ext': build_ext},
)
