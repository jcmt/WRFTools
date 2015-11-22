'''
wrftools ...
'''

classifiers = '''
              Development Status :: alpha
              Environment :: Console
              Intended Audience :: Science/Research
              Intended Audience :: Developers
              License :: GNU
              Operating System :: OS Independent
              Programming Language :: Python
              Topic :: Scientific/Engineering
              Topic :: Software Development :: Libraries :: Python Modules
              '''

from numpy.distutils.core import Extension
from numpy.distutils.command.install import install
from glob import glob

class my_install(install):
    def run(self):
        install.run(self)

        print '''
        enjoy wrftools
        '''


doclines = __doc__.split("\n")

if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(name='wrftools',
          version=0,
          description=doclines[0],
          long_description="\n".join(doclines[2:]),
          author='Joao Teixeira',
          author_email='jcmt87@gmail.com',
          url='NA',
          packages=['wrftools'],
          license='GNU',
          platforms=['any'],
          ext_modules=[],
          data_files=[('wrftools/cmap', glob('wrftools/cmap/*')),
                      ('wrftools/colormaps', glob('wrftools/colormaps/*')),
                      ('wrftools/', glob('wrftools/*.so'))],
          classifiers = filter(None, classifiers.split("\n")),
          cmdclass={'install': my_install},
          )

