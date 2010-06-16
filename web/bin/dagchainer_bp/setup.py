from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


version = '0.1'


setup(name='dagtools',
      version=version,
      description="",
      long_description="""\
""",
      cmdclass= {'build_ext': build_ext },
      ext_modules=[ Extension("cdagline",
                      sources=["dagtools/cdagline.c"],)],
      keywords='bio',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='BSD',
      test_suite='nose.collector',
      #package_dir = {'': 'dagtools'},
      zip_safe=False,
      packages = ['dagtools'],
      entry_points={
                 'console_scripts': ['dagtools = dagtools:main']
      },


  )
