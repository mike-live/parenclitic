from setuptools import setup

setup(name='parenclitic',
      version='0.1.1',
      description='Parenclitic approach with kernels inside',
      url='https://github.com/mike-live/parenclitic',
      author='Mikhail Krivonosov',
      author_email='mike_live@mail.ru',
      license='MIT',
      packages=['parenclitic'],
      install_requires=[
          'numpy',
          'python-igraph',
          'pandas',
          'sklearn',
          'scipy',
          'multiprocessing'
      ],
      zip_safe=False)