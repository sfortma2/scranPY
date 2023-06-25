from setuptools import setup

setup(
    name='scranPY',
    version='0.1',
    description='A python implementation of scran computeSumFactors: normalization by deconvolution for single-cell RNA-sequencing',
    author='Seth Fortmann',
    author_email='seth.fortmann@gmail.com',
    url='https://github.com/sfortma2/scranPY',
    packages=['scranPY'],
    install_requires=[
        # List any dependencies your module requires
        'scanpy',
        'anndata',
        'numpy',
        'pandas',
        'multiprocessing',
        'scipy',
        'cvxpy',
        'matplotlib',
    ],
)
