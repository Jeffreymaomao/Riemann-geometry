from setuptools import setup, find_packages

setup(
    name='riemann_geometry',
    version='0.1.0',
    packages=find_packages(),
    description='A package for calculating Riemannian geometry tensors',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='ChangMao Yang',
    author_email='jeffrey0613mao@gmail.com',
    install_requires=[
        'sympy',
        'tqdm'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
