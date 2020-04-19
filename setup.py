import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='razi',
    version='0.0.0',
    description='Using SQLAlchemy with chemical databases',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Environment :: Plugins",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    url='http://razi.readthedocs.org/',
    keywords='chemistry cheminformatics sqlalchemy orm',
    packages=setuptools.find_packages(exclude=[
        'docs', 'docs.*',
        'tests', 'tests.*',
    ]),
    install_requires=[
        'SQLAlchemy>=0.7.0',
    ],
    python_requires='>=3.6',
)
