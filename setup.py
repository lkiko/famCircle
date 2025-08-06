# -*- coding: UTF-8 -*-

from setuptools import find_packages, setup

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

required = ['pandas', 'numpy', 'biopython', 'matplotlib', 'scipy','interval','seaborn','tqdm']

setup(
    name="famCircle",
    version="0.2.6",
    author="charles_kiko",
    author_email="charles_kiko@163.com",
    description="famCircle",
    license="BSD License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    package_data={'': ['*.conf','*.ini', '*.csv','*.ttf']},
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'famCircle = famCircle.run:main',
        ]
    },
    zip_safe=True,
    install_requires=required
)