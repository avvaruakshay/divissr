#1 /usr/bin/env python

import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'divissr',
    version = '0.1.0',
    author = 'Akshay Avvaru',
    author_email = 'avvaru@ccmb.res.in',
    description = 'DiviSSR is a DNA tandem repeat identification tool',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/avvaruakshay/divissr.git',
    project_urls = {
        'Bug Tracker': 'https://github.com/avvaruakshay/divissr/issues'
    },
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    packages = setuptools.find_packages(),
    entry_points={
        'console_scripts': ['divissr=divissr.core:main']
    },
    include_package_data=True, # change path according to package name in MANIFEST.in
)