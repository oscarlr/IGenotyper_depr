from setuptools import setup, find_packages

setup(
    name='IGenotyper_clean',
    description='',
    packages=find_packages(),
    package_data={'IGenotyper': ['data/*']},
    include_package_data=True,
    entry_points = {
        'console_scripts': ['IG = IGenotyper.IGenotyper:main'],
        'console_scripts': ['IG-make-ref = IGenotyper.make_ref:main'],
        },
    platforms='any'
)
