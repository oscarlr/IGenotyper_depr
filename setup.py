from setuptools import setup, find_packages

setup(
    name='IGenotyper_clean',
    description='',
    packages=find_packages(),
    package_data={'IGenotyper': ['data/*','bash_scripts/*','python_scripts/*','r_scripts/*']},
    include_package_data=True,
    entry_points = {
        'console_scripts': ['IG-clean = IGenotyper.main:main',
                            'IG-make-ref = IGenotyper.make_ref:main',
                            'IG-phase-reads =  IGenotyper.python_scripts.phase_reads:main',
                            'IG-filter-vcf =  IGenotyper.python_scripts.filter_vcf:main']
        },
    platforms='any'
)
