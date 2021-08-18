"""
seekr2
Simulation-Enabled Estimation of Kinetic Rates - Version 2
"""
import sys
from setuptools import setup, find_packages
import versioneer
import subprocess
import fileinput
import re

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='seekr2',
    author='Lane Votapka',
    author_email='lvotapka100@gmail.com',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    install_requires=["numpy", "scipy", "parmed", "pytest", "matplotlib", 
                       "nptyping", "mdtraj", "bubblebuster"],       
    platforms=['Linux',
    #            'Mac OS-X',
                'Unix',]
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.7",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
if "install" in sys.argv:
    try:
        import abserdes
        
    except ImportError:
        abserdes_repo_url = 'https://github.com/astokely/abserdes.git'
        process = subprocess.Popen([
            "git",
            "ls-remote",
            'https://github.com/astokely/abserdes.git'
        ], stdout=subprocess.PIPE)
        stdout = process.communicate()[0]
        sha_tarball = re.split(r'\t+', stdout.decode('ascii'))[0] + ".tar.gz"
        abserdes_tarball_link = [line for line in open("requirements.txt" , "r+") if 'abserdes' in line][0]
        for line in fileinput.input("requirements.txt", inplace = 1):
            print(line.replace(
                abserdes_tarball_link.rsplit('/', 1)[-1],
                sha_tarball
            ))
        
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])