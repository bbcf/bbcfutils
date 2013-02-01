from distutils.core import setup

setup(
        name             =   'maplot',
        version          =   '1.0.3',
        description      =   'Python package for creating interactive MA-plots',
        long_description =   open('README').read(),
        license          =   'GNU General Public License 3.0',
        url              =   'http://delafont.github.com/maplot',
        author           =   'Julien Delafontaine',
        author_email     =   'julien.delafontaine@yandex.com',
        classifiers      =   ['Topic :: Scientific/Engineering :: Bio-Informatics'],
        packages         =   ['maplot',],
        install_requires =   ['scipy','numpy','matplotlib'],
        scripts          =   ['maplot/maplot'],
    )
