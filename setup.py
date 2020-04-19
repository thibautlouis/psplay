# -*- coding: utf-8 -*-

import json
import os
import sys
from distutils import log
from glob import glob
from subprocess import check_call

from setuptools import Command, find_packages, setup
from setuptools.command.build_py import build_py
from setuptools.command.egg_info import egg_info
from setuptools.command.sdist import sdist

here = os.path.dirname(os.path.abspath(__file__))
node_root = os.path.join(here, 'js')
is_repo = os.path.exists(os.path.join(here, '.git'))

npm_path = os.pathsep.join([
    os.path.join(node_root, 'node_modules', '.bin'),
    os.environ.get('PATH', os.defpath),
])

log.set_verbosity(log.DEBUG)
log.info('setup.py entered')
log.info('$PATH=%s' % os.environ['PATH'])

def js_prerelease(command, strict=False):
    """Decorator for building minified js/css prior to another command."""
    class DecoratedCommand(command):
        def run(self):
            jsdeps = self.distribution.get_command_obj('jsdeps')
            if not is_repo and all(os.path.exists(t) for t in jsdeps.targets):
                # sdist, nothing to do
                command.run(self)
                return

            try:
                self.distribution.run_command('jsdeps')
            except Exception as e:
                missing = [t for t in jsdeps.targets if not os.path.exists(t)]
                if strict or missing:
                    log.warn('rebuilding js and css failed')
                    if missing:
                        log.error('missing files: %s' % missing)
                    raise e
                else:
                    log.warn('rebuilding js and css failed (not a problem)')
                    log.warn(str(e))
            command.run(self)
            update_package_data(self.distribution)
    return DecoratedCommand


def get_data_files():
    with open(os.path.join(node_root, 'package.json')) as f:
        package_json = json.load(f)
    tgz = '%s-%s.tgz' % (package_json['name'], package_json['version'])

    return [
        ('share/jupyter/nbextensions/jupyter-leaflet-car', glob('psplay/static/*')),
        ('etc/jupyter/nbconfig/notebook.d', ['jupyter-leaflet-car.json']),
        ('share/jupyter/lab/extensions', ['js/' + tgz]),
    ]


def update_package_data(distribution):
    """Update package_data to catch changes during setup."""
    build_py = distribution.get_command_obj('build_py')
    distribution.data_files = get_data_files()
    # re-init build_py options which load package_data
    build_py.finalize_options()


class NPM(Command):
    description = 'install package.json dependencies using npm'

    user_options = []

    node_modules = os.path.join(node_root, 'node_modules')

    targets = [
        os.path.join(here, 'psplay', 'static', 'extension.js'),
        os.path.join(here, 'psplay', 'static', 'index.js')
    ]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def has_npm(self):
        try:
            check_call(['npm', '--version'])
            return True
        except:
            return False

    def run(self):
        has_npm = self.has_npm()
        if not has_npm:
            log.error("`npm` unavailable.  If you're running this command using sudo, make sure `npm` is available to sudo")

        env = os.environ.copy()
        env['PATH'] = npm_path

        if self.has_npm():
            log.info("Installing build dependencies with npm.  This may take a while...")
            check_call(['npm', 'install'], cwd=node_root, stdout=sys.stdout, stderr=sys.stderr)
            check_call(['npm', 'pack'], cwd=node_root, stdout=sys.stdout, stderr=sys.stderr)

        for t in self.targets:
            if not os.path.exists(t):
                msg = 'Missing file: %s' % t
                if not has_npm:
                    msg += '\nnpm is required to build a development version of widgetsnbextension'
                raise ValueError(msg)

        # update package data in case this created new files
        update_package_data(self.distribution)


setup_args = {
    'name': 'psplay',
    'version': '0.0.1',
    'description': 'A Jupyter widget for playing with CAR Leaflet maps',
    'license': 'BSD License',
    'include_package_data': True,
    'data_files': get_data_files(),
    'install_requires': [
        'ipywidgets>=7.5.0,<8',
        'traittypes>=0.2.1,<3',
        'ipyleaflet==0.12.3',
        'psply>=1.3.0'
    ],
    'packages': find_packages(),
    'zip_safe': False,
    'cmdclass': {
        'build_py': js_prerelease(build_py),
        'egg_info': js_prerelease(egg_info),
        'sdist': js_prerelease(sdist, strict=True),
        'jsdeps': NPM,
    },
    'author': 'Simons Observatory Collaboration Power Spectrum Task Force',
    'url': 'https://github.com/simonsobs/psplay',
    'keywords': ['CMB', 'spectra', 'ipython', 'jupyter', 'widgets', 'graphics'],
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Multimedia :: Graphics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
    ],
}

setup(**setup_args)
