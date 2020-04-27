"""
motifscan.config
----------------

MotifScan configuration handler, this configuration manages the paths of
genome and motif data files in the file system.
"""

import os
from configparser import ConfigParser, Error

from motifscan.exceptions import InvalidConfigFileError, GenomeNotFoundError, \
    MotifSetNotFoundError

user_rc_path = os.path.expanduser("~/.motifscanrc")
user_genome_dir = os.path.expanduser("~/.motifscan/genomes/")
user_motif_dir = os.path.expanduser("~/.motifscan/motifs/")


class Config:
    """Configuration handler for MotifScan."""

    _sections = ['motifscan', 'genome', 'motif']

    def __init__(self, path=None):
        self.path = path or user_rc_path
        self._config = ConfigParser(allow_no_value=False)
        try:
            self._config.read(self.path)
        except Error as e:
            raise InvalidConfigFileError(self.path) from e
        # set default for all sections
        for section in self._sections:
            if not self._config.has_section(section):
                self._config.add_section(section)
        # set the default genome root directory
        if not self._config.has_option('motifscan', 'genome_dir'):
            self.set_genome_dir(user_genome_dir)
        # set the default motif root directory
        if not self._config.has_option('motifscan', 'motif_dir'):
            self.set_motif_dir(user_motif_dir)

    def get_genome_dir(self):
        """Get the genome root directory."""
        return self._config.get('motifscan', 'genome_dir')

    def set_genome_dir(self, path):
        """Set the specified path as the genome root directory."""
        self._config.set('motifscan', 'genome_dir', path)

    def get_motif_dir(self):
        """Get the motif root directory."""
        return self._config.get('motifscan', 'motif_dir')

    def set_motif_dir(self, path):
        """Set the specified path as the motif root directory."""
        self._config.set('motifscan', 'motif_dir', path)

    def list_genome_assemblies(self):
        """List configured (installed) genome assemblies."""
        for name, path in self._config.items('genome'):
            yield name, path

    def has_genome_assembly(self, name):
        """Returns if the specified genome assembly is configured."""
        return self._config.has_option('genome', name)

    def get_genome_path(self, name):
        """Get the genome path of the specified genome assembly."""
        if self._config.has_option('genome', name):
            return self._config.get('genome', name)
        else:
            raise GenomeNotFoundError(name)

    def set_genome_path(self, name, path):
        """Set the genome path for the specified genome assembly."""
        self._config.set('genome', name, path)

    def remove_genome_path(self, name):
        """Remove the specified genome assembly out of the configuration."""
        if self._config.has_option('genome', name):
            return self._config.remove_option('genome', name)
        else:
            raise GenomeNotFoundError(name)

    def list_motif_sets(self):
        """List configured (installed) motif PFMs sets."""
        for name, path in self._config.items('motif'):
            yield name, path

    def has_motif_set(self, name):
        """Returns if the specified motif PFMs set is configured."""
        return self._config.has_option('motif', name)

    def get_motif_path(self, name):
        """Get the motif path of the specified motif PFMs set."""
        if self._config.has_option('motif', name):
            return self._config.get('motif', name)
        else:
            raise MotifSetNotFoundError(name)

    def set_motif_path(self, name, path):
        """Set the motif path for the specified motif PFMs set."""
        self._config.set('motif', name, path)

    def remove_motif_path(self, name):
        """Remove the specified motif PFMs set out of the configuration."""
        if self._config.has_option('motif', name):
            return self._config.remove_option('motif', name)
        else:
            raise MotifSetNotFoundError(name)

    def write(self, path=None):
        """Save the configuration."""
        path = path or self.path
        with open(path, 'w') as f_config:
            self._config.write(f_config)
