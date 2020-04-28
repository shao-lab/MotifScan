"""
motifscan.exceptions
--------------------

This module contains the exceptions of MotifScan.
"""


class MotifScanError(Exception):
    """Base class for MotifScan errors."""


class InvalidConfigFileError(MotifScanError):
    def __init__(self, path):
        msg = f"Invalid config file: {path}"
        super(MotifScanError, self).__init__(msg)


class RemoteGenomeNotFoundError(MotifScanError):
    def __init__(self, database, assembly):
        msg = f"No genome assembly {assembly!r} in the {database} database"
        super(MotifScanError, self).__init__(msg)


class RemoteGenomeFileNotFoundError(MotifScanError):
    def __init__(self, database, assembly, which):
        msg = f"No {which} file for {assembly!r} in the {database} database"
        super(MotifScanError, self).__init__(msg)


class GenomeNotFoundError(MotifScanError):
    def __init__(self, name):
        msg = f"No such genome assembly: {name!r}"
        super(MotifScanError, self).__init__(msg)


class GenomeFileNotFoundError(MotifScanError):
    def __init__(self, name, which):
        msg = f"No {which} file for assembly {name!r}"
        super(MotifScanError, self).__init__(msg)


class BackgroundFormatError(MotifScanError):
    def __init__(self, line_num, line):
        msg = f"Invalid background format at line {line_num}: {line!r}"
        super(MotifScanError, self).__init__(msg)


class RemoteMotifPFMsNotFoundError(MotifScanError):
    def __init__(self, database, pfms):
        msg = f"No motif PFMs {pfms!r} in the {database} database"
        super(MotifScanError, self).__init__(msg)


class MotifSetNotFoundError(MotifScanError):
    def __init__(self, name):
        msg = f"No such motif set: {name!r}"
        super(MotifScanError, self).__init__(msg)


class PfmsFileNotFoundError(MotifScanError):
    def __init__(self, name):
        msg = f"No PFMs file for motif set {name!r}"
        super(MotifScanError, self).__init__(msg)


class PwmsFileNotFoundError(MotifScanError):
    def __init__(self, name, genome):
        msg = f"No PWMs file for motif set {name!r} under genome {genome}"
        super(MotifScanError, self).__init__(msg)


class PfmsJasparFormatError(MotifScanError):
    def __init__(self, line_num, line):
        msg = f"Invalid JASPAR PFMs format at line {line_num}: {line!r}"
        super(MotifScanError, self).__init__(msg)


class PwmsMotifScanFormatError(MotifScanError):
    def __init__(self, line_num, line):
        msg = f"Invalid MotifScan PWMs format at line {line_num}: {line!r}"
        super(MotifScanError, self).__init__(msg)


class RegionFileFormatError(MotifScanError):
    def __init__(self, format, line_num, line):
        msg = f"Invalid {format} format at line {line_num}: {line!r}"
        super().__init__(msg)
