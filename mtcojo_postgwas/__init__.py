"""
mtcojo_postgwas package.

Public functions are loaded lazily so lightweight commands such as
`mtcojo-postgwas --help` do not import plotting, LDSC, or VCF dependencies.
"""

__version__ = "1.0.0"
__author__ = "JJOHN41"

_EXPORTS = {
    "convert_vcf_single_pass": ("mtcojo_postgwas.io.vcf_converter", "convert_vcf_single_pass"),
    "sanitize_bim": ("mtcojo_postgwas.io.bim_sanitizer", "sanitize_bim"),
    "detect_bim_id_format": ("mtcojo_postgwas.io.bim_sanitizer", "detect_bim_id_format"),
    "run_gcta_mtcojo": ("mtcojo_postgwas.stages.gcta", "run_gcta_mtcojo"),
    "run_postgwas_harmonisation": ("mtcojo_postgwas.stages.postgwas", "run_postgwas_harmonisation"),
    "run_ldsc_pipeline": ("mtcojo_postgwas.stages.ldsc", "run_ldsc_pipeline"),
    "get_logger": ("mtcojo_postgwas.core.logger", "get_logger"),
}

__all__ = ["__version__", "__author__", *_EXPORTS]


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    from importlib import import_module

    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value
