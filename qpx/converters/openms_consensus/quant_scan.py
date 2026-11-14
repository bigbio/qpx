"""The scan and isolation purity an isobaric consensus feature was quantified from."""

from qpx.core.scan import scan_from_native_id


def quant_scan_fields(cf) -> dict:
    """Purity and quantifying scan IsobaricAnalyzer records on a consensus feature; empty for other inputs.

    ``quant_precursor_purity`` describes the isolation of the quantifying scan, so it
    is preferred. ``id_scan_id`` is written beside ``scan_id`` when the two scans can differ.
    """
    fields = {}
    for key in ("quant_precursor_purity", "precursor_purity"):
        if cf.metaValueExists(key):
            fields["parent_ion_fraction"] = float(cf.getMetaValue(key))
            break
    if not cf.metaValueExists("scan_id"):
        return fields
    quant_scan = str(cf.getMetaValue("scan_id"))
    fields["quant_scan"] = scan_from_native_id(quant_scan)
    if cf.metaValueExists("id_scan_id"):
        fields["quant_ms_level"] = 2 if str(cf.getMetaValue("id_scan_id")) == quant_scan else 3
    return fields
