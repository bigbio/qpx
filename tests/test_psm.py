from pathlib import Path
import json

from quantmsio.core.psm import Psm

TEST_DATA_ROOT = Path(__file__).parent / "examples"


def test_convert_mztab_to_feature():
    mztab_path = TEST_DATA_ROOT / "DDA-lfq/PXD040438.mzTab"
    psm = Psm(mztab_path)
    for _ in psm.generate_report():
        print("ok")


def test_psm_schema_contains_new_array_fields():
    """Test that the PSM schema contains the new fragment annotation array fields."""
    schema_path = Path(__file__).parent.parent / "docs" / "psm.avsc"
    
    assert schema_path.exists(), f"Schema file not found: {schema_path}"
    
    with open(schema_path, 'r') as f:
        schema = json.load(f)
    
    # Find the PSM record fields
    psm_fields = None
    for field in schema.get("fields", []):
        if field.get("name") == "psms":
            psm_type = field.get("type", {})
            if isinstance(psm_type, dict) and psm_type.get("type") == "array":
                items = psm_type.get("items", {})
                if isinstance(items, dict) and items.get("name") == "psm":
                    psm_fields = items.get("fields", [])
                    break
    
    assert psm_fields is not None, "Could not find PSM fields in schema"
    
    # Check for new fields
    field_names = [field.get("name") for field in psm_fields]
    
    # Verify the new fragment annotation fields are present
    assert "charge_array" in field_names, "charge_array field missing from PSM schema"
    assert "ion_type_array" in field_names, "ion_type_array field missing from PSM schema"
    assert "ion_mobility_array" in field_names, "ion_mobility_array field missing from PSM schema"
    
    # Verify field types are correct
    for field in psm_fields:
        field_name = field.get("name")
        if field_name == "charge_array":
            field_type = field.get("type")
            assert field_type == ["null", {"type": "array", "items": "int32"}], f"charge_array has incorrect type: {field_type}"
        elif field_name == "ion_type_array":
            field_type = field.get("type")
            assert field_type == ["null", {"type": "array", "items": "string"}], f"ion_type_array has incorrect type: {field_type}"
        elif field_name == "ion_mobility_array":
            field_type = field.get("type")
            assert field_type == ["null", {"type": "array", "items": "float32"}], f"ion_mobility_array has incorrect type: {field_type}"
