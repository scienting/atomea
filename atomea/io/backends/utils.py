def apply_dtype_map(fields, dtype_map):
    new_fields = []
    for key, dtype in fields.items():
        new_fields.append((key, dtype_map[dtype]))
    return new_fields
