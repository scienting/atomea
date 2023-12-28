def apply_dtype_map(fields, dtype_map):
    new_fields = []
    for field in fields:
        new_fields.append((field[0], dtype_map[field[1]]))
    return new_fields
