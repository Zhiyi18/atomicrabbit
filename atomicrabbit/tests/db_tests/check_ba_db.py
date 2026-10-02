import sqlite3

conn = sqlite3.connect("atomicrabbit/db/spectroscopy.db")
conn.row_factory = sqlite3.Row

print("=== Ba atoms in DB ===")
rows = conn.execute("SELECT atom_id, element, mass_number, nuclear_spin FROM atoms WHERE element='Ba' ORDER BY mass_number").fetchall()
for r in rows:
    print(dict(r))

print()
print("=== Ba energy levels in DB ===")
rows = conn.execute("""
    SELECT a.mass_number, el.configuration, el.term_symbol, el.J, el.source_type, el.is_recommended
    FROM energy_levels el JOIN atoms a ON el.atom_id=a.atom_id
    WHERE a.element='Ba' ORDER BY a.mass_number, el.energy_cm1
""").fetchall()
for r in rows:
    print(dict(r))

conn.close()
